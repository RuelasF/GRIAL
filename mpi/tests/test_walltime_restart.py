#!/usr/bin/env python3
"""Integration regression: complete RK stop, atomic checkpoint rejection, continuation.

Uses only Python's standard library. Each invocation preserves an isolated result
folder; --output selects an empty directory for campaign evidence.
"""
import argparse
from array import array
import json
import os
from pathlib import Path
import struct
import subprocess
import tempfile


def records(path):
    data = path.read_bytes()
    result, offset = [], 0
    while offset < len(data):
        n = struct.unpack_from('=i', data, offset)[0]
        assert n >= 0, 'Test reader supports small native-endian records'
        offset += 4
        result.append(data[offset:offset+n])
        offset += n
        assert struct.unpack_from('=i', data, offset)[0] == n
        offset += 4
    assert offset == len(data)
    return result


def state(path):
    r = records(path)
    v8 = r[0].strip() in (b'GRHDCP_V8',b'GRHDCP_V9')
    t, step, payload = r[-4:-1] if v8 else r[-3:]
    return struct.unpack('=d', t)[0], struct.unpack('=i', step)[0], payload


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--exe', default='./grial')
    ap.add_argument('--mpi-one', action='store_true', help='launch the one-process MPI executable with mpirun')
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    exe = str(Path(args.exe).resolve())
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-walltime-'))
    out.mkdir(parents=True, exist_ok=True)
    assert not (out/'segmented').exists(), 'Use a new output folder'
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE')
    checks = []

    def run(name, text, *, check=False, success=True):
        work = out/name
        work.mkdir(parents=True)
        par = work/'case.par'
        par.write_text(text)
        # A singleton bootstraps its own PMIx server inside the 0.30 s budget.
        # Production MPI launches supply that server before the executable starts.
        prefix = ['mpirun', '--oversubscribe', '--bind-to', 'none', '-n', '1'] if args.mpi_one else []
        cmd = prefix + [exe] + (['--check'] if check else []) + [str(par)]
        p = subprocess.run(cmd, cwd=work, env=env, capture_output=True, timeout=90)
        (work/'run.log').write_bytes(p.stdout+p.stderr)
        assert (p.returncode == 0) == success, (name, p.returncode, (p.stdout+p.stderr)[-2500:])
        checks.append({'name': name, 'returncode': p.returncode})
        return work

    base = '''problem=sod
nx=800
reconstruction=weno5
riemann_solver=hlle
shock_sensor=true
final_time=0.2
cfl=0.3
save_interval=0.05
checkpoint_interval=10
output_folder=native
output_prefix=sod
numerical_diagnostics=true
diagnostic_stride=5
rho_floor=2e-11
p_floor=7e-15
'''
    limited = base + 'walltime_limit_seconds=0.30\nwalltime_margin_seconds=0.05\nwalltime_checkpoint_seconds=0.05\n'
    coldcheck = run('cold_check', limited, check=True)
    assert not (coldcheck/'native').exists()
    seg = run('segmented', limited)
    def status(work):
        path = work/'native/run_status.dat'
        if not path.exists():
            path = seg/'native/run_status.dat'
        return dict(line.split('=', 1) for line in path.read_text().splitlines())
    stopped = status(seg)
    assert stopped['status'] == 'stopped_for_walltime', stopped
    checkpoint = seg/stopped['checkpoint']
    t0, n0, u0 = state(checkpoint)
    assert 0 < n0 and 0 < t0 < 0.2, stopped
    assert len(list((seg/'native').glob('*.rst'))) >= 2, 'Wall-clock checkpoint did not fire'
    assert not list((seg/'native').glob('*.tmp'))
    r = records(checkpoint)
    assert struct.unpack('=4d', r[24])[:2] == (2e-11, 7e-15)
    restart = f'restart=true\nrestart_file={checkpoint}\nwalltime_limit_seconds=0\noutput_folder=../segmented/native\n'
    restartcheck = run('restart_check', restart, check=True)
    assert not (restartcheck/'native').exists()

    # Include every prior step; compare all five conserved fields at the endpoint.
    continuous = run('continuous', base)
    done = status(continuous)
    assert done['status'] == 'completed'
    reference_final = continuous/done['checkpoint']
    resumed = run('resumed', restart)
    resumed_status = status(resumed)
    assert resumed_status['status'] == 'completed'
    actual_final = resumed/resumed_status['checkpoint']
    tf, nf, expected = state(reference_final)
    tr, nr, actual = state(actual_final)
    assert (tf, nf) == (tr, nr), ((tf, nf), (tr, nr))
    expected_values, actual_values = array('d'), array('d')
    expected_values.frombytes(expected); actual_values.frombytes(actual)
    error = max(abs(a-b) for a, b in zip(expected_values, actual_values))
    assert actual == expected, f'Segmented conserved state differs by {error}'
    assert struct.unpack('=4d', records(actual_final)[24]) == struct.unpack('=4d', r[24])

    # All failed --check paths must reject missing records before evolving.
    # Initialization is included: an exhausted startup budget produces step zero.
    startup = base + 'walltime_limit_seconds=0.000002\nwalltime_margin_seconds=0.000001\n'
    startup_work = run('startup_stop', startup)
    startup_status = status(startup_work)
    assert startup_status['status'] == 'stopped_for_walltime'
    assert state(startup_work/startup_status['checkpoint'])[:2] == (0.0, 0)
    run('startup_restart_check', 'restart=true\nrestart_file='+str(startup_work/startup_status['checkpoint'])+'\n', check=True)

    # Original V7 payloads remain readable and explicitly announce legacy floors.
    legacy = out/'legacy_v7.rst'
    legacy_records = [b'GRHDCP_V7   '] + r[1:24] + r[-4:-1]
    legacy.write_bytes(b''.join(struct.pack('=i', len(record))+record+struct.pack('=i', len(record))
                                for record in legacy_records))
    legacy_work = run('legacy_v7_check', f'restart=true\nrestart_file={legacy}\n', check=True)
    assert b'legacy checkpoint lacks floors' in (legacy_work/'run.log').read_bytes()

    legacy8 = out/'legacy_v8.rst'
    legacy8_records = [b'GRHDCP_V8   '] + r[1:27] + r[-4:-1] + [r[-1].replace(b'GRIAL_END_V9',b'GRIAL_END_V8',1)]
    legacy8.write_bytes(b''.join(struct.pack('=i',len(record))+record+struct.pack('=i',len(record))
                                 for record in legacy8_records))
    run('legacy_v8_check',f'restart=true\nrestart_file={legacy8}\n',check=True)
    raw = checkpoint.read_bytes()
    variants = {
        'empty': b'',
        'header_only': raw[:100],
        'no_footer': raw[:-len(r[-1])-8],
        'partial_footer': raw[:-7],
        'partial_state': raw[:len(raw)//2],
        'unknown_version': raw.replace(r[0].strip(),b'GRHDCP_VX',1),
        'wrong_footer': raw.replace(r[-1][:12],b'GRIAL_BAD_V9',1),
    }
    for name, data in variants.items():
        invalid = out/(name+'.rst'); invalid.write_bytes(data)
        work = run('reject_'+name, f'restart=true\nrestart_file={invalid}\n', check=True, success=False)
        assert not (work/'native').exists()
    for key, value in [('rho_floor','0'),('p_floor','nan'),('walltime_limit_seconds','-1'),
                       ('walltime_checkpoint_seconds','-1'),('walltime_margin_seconds','0.4')]:
        prefix = limited if key.startswith('walltime') else base
        text = '\n'.join(line for line in prefix.splitlines() if not line.startswith(key+'='))
        run('invalid_'+key, text+'\n'+key+'='+value+'\n', check=True, success=False)
    report = {'passed': True, 'checks': checks, 'stop_step': n0, 'stop_time': t0,
              'final_step': nf, 'final_time': tf, 'max_abs_conserved_error': error,
              'same_payload_bytes': actual == expected, 'evidence': str(out)}
    (out/'VERIFIED.json').write_text(json.dumps(report, indent=2)+'\n')
    print('WALLTIME RESTART TESTS PASSED:', out)


if __name__ == '__main__':
    main()
