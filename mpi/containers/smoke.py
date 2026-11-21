#!/usr/bin/env python3
"""Check a portable CPU package with real OpenMP and MPI Sod evolutions."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
sys.path.insert(0, str(ROOT / 'tools'))
from test_observer_mpi import field_state
from parallel_vtk import collect_frames, read_frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--serial', type=Path, default=shutil.which('grial'))
    parser.add_argument('--mpi', type=Path, default=shutil.which('grial_mpi'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if not args.serial or not args.mpi:
        parser.error('Provide --serial/--mpi or put grial and grial_mpi in PATH')
    serial, mpi = Path(args.serial).resolve(), Path(args.mpi).resolve()
    for exe in [serial, mpi]:
        if not exe.is_file() or not os.access(exe, os.X_OK):
            parser.error('Missing executable: ' + str(exe))
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-container-'))
    if args.output:
        out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE', OMP_WAIT_POLICY='PASSIVE')
    reference = None
    checks = []
    try:
        for ranks in [0, 1, 2, 4]:
            label = 'openmp1' if ranks == 0 else f'mpi{ranks}'
            work = out / label
            work.mkdir()
            par = work / 'case.par'
            shutil.copyfile(ROOT / 'containers/sod.par', par)
            # Small packaging control only; oversubscription is not a benchmark.
            prefix = [] if ranks == 0 else ['mpirun', '--oversubscribe', '--bind-to', 'none', '-n', str(ranks)]
            program = serial if ranks == 0 else mpi
            for name, options in [('check', ['--check']), ('run', [])]:
                arguments = [*prefix, str(program), '--distributed', *options, str(par), '--grid', str(max(1, ranks)), '1']
                with (work / (name + '.log')).open('wb') as log:
                    result = subprocess.run(arguments, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
                if result.returncode:
                    raise RuntimeError(f'{label}/{name}: exit {result.returncode}\n' + (work / (name + '.log')).read_text()[-5000:])
                if name == 'check':
                    assert not (work / 'native').exists(), 'Parameter check wrote output'
            path = work / 'native/weno5_hlle/run_status.dat'
            receipt = dict(line.split('=', 1) for line in path.read_text().splitlines())
            assert receipt['status'] == 'completed' and int(receipt['ranks']) == max(1, ranks), receipt
            assert abs(float(receipt['physical_time']) - 0.02) < 1e-14, receipt
            values, events, shape = field_state(work / receipt['checkpoint'])
            assert shape == (64, 1, 1) and all(map(math.isfinite, values))
            if reference is None:
                reference = values, events, receipt
            prior, old_events, old_receipt = reference
            assert receipt['step'] == old_receipt['step']
            error = max(abs(u-v) / max(1.0, abs(u), abs(v)) for u, v in zip(prior, values))
            assert error <= 1e-11, (label, error)
            assert events[7:13] == old_events[7:13], 'Physical recovery counters differ'
            fields = work / 'native/weno5_hlle/mpi_fields'
            frames = collect_frames(fields, fields / 'frames.pvd')
            final = read_frame(fields / f'step_{int(receipt["step"]):010}', reconstruct=True)
            assert frames and tuple(final['shape']) == shape
            checks.append(dict(profile=label, ranks=max(1, ranks), threads=1, field_error=error,
                               fields_exact=values.tobytes() == prior.tobytes(), step=int(receipt['step']),
                               physical_time=float(receipt['physical_time']), visualization_valid=True))
        report = dict(passed=True, scope='short single-host CPU packaging check; not V01 or a benchmark',
                      tolerance=1e-11, checks=checks,
                      serial_sha256=hashlib.sha256(serial.read_bytes()).hexdigest(),
                      mpi_sha256=hashlib.sha256(mpi.read_bytes()).hexdigest())
    except Exception as error:
        (out / 'SUMMARY.json').write_text(json.dumps(dict(passed=False, error=str(error), checks=checks), indent=2) + '\n')
        raise
    (out / 'SUMMARY.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
