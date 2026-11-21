#!/usr/bin/env python3
"""PPI coefficients, Finn--Evans tensor and accretion reduced to rank zero."""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import time


def rows(path):
    return [[float(v) for v in line.split()] for line in path.read_text().splitlines()
            if line.strip() and not line.startswith('#')]


def inspect_outputs(work):
    checked = []
    for kind, count in [('synthetic', 5), ('evolution', 8)]:
        for fixture in range(1, count+1):
            directory = work / f'{kind}_{fixture}'
            files = [('ppi_modes.dat', 18), ('GW_signal.dat', 3), ('m_dot.dat', 2)]
            if kind == 'evolution' and fixture > 4:
                files = files[1:]
            expected = {f'mpi_{name}' for name, _ in files}
            if kind == 'evolution':
                expected |= {'mpi_global_diagnostics.dat', 'mpi_conservation_budget.dat', 'mpi_numerical_diagnostics.dat'}
            assert {p.name for p in (directory/'distributed/d02').iterdir()} == expected, directory
            integrals = rows(directory/'serial/d02/reference_integrals.dat')
            for name, columns in files:
                full = rows(directory/'serial/d02'/name)
                local = rows(directory/'distributed/d02'/f'mpi_{name}')
                assert len(full) == len(local) == (12 if kind == 'synthetic' else 1)
                assert all(len(row) == columns for row in local)
                assert all(math.isfinite(x) for row in local for x in row)
                assert [r[0] for r in full] == [r[0] for r in local]
                assert all(a[0] < b[0] for a, b in zip(local, local[1:]))
                error = 0.0
                if name == 'ppi_modes.dat':
                    scale = max(max(abs(r[1]) for r in full), 1e-10)
                    for a, b in zip(local, full):
                        error = max(error, abs(a[1]-b[1])/scale)
                        for m in range(1, 5):
                            base = 2 + 4*(m-1)
                            error = max(error, *(abs(a[j]-b[j])/scale for j in [base, base+1]),
                                        abs(a[base+2]-b[base+2]))
                            if b[base+2] > 1e-10:
                                error = max(error, 2*abs(math.sin(0.5*m*(a[base+3]-b[base+3]))))
                else:
                    scale = max(max(abs(x) for row in full for x in row[1:]), 1e-10)
                    if name == 'GW_signal.dat':
                        # Axisymmetry cancels TT components; use the unreduced tensor's scale.
                        scale = max(scale, max(abs(x) for row in integrals for x in row[10:13])/1000.0)
                    error = max(abs(a-b)/scale for ra, rb in zip(local, full) for a, b in zip(ra[1:], rb[1:]))
                assert error < 1e-11, (directory, name, error)
                checked.append(dict(case=f'{kind}_{fixture}', file=name, samples=len(local), error=error))
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--mpiexec', default='mpirun')
    parser.add_argument('--grid', type=int, nargs=2)
    args = parser.parse_args()
    probe = args.probe.resolve()
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-d02-'))
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE', OMP_WAIT_POLICY='PASSIVE')
    env.pop('OMPI_MPI_THREAD_LEVEL', None)
    grids = [args.grid] if args.grid else [(1, 1), (2, 1), (1, 2), (4, 1), (1, 4), (2, 2), (4, 2), (2, 4)]
    invalid = {
        'bad-flags': 'observable controls differ between ranks',
        'bad-spin': 'observable controls differ between ranks',
        'bad-probe': 'invalid distributed accretion-probe ownership',
        'bad-state': 'nonfinite physical observable state or metric',
        'bad-time': 'observable time differs from RK epoch',
        'bad-midstep': 'cannot sample observables inside an RK step',
        'bad-root-destination': 'diagnostic reduction destination differs between ranks',
        'bad-duplicate': 'cannot publish root observable',
        'bad-existing': 'cannot publish root observable',
        'bad-io': 'cannot publish root observable',
    }
    runs = [(pr, pt, '') for pr, pt in grids]
    if not args.grid:
        runs += [(pr, pt, mode) for mode in invalid for pr, pt in [(1, 2), (2, 2)]]
    results = []
    for pr, pt, mode in runs:
        work = out/f'{mode or "observables"}_{pr}x{pt}'
        work.mkdir()
        command = [args.mpiexec, '--oversubscribe', '--bind-to', 'none', '-n', str(pr*pt),
                   str(probe), str(pr), str(pt), mode or 'report']
        started = time.monotonic()
        process = subprocess.Popen(command, cwd=work, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            log, _ = process.communicate(timeout=360 if not mode else 30)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            log, _ = process.communicate()
            (work/'run.log').write_bytes(log)
            raise AssertionError(f'MPI observable timeout: {work}')
        (work/'run.log').write_bytes(log)
        decoded = log.decode(errors='replace')
        assert (process.returncode == 0) == (not mode), (work, process.returncode, decoded[-6000:])
        checks = []; ranks = []
        if mode:
            assert invalid[mode] in decoded, (work, decoded[-6000:])
            if mode == 'bad-existing':
                assert (work/'invalid_output/d02/mpi_ppi_modes.dat').read_text().strip() == 'preserve'
        else:
            assert 'PARALLEL OBSERVABLES TESTS PASSED' in decoded
            checks = inspect_outputs(work)
            for rank in range(pr*pt):
                values = (work/f'observables_rank_{rank}.txt').read_text().split()
                r, size, assertions, cases = map(int, values[:4])
                field, observable, analytic = map(float, values[4:])
                assert r == rank and size == pr*pt and assertions > 100 and cases == 13
                assert field == 0 and max(observable, analytic) < 1e-11
                ranks.append(dict(rank=rank, assertions=assertions, cases=cases, field=field,
                                  observable=observable, analytic=analytic))
        results.append(dict(grid=[pr, pt], mode=mode or 'observables', ranks=ranks, files=checks,
                            command=command, seconds=time.monotonic()-started))
        print(work.name+': PASS', flush=True)
    (out/'VERIFIED.json').write_text(json.dumps(dict(passed=True, checks=results,
        threshold=1e-11, writer_rank=0, local_integrals=True, scope='CPU D02; no distributed restart/production/scaling'), indent=2)+'\n')
    print(f'OBSERVABLES MPI MATRIX PASSED: {out}')


if __name__ == '__main__':
    main()
