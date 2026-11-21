#!/usr/bin/env python3
"""CPU r/theta blocks against a full-domain oracle, with complete phi."""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--mpiexec', default='mpirun')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--invalid-only', action='store_true', help='Run only the collective rejection controls')
    args = parser.parse_args()
    probe = args.probe.resolve()
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-t01-'))
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE')
    env.pop('OMPI_MPI_THREAD_LEVEL', None)
    grids = [(1, 1), (2, 1), (1, 2), (4, 1), (1, 4), (2, 2), (4, 2), (2, 4)]
    invalid = {
        'bad-small': 'polar block too small',
        'bad-grid': 'process grid must cover',
        'bad-periodic': 'polar periodicity does not match',
        'bad-config': 'inconsistent halo configuration',
        'bad-case': 'polar decomposition is unsupported',
        'bad-odd': 'polar boundaries require even global nz',
        'bad-layout': 'inconsistent domain layout',
        'bad-metadata': 'global grid metadata changed',
    }
    results = []
    runs = [] if args.invalid_only else [(pr, pt, '') for pr, pt in grids]
    runs += [(pr, pt, mode) for mode in invalid for pr, pt in [(1, 2), (2, 2)]]
    for pr, pt, mode in runs:
        n = pr * pt
        name = f'{mode or "evolution"}_{pr}x{pt}'
        work = out / name
        work.mkdir()
        command = [args.mpiexec, '--oversubscribe', '--bind-to', 'none', '-n', str(n),
                   str(probe), str(pr), str(pt), mode or 'report']
        started = time.monotonic()
        process = subprocess.Popen(command, cwd=work, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            output, _ = process.communicate(timeout=600 if not mode else 30)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            output, _ = process.communicate()
            (work / 'run.log').write_bytes(output)
            raise AssertionError('MPI rtheta test timed out: ' + name)
        (work / 'run.log').write_bytes(output)
        decoded = output.decode(errors='replace')
        assert (process.returncode == 0) == (not mode), (name, process.returncode, decoded[-5000:])
        ranks = []
        expected = {'run.log'}
        if not mode:
            for r in range(n):
                expected.update((f'rtheta_rank_{r}.txt', f'rtheta_cases_rank_{r}.txt'))
                row = (work / f'rtheta_rank_{r}.txt').read_text().split()
                assert len(row) == 10, row
                rank, size, checks, cases, steps = map(int, row[:5])
                fields, dt, balance, halos, host = map(float, row[5:])
                assert rank == r and size == n and checks > 1000 and cases == 80 and steps >= 80, row
                assert all(math.isfinite(v) for v in [fields, dt, balance, halos, host]), row
                assert fields == dt == halos == 0.0 and balance < 1e-11 and host <= 2e-13, row
                lines = (work / f'rtheta_cases_rank_{r}.txt').read_text().splitlines()
                assert len(lines) == 81, (name, r, len(lines))
                ranks.append({'rank': rank, 'checks': checks, 'cases': cases, 'steps': steps,
                              'field_error': fields, 'dt_error': dt, 'balance_error': balance,
                              'halo_error': halos, 'host_rhs_error': host})
            assert 'RTHETA TESTS PASSED' in decoded
        else:
            assert invalid[mode] in decoded, (name, decoded[-5000:])
        assert {p.name for p in work.iterdir()} == expected, ('Unexpected output', name)
        results.append({'grid': [pr, pt], 'mode': mode or 'evolution', 'ranks': ranks,
                        'seconds': time.monotonic() - started, 'command': command})
        print(name + ': PASS', flush=True)
    (out / 'VERIFIED.json').write_text(json.dumps({'cases': results, 'count': len(results),
        'scope': 'CPU local domain/cache/halos/RK; phi complete; not production IO or scaling'}, indent=2) + '\n')
    print(f'{len(results)} RTHETA MPI CHECKS PASSED; {out}')


if __name__ == '__main__':
    main()
