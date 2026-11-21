#!/usr/bin/env python3
"""M03 native primitive halo tests; no distributed time evolution."""
import argparse
import json
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
    args = parser.parse_args()
    probe = args.probe.resolve()
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-m03-'))
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE')
    env.pop('OMPI_MPI_THREAD_LEVEL', None)
    invalid = {'bad-layout': 'collective validation failed',
               'bad-periodic': 'halo periodicity does not match',
               'bad-global-shape': 'inconsistent halo configuration',
               'bad-array': 'wrong shape',
               'bad-thread': 'inside an OpenMP region',
               'bad-lifecycle': 'not initialized',
               'bad-count': 'MPI count range',
               'bad-cache': 'require the local metric cache',
               'bad-injector': 'requires its injector',
               'bad-polar': 'require even global nz'}
    cases = [(n, '') for n in [1, 2, 4]] + [(n, mode) for mode in invalid for n in [2, 4]]
    results = []
    for n, mode in cases:
        name = f'{mode or "halos"}_{n}'
        work = out / name
        work.mkdir()
        command = [args.mpiexec, '--oversubscribe', '-n', str(n), str(probe), mode or 'report']
        start = time.monotonic()
        process = subprocess.Popen(command, cwd=work, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            output, _ = process.communicate(timeout=45)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            output, _ = process.communicate()
            (work / 'run.log').write_bytes(output)
            raise AssertionError('MPI halo test timed out: ' + name)
        (work / 'run.log').write_bytes(output)
        decoded = output.decode(errors='replace')
        assert (process.returncode == 0) == (not mode), (name, process.returncode, decoded[-5000:])
        expected_files = {'run.log'} | ({f'halo_rank_{r}.txt' for r in range(n)} if not mode else set())
        assert {p.name for p in work.iterdir()} == expected_files, 'Unexpected output from halo probe'
        ranks = []
        if not mode:
            for expected_rank in range(n):
                row = (work / f'halo_rank_{expected_rank}.txt').read_text().split()
                assert len(row) == 5, row
                r, s, c, fixtures = map(int, row[:4])
                error = float(row[4])
                assert r == expected_rank and s == n and c > 50000, row
                assert fixtures == (93 if n == 1 else 90) and error == 0, row
                ranks.append(dict(rank=r, processes=s, assertions=c, fixtures=fixtures, max_abs_error=error))
        else:
            assert invalid[mode] in decoded, (mode, decoded[-4000:])
        results.append(dict(name=name, processes=n, returncode=process.returncode,
                            seconds=time.monotonic()-start, ranks=ranks))
    (out / 'VERIFIED.json').write_text(json.dumps(dict(passed=True, checks=results,
        halo_depths=[1, 2, 3], physical_cases=list(range(1, 16)), aliases='names and numeric strings',
        epochs=2, distributed_evolution=False), indent=2) + '\n')
    print('MPI HALO TESTS PASSED:', out)


if __name__ == '__main__':
    main()
