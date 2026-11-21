#!/usr/bin/env python3
"""M02 native MPI grid checks; no fluid time integration or halo exchanges."""
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
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-m02-'))
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE')
    env.pop('OMPI_MPI_THREAD_LEVEL', None)
    results = []
    cases = [(n, '') for n in [1, 2, 4]] + [(n, 'bad-one-rank') for n in [2, 4]] + [(4, 'bad-small')]
    for n, mode in cases:
        name = f'{mode or "grid"}_{n}'
        work = out / name
        work.mkdir()
        command = [args.mpiexec, '--oversubscribe', '-n', str(n), str(probe)]
        command.append(mode or 'report')
        start = time.monotonic()
        process = subprocess.Popen(command, cwd=work, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            output, _ = process.communicate(timeout=30)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            output, _ = process.communicate()
            (work / 'run.log').write_bytes(output)
            raise AssertionError('MPI grid test timed out: ' + name)
        (work / 'run.log').write_bytes(output)
        decoded = output.decode(errors='replace')
        assert (process.returncode == 0) == (not mode), (name, process.returncode, decoded[-4000:])
        expected_files = {'run.log'} | ({f'domain_rank_{r}.txt' for r in range(n)} if not mode else set())
        assert {p.name for p in work.iterdir()} == expected_files, 'Unexpected output from grid probe'
        ranks = []
        roundoff = []
        if not mode:
            # MPI output forwarding can split/interleave long Fortran stdout records.
            for expected_rank in range(n):
                row = (work / f'domain_rank_{expected_rank}.txt').read_text().split()
                assert len(row) == 7, row
                r, s, c = map(int, row[:3])
                p, u, pu, uu = map(float, row[3:])
                assert r == expected_rank and s == n and c > 1000, row
                assert 0 <= pu <= 4 and 0 <= uu <= 16, row
                ranks.append(dict(rank=r, processes=s, assertions=c))
                roundoff.append(dict(rank=r, primitive_error=p, conserved_error=u,
                                     primitive_ulps=pu, conserved_ulps=uu))
        else:
            assert 'collective validation failed' in decoded, decoded
        results.append(dict(name=name, processes=n, returncode=process.returncode,
                            seconds=time.monotonic()-start, ranks=ranks, khi_roundoff=roundoff))
    (out / 'VERIFIED.json').write_text(json.dumps(dict(passed=True, checks=results,
        grid_fixtures_per_rank=27, distributed_evolution=False), indent=2) + '\n')
    print('MPI DOMAIN TESTS PASSED:', out)


if __name__ == '__main__':
    main()
