#!/usr/bin/env python3
"""M04 native reconstruction/flux/static-RHS tests, without time evolution."""
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
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-m04-'))
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE')
    env.pop('OMPI_MPI_THREAD_LEVEL', None)
    results = []
    cases = [(n, '') for n in [1, 2, 4]] + [(n, mode) for mode in ['bad-offset', 'bad-total'] for n in [2, 4]]
    for n, mode in cases:
        name = f'{mode or "interfaces"}_{n}'
        work = out / name
        work.mkdir()
        command = [args.mpiexec, '--oversubscribe', '-n', str(n), str(probe), mode or 'report']
        start = time.monotonic()
        process = subprocess.Popen(command, cwd=work, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            output, _ = process.communicate(timeout=60)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            output, _ = process.communicate()
            (work / 'run.log').write_bytes(output)
            raise AssertionError('MPI interface test timed out: ' + name)
        (work / 'run.log').write_bytes(output)
        decoded = output.decode(errors='replace')
        assert (process.returncode == 0) == (not mode), (name, process.returncode, decoded[-5000:])
        expected = {'run.log'} | ({f'interface_rank_{r}.txt' for r in range(n)} if not mode else set())
        assert {p.name for p in work.iterdir()} == expected, 'Unexpected interface probe output'
        ranks = []
        if not mode:
            for expected_rank in range(n):
                row = (work / f'interface_rank_{expected_rank}.txt').read_text().split()
                assert len(row) == 11, row
                r, size, checks, fixtures = map(int, row[:4])
                errors = list(map(float, row[4:]))
                assert r == expected_rank and size == n and checks > 10000 and fixtures == 25, row
                assert errors[:4] == [0.0]*4, row
                assert 0 <= errors[4] <= 5e-15 and all(0 <= z <= 2e-13 for z in errors[5:]), row
                ranks.append(dict(rank=r, processes=size, assertions=checks, fixtures=fixtures,
                    state_error=errors[0], flux_error=errors[1], rhs_error=errors[2], shared_face_error=errors[3],
                    face_state_scaled_error=errors[4], face_flux_scaled_error=errors[5], face_rhs_scaled_error=errors[6]))
        else:
            assert 'invalid reconstruction domain metadata' in decoded, decoded[-4000:]
        results.append(dict(name=name, processes=n, returncode=process.returncode,
                            seconds=time.monotonic()-start, ranks=ranks))
    (out / 'VERIFIED.json').write_text(json.dumps(dict(passed=True, checks=results,
        configurations=14, geometries=5, profiles=5, solvers=['HLLE', 'HLLC'],
        time_evolution=False, gpu_execution=False), indent=2) + '\n')
    print('MPI INTERFACE MATRIX PASSED:', out)


if __name__ == '__main__':
    main()
