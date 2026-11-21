#!/usr/bin/env python3
"""Native CPU MPI full-step evolution controls; no production writers or restart."""
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
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-m06-'))
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE')
    env.pop('OMPI_MPI_THREAD_LEVEL', None)
    invalid = {
        'rk-zero': 'invalid RK time step', 'rk-nan': 'invalid RK time step',
        'rk-step': 'inconsistent RK time step', 'rk-state': 'RK state arrays',
        'rk-halo': 'halo exchange is not initialized', 'rk-shape': 'RK state arrays',
        'rk-budget': 'radial mass budget requires D01', 'rk-thread': 'inside an OpenMP region',
    }
    results = []
    for n, mode in [(n, '') for n in [1, 2, 4]] + [(n, mode) for mode in invalid for n in [2, 4]]:
        name = f'{mode or "evolution"}_{n}'
        work = out / name
        work.mkdir()
        command = [args.mpiexec, '--oversubscribe', '-n', str(n), str(probe), mode or 'report']
        started = time.monotonic()
        process = subprocess.Popen(command, cwd=work, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            output, _ = process.communicate(timeout=360 if not mode else 30)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            output, _ = process.communicate()
            (work / 'run.log').write_bytes(output)
            raise AssertionError('MPI RK test timed out: ' + name)
        (work / 'run.log').write_bytes(output)
        decoded = output.decode(errors='replace')
        assert (process.returncode == 0) == (not mode), (name, process.returncode, decoded[-5000:])
        expected = {'run.log'} | ({name for r in range(n)
                                  for name in [f'rk_rank_{r}.txt', f'rk_cases_rank_{r}.txt']} if not mode else set())
        assert {p.name for p in work.iterdir()} == expected, 'Unexpected RK output'
        ranks = []
        if not mode:
            for expected_rank in range(n):
                row = (work / f'rk_rank_{expected_rank}.txt').read_text().split()
                assert len(row) == 8, row
                r, size, checks, cases, steps = map(int, row[:5])
                fields, dt, balance = map(float, row[5:])
                assert r == expected_rank and size == n and checks > 1000 and cases == 100 and steps > 1000, row
                assert fields == dt == 0 and balance < 1e-11, row
                ranks.append(dict(rank=r, processes=size, assertions=checks, evolution_cases=cases,
                                  steps=steps, field_error=fields, dt_error=dt, periodic_balance_error=balance))
                rows = (work / f'rk_cases_rank_{r}.txt').read_text().splitlines()[1:]
                assert len(rows) == 100
                for line in rows:
                    row = line.split()
                    assert len(row) == 34, row
                    integers = list(map(int, row[:24]))
                    values = list(map(float, row[24:]))
                    cells = integers[4]*integers[5]*integers[6]
                    assert integers[15] == 3*integers[7]*cells, 'C2P count excludes halo cells'
                    assert values[1] == values[2] == 0 and values[3] > 0 and values[4] > 0
        else:
            assert invalid[mode] in decoded, (mode, decoded[-4000:])
        results.append(dict(name=name, processes=n, returncode=process.returncode,
                            seconds=time.monotonic()-started, ranks=ranks))
    (out / 'VERIFIED.json').write_text(json.dumps(dict(passed=True, checks=results,
        fixtures=10, configurations=10, distributed_fluid_evolution=True,
        production_io=False, distributed_checkpoint=False), indent=2) + '\n')
    print('MPI RK MATRIX PASSED:', out)


if __name__ == '__main__':
    main()
