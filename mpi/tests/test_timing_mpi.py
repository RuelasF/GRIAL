#!/usr/bin/env python3
"""M05 local MPI CFL/epoch/event/stop tests; no distributed fluid evolution."""
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
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-m05-'))
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE')
    env.pop('OMPI_MPI_THREAD_LEVEL', None)
    invalid = {
        'cfl-spacing': 'CFL context', 'cfl-nan-spacing': 'CFL context',
        'cfl-courant': 'inconsistent CFL context', 'cfl-nan-courant': 'CFL context',
        'cfl-axis': 'inconsistent CFL context', 'cfl-cache': 'CFL context', 'cfl-shape': 'CFL context',
        'cfl-nan-metric': 'invalid metric data in CFL cells', 'cfl-lapse': 'invalid metric data in CFL cells',
        'cfl-inverse': 'invalid metric data in CFL cells', 'cfl-shift': 'invalid metric data in CFL cells',
        'cfl-thread': 'inside an OpenMP region',
        'time-epoch': 'inconsistent time epoch', 'time-step': 'inconsistent time epoch',
        'time-target': 'inconsistent time epoch', 'time-calendar': 'inconsistent time epoch',
        'time-enabled': 'inconsistent time epoch', 'time-nan': 'invalid time epoch',
        'time-calendar-nan': 'invalid time epoch',
        'clip-stagnant': 'does not advance physical time', 'clip-nan': 'invalid temporal step',
        'rate-overflow': 'CFL propagation rate overflow',
    }
    cases = [(n, '') for n in [1, 2, 4]] + [(n, mode) for mode in invalid for n in [2, 4]]
    results = []
    for n, mode in cases:
        name = f'{mode or "timing"}_{n}'
        work = out / name
        work.mkdir()
        command = [args.mpiexec, '--oversubscribe', '-n', str(n), str(probe), mode or 'report']
        started = time.monotonic()
        process = subprocess.Popen(command, cwd=work, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            output, _ = process.communicate(timeout=30)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            output, _ = process.communicate()
            (work / 'run.log').write_bytes(output)
            raise AssertionError('MPI temporal test timed out: ' + name)
        (work / 'run.log').write_bytes(output)
        decoded = output.decode(errors='replace')
        assert (process.returncode == 0) == (not mode), (name, process.returncode, decoded[-5000:])
        expected = {'run.log'} | ({f'timing_rank_{r}.txt' for r in range(n)} if not mode else set())
        assert {p.name for p in work.iterdir()} == expected, 'Unexpected temporal probe output'
        ranks = []
        if not mode:
            for expected_rank in range(n):
                row = (work / f'timing_rank_{expected_rank}.txt').read_text().split()
                assert len(row) == 7, row
                r, size, checks, steps = map(int, row[:4])
                error, wrong_2d, wrong_3d = map(float, row[4:])
                assert r == expected_rank and size == n and checks >= 50 and steps == 8 and error == 0, row
                if n > 1:
                    assert wrong_2d > 0.4/100 and wrong_3d > 0.4/340, row
                ranks.append(dict(rank=r, processes=size, assertions=checks, timeline_steps=steps,
                    cfl_error=error, minimum_local_dt_2d=wrong_2d, minimum_local_dt_3d=wrong_3d))
        else:
            assert invalid[mode] in decoded, (mode, decoded[-4000:])
        results.append(dict(name=name, processes=n, returncode=process.returncode,
                            seconds=time.monotonic()-started, ranks=ranks))
    (out / 'VERIFIED.json').write_text(json.dumps(dict(passed=True, checks=results,
        synthetic_cfl_cases=6, metric_geometries=5, correct_dt_2d=0.4/100, correct_dt_3d=0.4/340,
        distributed_fluid_evolution=False), indent=2) + '\n')
    print('MPI TIMING MATRIX PASSED:', out)


if __name__ == '__main__':
    main()
