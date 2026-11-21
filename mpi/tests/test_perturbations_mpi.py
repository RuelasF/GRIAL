#!/usr/bin/env python3
"""Global-body selection, historical noise ordering and once-only event matrix."""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import tempfile


def rows(path):
    return [[float(x) for x in line.split()] for line in path.read_text().splitlines()
            if line.strip() and not line.startswith('#')]


def inspect(work, count):
    summary = []
    for rank in range(count):
        values = rows(work/f'perturbations_rank_{rank}.txt')[0]
        assert values[:2] == [rank, count] and values[3] == 54, values
        assert all(math.isfinite(x) and 0 <= x < 1e-11 for x in values[4:]), values
        cases = rows(work/f'perturbation_cases_rank_{rank}.txt')
        assert len(cases) == 54
        assert all(len(row) == 13 and all(math.isfinite(x) for x in row) for row in cases)
        assert [row[0] for row in cases] == list(range(1, 55))
        assert all(row[4] == (0 if row[5] == 2 else 1) for row in cases)
        assert all(max(row[8:]) < 1e-11 for row in cases)
        summary.append(dict(rank=rank, cases=54, field_error=values[4], budget_error=values[5],
                            forcing_error=values[6], observable_error=values[7],
                            exactly_equal_cases=sum(row[6] == 1 for row in cases)))
    event_count = 0
    for index in range(1, 55):
        directory = work/f'case_{index}'/'d03'
        event = directory/'mpi_perturbation_event.dat'
        if event.exists():
            data = rows(event)
            assert len(data) == 1 and len(data[0]) == 13
            assert all(math.isfinite(x) for x in data[0]) and data[0][7] > 0
            event_count += 1
        data = rows(directory/'mpi_conservation_budget.dat')
        assert all(len(row) == 62 for row in data)
        assert len(data) in [1, 2]
        assert all(a[0] < b[0] for a, b in zip(data, data[1:]))
        for line in data:
            scale = max(abs(x) for x in line[2:12])
            assert max(abs(x) for x in line[-5:])/scale < 1e-11
    assert event_count == 51, event_count
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--grid', type=int, nargs=2)
    args = parser.parse_args()
    probe = args.probe.resolve()
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-d03-'))
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE', OMP_WAIT_POLICY='PASSIVE')
    env.pop('OMPI_MPI_THREAD_LEVEL', None)
    grids = [args.grid] if args.grid else [(1, 1), (2, 1), (1, 2), (4, 1), (1, 4), (2, 2), (4, 2), (2, 4)]
    invalid = {
        'bad-amplitude': 'invalid perturbation controls',
        'bad-seed': 'perturbation controls differ between ranks',
        'bad-flag': 'perturbation controls differ between ranks',
        'bad-clock': 'perturbation controls differ between ranks',
        'bad-type': 'invalid perturbation controls',
        'bad-layout': 'inconsistent domain layout',
        'bad-state': 'unphysical state before perturbation',
        'bad-up': 'nonfinite physical observable state or metric',
        'bad-midstep': 'cannot sample observables inside an RK step',
        'bad-halo': 'not initialized',
        'bad-io': 'cannot publish unique distributed perturbation event',
    }
    runs = [(pr, pt, '') for pr, pt in grids]
    if not args.grid:
        runs += [(pr, pt, mode) for mode in invalid for pr, pt in [(1, 2), (2, 2)]]
    results = []
    for pr, pt, mode in runs:
        work = out/f'{mode or "perturbations"}_{pr}x{pt}';work.mkdir()
        command = ['mpirun', '--oversubscribe', '--bind-to', 'none', '-np', str(pr*pt),
                   str(probe), str(pr), str(pt), mode or 'report']
        with (work/'run.log').open('wb') as log:
            process = subprocess.Popen(command, cwd=work, env=env, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = process.wait(timeout=300 if not mode else 45)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL);process.wait()
                raise RuntimeError('D03 matrix timed out: '+str(work))
        if mode:
            text = (work/'run.log').read_text(errors='replace')
            assert code != 0 and invalid[mode] in text, (work, code, text[-4000:])
            if mode == 'bad-io':
                assert (work/'invalid/d03/mpi_perturbation_event.dat').read_text() == 'preserve previous evidence\n'
            result = dict(grid=[pr, pt], mode=mode, returncode=code, rejected=True)
        else:
            assert code == 0, (work, code, (work/'run.log').read_text()[-4000:])
            result = dict(grid=[pr, pt], passed=True, ranks=inspect(work, pr*pt))
        results.append(result)
        (out/'PROGRESS.json').write_text(json.dumps(results, indent=2)+'\n')
    (out/'VERIFIED.json').write_text(json.dumps(dict(passed=True, results=results), indent=2)+'\n')
    print('PERTURBATIONS MPI MATRIX PASSED:', out)


if __name__ == '__main__':
    main()
