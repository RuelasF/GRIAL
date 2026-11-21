#!/usr/bin/env python3
"""Global CPU state, RK budgets and root-only diagnostic publication."""
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
    return [line.split() for line in path.read_text().splitlines() if line and not line.startswith('#')]


def inspect_outputs(work, pr, pt, baseline):
    checks = []
    expected_cases = [(f, c) for f in range(1, 9) for c in range(1, 5)] + [(3, 5), (3, 6)]
    actual = {}
    for fixture, config in expected_cases:
        directory = work / f'case_{fixture}_{config}' / 'd01'
        paths = ['mpi_global_diagnostics.dat', 'mpi_conservation_budget.dat', 'mpi_numerical_diagnostics.dat']
        assert {p.name for p in directory.iterdir()} == set(paths), directory
        data = [rows(directory / name) for name in paths]
        assert all(len(lines) == 2 for lines in data), 'Exactly one initial/final global sample'
        assert [len(lines[0]) for lines in data] == [28, 62, 26]
        assert all(float(lines[0][0]) == 0 and int(lines[0][1]) == 0 for lines in data)
        assert len({tuple(lines[-1][:2]) for lines in data}) == 1, directory
        state, budget, numerical = [lines[-1] for lines in data]
        steps = int(state[1])
        shape = [37, 19, 8]
        if fixture in [2, 3, 4]:
            shape = [400 if config >= 5 else 73, 1, 1]
        if fixture == 8:
            shape = [37, 1, 19]
        parts = [pr, pt] if shape[1] > 1 else [pr * pt, 1]
        cells = math.prod(shape)
        unique = sum((shape[d] + 1) * math.prod(shape[a] for a in range(3) if a != d)
                     for d in range(3) if shape[d] > 1)
        work_faces = sum((shape[d] + (parts[d] if d < 2 else 1)) *
                         math.prod(shape[a] for a in range(3) if a != d)
                         for d in range(3) if shape[d] > 1)
        counts = list(map(int, numerical[2:18]))
        assert int(state[25]) == cells and int(state[27]) == 0
        assert counts[7] == 3 * steps * cells
        assert sum(counts[:5]) == 3 * steps * work_faces
        assert list(map(int, numerical[23:])) == [unique, work_faces, 3 * steps * (work_faces - unique)]
        values = list(map(float, budget[2:]))
        assert all(math.isfinite(v) for v in values)
        scale = max(map(abs, values[:10]))
        assert all(value==0 for value in values[50:55]), 'unforced D01 case has a perturbation term'
        residual = max(map(abs, values[-5:])) / scale
        assert residual < 1e-11, (fixture, config, residual)
        mass_residual = abs(values[-5]) / max(abs(values[0]), abs(values[5]))
        assert mass_residual < 1e-11, (fixture, config, 'mass', mass_residual)
        j_residual = None
        if fixture in [5, 6]:
            j_residual = abs(values[-1]) / max(abs(values[4]), abs(values[9]))
            assert j_residual < 1e-11, (fixture, config, 'Kerr angular momentum', j_residual)
        error = 0.0
        if baseline:
            previous = baseline[(fixture, config)]
            assert state[:2] == previous[0][-1][:2]
            for kind, (now, before) in enumerate(zip(data, previous)):
                # Operation counts on shared faces vary; physical C2P counters do not.
                if kind == 2:
                    assert counts[7:13] == list(map(int, before[-1][2:18]))[7:13]
                    continue
                for row_now, row_before in zip(now, before):
                    count = 23 if kind == 0 else 60
                    for a, b in zip(row_now[2:2+count], row_before[2:2+count]):
                        error = max(error, abs(float(a) - float(b)) / scale)
            assert error < 1e-11, (fixture, config, error)
        actual[(fixture, config)] = data
        check = dict(fixture=fixture, config=config, shape=shape, steps=steps,
                     field_norms_vs_mpi1=error, budget_residual=residual,
                     mass_residual=mass_residual, angular_momentum_residual=j_residual,
                     c2p_calls=counts[7], extra_interface_evaluations=3*steps*(work_faces-unique))
        if fixture == 3 and config >= 5:
            check['mass_change_relative'] = (values[5] - values[0]) / abs(values[0])
            check['mass_budget_residual_relative'] = values[-5] / abs(values[0])
        checks.append(check)
    return checks, actual


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--mpiexec', default='mpirun')
    parser.add_argument('--output', type=Path)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--invalid-only', action='store_true')
    group.add_argument('--grid', nargs=2, type=int)
    args = parser.parse_args()
    probe = args.probe.resolve()
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-d01-'))
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE', OMP_WAIT_POLICY='PASSIVE')
    env.pop('OMPI_MPI_THREAD_LEVEL', None)
    grids = [tuple(args.grid)] if args.grid else [(1, 1), (2, 1), (1, 2), (4, 1), (1, 4), (2, 2), (4, 2), (2, 4)]
    invalid = {
        'bad-layout': 'diagnostic layout does not match',
        'bad-geometry': 'inconsistent diagnostic configuration',
        'bad-floors': 'inconsistent diagnostic spacing or floors',
        'bad-counters-disabled': 'require CPU numerical counters',
        'bad-restart': 'diagnostic restart requires history',
        'bad-inactive': 'parallel diagnostics are not initialized',
        'bad-stage': 'RK RHS budget stage out of sequence',
        'bad-midstep': 'cannot collect a partial RK budget',
        'bad-step': 'inconsistent RK time step',
        'bad-shape': 'diagnostic reference has wrong shape',
        'bad-duplicate': 'distributed diagnostic output failed',
        'bad-epoch': 'inconsistent diagnostic time epoch',
        'bad-grid': 'parallel diagnostic grid or controls changed',
        'bad-state': 'nonfinite diagnostic state or metric',
        'bad-count-range': 'diagnostic counters exceed the int64 reduction range',
    }
    runs = [] if args.invalid_only else [(pr, pt, '') for pr, pt in grids]
    if not args.grid:
        runs += [(pr, pt, mode) for mode in invalid for pr, pt in [(1, 2), (2, 2)]]
    results = []
    baseline = None
    for pr, pt, mode in runs:
        n = pr * pt
        work = out / f'{mode or "evolution"}_{pr}x{pt}'
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
            raise AssertionError('MPI diagnostic test timed out: ' + work.name)
        (work / 'run.log').write_bytes(output)
        decoded = output.decode(errors='replace')
        assert (process.returncode == 0) == (not mode), (work.name, process.returncode, decoded[-4500:])
        ranks = []
        cases = []
        if mode:
            assert invalid[mode] in decoded, (mode, decoded[-4500:])
        else:
            for r in range(n):
                row = (work / f'diagnostics_rank_{r}.txt').read_text().split()
                rank, size, checks, count, steps = map(int, row[:5])
                field, balance, state = map(float, row[5:])
                assert rank == r and size == n and checks > 1000 and count == 34 and steps > 1000, row
                assert field == 0 and balance < 1e-11 and state < 1e-11, row
                ranks.append(dict(rank=rank, checks=checks, cases=count, steps=steps,
                                  field_error=field, balance_error=balance, state_error=state))
            cases, actual = inspect_outputs(work, pr, pt, baseline)
            if (pr, pt) == (1, 1):
                baseline = actual
            assert 'PARALLEL DIAGNOSTICS TESTS PASSED' in decoded
        results.append(dict(grid=[pr, pt], mode=mode or 'evolution', ranks=ranks, cases=cases,
                            command=command, seconds=time.monotonic()-started))
        print(work.name + ': PASS', flush=True)
    (out / 'VERIFIED.json').write_text(json.dumps(dict(passed=True, checks=results,
        int64_control=True, scope='CPU global diagnostics and RK budgets; root-only text output; no restart/production IO/scaling'), indent=2)+'\n')
    print(f'DIAGNOSTICS MPI MATRIX PASSED: {out}')


if __name__ == '__main__':
    main()
