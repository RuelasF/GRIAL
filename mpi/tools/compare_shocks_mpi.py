#!/usr/bin/env python3
"""Compare developed SRHD shocks across serial and 1/2/4 MPI processes."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import numpy as np

METHODS = ('godunov', 'tvd', 'weno3', 'mp5', 'weno5')
FIELD_NAMES = ('rho', 'p', 'vx', 'vy', 'vz', 'D', 'tau', 'Sx', 'Sy', 'Sz')
CASES = {'sod': (1., .125, 1., .1, 1.4), 'strong': (1., 1., 1000., .01, 5/3)}
TIME = .4


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_run(work, processes, nx):
    chunks, traces, events, repairs = [], [], [], []
    for rank in range(processes):
        metadata = (work/f'profile_metadata_rank_{rank}.dat').read_text().splitlines()
        integers = list(map(int, metadata[0].split()))
        values = list(map(float, metadata[1].split()))
        counts = np.array(list(map(int, metadata[2].split())), dtype=np.int64)
        assert len(integers) == 8 and len(values) == 8 and len(counts) == 16
        r, size, global_n, offset, local_n, steps, method, solver = integers
        assert (r, size, global_n) == (rank, processes, nx)
        assert values[0] == TIME and values[2] == .3
        assert counts[7] == 3*steps*local_n
        data = np.loadtxt(work/f'profile_rank_{rank}.dat', ndmin=2)
        trace = np.loadtxt(work/f'profile_trace_rank_{rank}.dat', ndmin=2)
        assert data.shape == (local_n, 12) and trace.shape == (steps, 3)
        assert np.array_equal(data[:, 0], np.arange(offset+1, offset+local_n+1))
        assert np.array_equal(trace[:, 0], np.arange(1, steps+1))
        assert trace[-1, 1] == TIME and np.all(np.diff(trace[:, 1]) > 0)
        assert np.all(np.isfinite(data)) and np.all(np.isfinite(trace)) and np.all(trace[:, 2] > 0)
        assert np.all(data[:, 2:4] > 0) and np.all(np.sum(data[:, 4:7]**2, axis=1) < 1)
        chunks.append(data); traces.append(trace); events.append(counts)
        repairs.append(np.array(list(map(float, metadata[3].split()))))
        if rank:
            assert traces[0].tobytes() == trace.tobytes(), 'different rank time traces'
    data = np.concatenate(chunks)
    data = data[np.argsort(data[:, 0])]
    assert data.shape == (nx, 12) and np.array_equal(data[:, 0], np.arange(1, nx+1))
    assert np.array_equal(data[:, 1], (np.arange(nx)+.5)*(1/nx))
    counts = np.sum(events, axis=0)
    assert counts[7] == 3*steps*nx
    return data, traces[0], counts, np.sum(repairs, axis=0)


def execute(exe, work, problem, method, solver, nx, route, processes, debug):
    work.mkdir(parents=True)
    command = [str(exe), problem, method, solver, str(nx), 'profile', route]
    if route == 'mpi':
        command = ['mpirun', '--oversubscribe', '--bind-to', 'none', '-n', str(processes), *command]
    env = dict(os.environ, OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE')
    env.pop('OMPI_MPI_THREAD_LEVEL', None)
    started = time.monotonic()
    process = subprocess.Popen(command, cwd=work, env=env, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, start_new_session=True)
    try:
        output, _ = process.communicate(timeout=180 if debug else 90)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        output, _ = process.communicate()
        (work/'run.log').write_bytes(output)
        raise AssertionError(f'Shock run timed out: {work}')
    (work/'run.log').write_bytes(output)
    assert process.returncode == 0 and b'SHOCK_PROFILE PASSED' in output, (work, output[-4000:])
    expected = {'run.log'} | {f'profile{kind}_rank_{r}.dat' for r in range(processes)
                                for kind in ['', '_trace', '_metadata']}
    if problem == 'shu':
        expected |= {f'profile_initial_rank_{r}.dat' for r in range(processes)}
    assert {p.name for p in work.iterdir()} == expected
    return (*read_run(work, processes, nx), time.monotonic()-started)


def write_table(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def plot(out, profiles, exact):
    os.environ['MPLCONFIGDIR'] = str(out/'matplotlib_cache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11})
    styles = {1: dict(color='#1764ad', lw=1.3),
              2: dict(color='#e1812c', lw=1.1, ls='--'),
              4: dict(color='#268559', lw=0, marker='o', ms=3, fillstyle='none', markevery=13)}
    labels = (r'$\rho$', r'$p$', r'$v_x$')
    positions = np.linspace(.0001, .9999, 3000)
    exact_lines = {}
    for problem in CASES:
        states = [exact[problem].sample((x-.5)/TIME) for x in positions]
        exact_lines[problem] = np.array([[state.rho, state.p, state.v] for state in states])

    def panel(ax, problem, method, solver, column):
        ax.plot(positions, exact_lines[problem][:, column], color='#242424', lw=1.1, label='Exacta', zorder=1)
        for n in (1, 2, 4):
            data = profiles['opt', problem, method, solver, 400, n]
            ax.plot(data[:, 1], data[:, 2+column], label=f'MPI {n}', **styles[n], zorder=2+n)
        ax.set(xlabel='x', ylabel=labels[column], xlim=(0, 1))
        ax.grid(alpha=.2)
        if problem == 'strong' and column == 1:
            ax.set_yscale('log')

    for problem in CASES:
        fig, axes = plt.subplots(2, 3, figsize=(13, 6.5), constrained_layout=True,
                                 gridspec_kw={'height_ratios': [3, 1.2]})
        for column in range(3):
            panel(axes[0, column], problem, 'weno5', 'hlle', column)
            base = profiles['opt', problem, 'weno5', 'hlle', 400, 1]
            maximum = 0.
            for n in (2, 4):
                data = profiles['opt', problem, 'weno5', 'hlle', 400, n]
                difference = data[:, 2+column]-base[:, 2+column]
                maximum = max(maximum, float(np.max(np.abs(difference))))
                axes[1, column].plot(data[:, 1], difference, label=f'MPI {n} − MPI 1', **styles[n])
            axes[1, column].set(xlabel='x', ylabel=r'$\Delta$'+labels[column], xlim=(0, 1))
            axes[1, column].grid(alpha=.2)
            if maximum == 0:
                axes[1, column].set(ylim=(-1, 1), yticks=[0])
            axes[1, column].text(.03, .85, f'máx |diferencia| = {maximum:.1e}',
                                transform=axes[1, column].transAxes, fontsize=9)
        axes[0, 0].legend(fontsize=9, ncol=2)
        axes[1, 0].legend(fontsize=8, loc='lower left')
        title = 'Sod' if problem == 'sod' else 'Choque fuerte'
        fig.suptitle(f'{title}: WENO5–HLLE, sensor activo, N=400, t=0.4\n'
                     'Perfiles y diferencias al cambiar el número de procesos MPI', fontsize=13)
        fig.savefig(out/f'{problem}_weno5_hlle_mpi.png', dpi=200)
        fig.savefig(out/f'{problem}_weno5_hlle_mpi.pdf')
        plt.close(fig)
        for solver in ('hlle', 'hllc'):
            fig, axes = plt.subplots(5, 3, figsize=(13, 14), constrained_layout=True)
            for row, method in enumerate(METHODS):
                for column in range(3):
                    panel(axes[row, column], problem, method, solver, column)
                    axes[row, column].set_title(f'{method.upper()} · {labels[column]}')
            axes[0, 0].legend(fontsize=8, ncol=2)
            fig.suptitle(f'{title}: {solver.upper()}, N=400, t=0.4 · sensor activo\n'
                         'Exacta y MPI 1/2/4: líneas continuas, punteadas y marcadores', fontsize=13)
            fig.savefig(out/f'{problem}_{solver}_reconstructores.png', dpi=180)
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--serial-exe', type=Path, required=True)
    parser.add_argument('--mpi-exe', type=Path, required=True)
    parser.add_argument('--debug-serial-exe', type=Path, required=True)
    parser.add_argument('--debug-mpi-exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--exact-reference', type=Path,
        default=Path(__file__).resolve().parents[2]/'docs/srhd_1d/exact_riemann.py')
    args = parser.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location('grial_exact_reference', args.exact_reference.resolve())
    reference = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = reference; spec.loader.exec_module(reference)
    exact = {}; means = {}; reference_metadata = {}
    for problem, (rl, rr, pl, pr, gamma) in CASES.items():
        solution = reference.Riemann(reference.State(rl, pl), reference.State(rr, pr), gamma)
        exact[problem] = solution
        assert min(solution.breaks(TIME)) > 0 and max(solution.breaks(TIME)) < 1
        if problem == 'strong':
            assert abs(solution.v_star-.960) < .002 and abs(solution.right_edges[0]-.986) < .002
        reference_metadata[problem] = dict(gamma=gamma, left=[rl, pl, 0.], right=[rr, pr, 0.],
            wave_positions=solution.breaks(TIME), p_star=solution.p_star, v_star=solution.v_star)
        for nx in (400, 401):
            means[problem, nx] = solution.cell_means(nx, TIME, quadrature_order=12)
            expected = np.array([(rl+rr)/2, (pl+pr)/2/(gamma-1), (pl-pr)*TIME, 0, 0])
            integrated = np.array([np.mean(means[problem, nx][name]) for name in ('D', 'tau', 'S')])
            assert np.max(np.abs(integrated-expected[:3])/np.maximum(1., np.abs(expected[:3]))) < 2e-9
    profiles, checks, errors = {}, [], []
    binary_hashes = {}
    for build in ('opt', 'debug'):
        serial_exe = (args.serial_exe if build == 'opt' else args.debug_serial_exe).resolve()
        mpi_exe = (args.mpi_exe if build == 'opt' else args.debug_mpi_exe).resolve()
        binary_hashes[build] = dict(serial=sha(serial_exe), mpi=sha(mpi_exe))
        configurations = [(p, m, s, 400) for p in CASES for m in
            (METHODS if build == 'opt' else ('mp5', 'weno5')) for s in
            (('hlle', 'hllc') if build == 'opt' else ('hlle',))]
        configurations += [(p, m, 'hlle', 401) for p in CASES for m in ('mp5', 'weno5')]
        for problem, method, solver, nx in configurations:
            folder = out/'runs'/build/f'{problem}_{method}_{solver}_n{nx}'
            baseline, trace, _, _, _ = execute(serial_exe, folder/'serial', problem, method, solver, nx,
                                               'serial', 1, build == 'debug')
            for n in (1, 2, 4):
                data, current_trace, counts, repair_delta, seconds = execute(mpi_exe, folder/f'mpi{n}',
                    problem, method, solver, nx, 'mpi', n, build == 'debug')
                difference = np.max(np.abs(data[:, 2:]-baseline[:, 2:]), axis=0)
                assert data.tobytes() == baseline.tobytes(), (build, problem, method, solver, nx, n, difference)
                assert current_trace.tobytes() == trace.tobytes(), 'different full time sequence'
                profiles[build, problem, method, solver, nx, n] = data
                rl, rr, pl, pr, gamma = CASES[problem]
                expected = np.array([(rl+rr)/2, (pl+pr)/2/(gamma-1), (pl-pr)*TIME, 0, 0])
                budget = np.mean(data[:, 7:12], axis=0)-expected
                scaled_budget = float(np.max(np.abs(budget)/np.maximum(1., np.abs(expected))))
                edge = data[[0, -1]]
                flux = np.column_stack([edge[:, 7]*edge[:, 4], edge[:, 9]-edge[:, 7]*edge[:, 4],
                    edge[:, 9]*edge[:, 4]+edge[:, 3], edge[:, 10]*edge[:, 4], edge[:, 11]*edge[:, 4]])
                constant_flux = np.zeros((2, 5)); constant_flux[:, 2] = [pl, pr]
                boundary_change = float(np.max(np.abs(flux-constant_flux)/np.maximum(1., np.abs(constant_flux))))
                # This diagnostic assumes constant boundary states. Keep its threshold and report failures.
                budget_within = scaled_budget < 1e-11
                if not budget_within:
                    assert boundary_change > 1e-11, ('unexplained constant-boundary budget residual',
                        build, problem, method, solver, nx, scaled_budget)
                assert np.sum(counts[8:14]) == 0 and np.all(repair_delta == 0), 'shock state repaired'
                for name, column, exact_name in [('D', 7, 'D'), ('tau', 8, 'tau'), ('Sx', 9, 'S')]:
                    error = data[:, column]-means[problem, nx][exact_name]
                    errors.append(dict(build=build, problem=problem, method=method, solver=solver, nx=nx,
                        processes=n, variable=name, L1=float(np.mean(np.abs(error))),
                        L2=float(np.sqrt(np.mean(error**2))), Linf=float(np.max(np.abs(error)))))
                checks.append(dict(build=build, problem=problem, method=method, solver=solver, nx=nx,
                    processes=n, steps=len(trace), time=TIME, field_max_abs_error=float(np.max(difference)),
                    time_sequence_identical=True, fields_byte_identical=True, seconds=seconds,
                    c2p_calls=int(counts[7]), c2p_repairs=int(np.sum(counts[8:14])),
                    events=counts.tolist(), repair_delta_u=repair_delta.tolist(),
                    analytic_budget_scaled_error=scaled_budget, analytic_budget_within_1e11=budget_within,
                    final_boundary_flux_scaled_change=boundary_change))
            print(f'{build} {problem} {method} {solver} N={nx}: MPI 1/2/4 identical; steps={len(trace)}', flush=True)
    write_table(out/'errors_exact.tsv', errors)
    write_table(out/'comparisons.tsv', [{k: v for k, v in row.items() if not isinstance(v, list)} for row in checks])
    (out/'REFERENCE.json').write_text(json.dumps(dict(time=TIME, domain=[0, 1], cfl=.3,
        shock_sensor=True, tvd_limiter='MC', odd_grid_initialization='exact initial cell fraction at x=0.5',
        exact_reference_sha256=sha(args.exact_reference), cases=reference_metadata), indent=2)+'\n')
    plot(out, profiles, exact)
    report = dict(passed=True, checks=checks, trajectories=128, mpi_comparisons=len(checks),
        serial_references=32, field_max_abs_error=max(row['field_max_abs_error'] for row in checks),
        analytic_budget_max_scaled_error=max(row['analytic_budget_scaled_error'] for row in checks),
        binary_sha256=binary_hashes, time=TIME, nx=[400, 401], processes=[1, 2, 4],
        c2p_repairs=sum(row['c2p_repairs'] for row in checks), production_mpi_io=False,
        scope='MPI versus serial field and time equivalence; analytic budget diagnostic reported separately',
        analytic_budget_threshold=1e-11,
        analytic_budget_flagged_comparisons=sum(not row['analytic_budget_within_1e11'] for row in checks),
        distributed_restart_tested=False, scalability_benchmark=False)
    assert len(checks) == 96
    (out/'VERIFIED.json').write_text(json.dumps(report, indent=2)+'\n')
    print('MPI DEVELOPED SHOCKS PASSED:', out)


if __name__ == '__main__':
    main()
