#!/usr/bin/env python3
"""Shu--Osher profiles and initial/final MPI differences with 1/2/4/8 ranks."""
import argparse
import json
import os
from pathlib import Path

import numpy as np

from compare_shocks_mpi import execute, sha, write_table

METHODS = ('godunov', 'tvd', 'weno3', 'mp5', 'weno5')
PROCESSES = (1, 2, 4, 8)
FIELDS = ('rho', 'p', 'vx', 'vy', 'vz', 'D', 'tau', 'Sx', 'Sy', 'Sz')


def initial_profile(work, processes, nx):
    chunks = [np.loadtxt(work/f'profile_initial_rank_{rank}.dat', ndmin=2) for rank in range(processes)]
    data = np.concatenate(chunks)
    data = data[np.argsort(data[:, 0])]
    assert data.shape == (nx, 12) and np.array_equal(data[:, 0], np.arange(1, nx+1))
    assert np.all(np.isfinite(data))
    expected_rho = np.where(data[:, 1] < .5, 5., 2.+.3*np.sin(50*data[:, 1]))
    expected_p = np.where(data[:, 1] < .5, 50., 5.)
    # A new IC audit, bounded before running; different vectorized libm lanes may round differently.
    assert np.max(np.abs(data[:, 2]-expected_rho)/np.maximum(1., np.abs(expected_rho))) <= 16*np.finfo(float).eps
    assert np.array_equal(data[:, 3], expected_p) and np.all(data[:, 4:7] == 0)
    return data


def make_plots(out, profiles):
    os.environ['MPLCONFIGDIR'] = str(out/'matplotlib_cache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10})
    styles = {1: dict(color='#1764ad', lw=1.4),
              2: dict(color='#e1812c', lw=1.2, ls='--'),
              4: dict(color='#268559', lw=1., ls='-.'),
              8: dict(color='#8755ad', lw=0, marker='s', ms=3, fillstyle='none', markevery=9)}
    labels = (r'$\rho$', r'$p$', r'$v_x$')

    def panel(ax, method, solver, column, zoom=False):
        baseline = profiles['opt', method, solver, 400, 0]
        ax.plot(baseline[:, 1], baseline[:, 2+column], ':', color='#262626', lw=1., label='Serial')
        for n in PROCESSES:
            a = profiles['opt', method, solver, 400, n]
            ax.plot(a[:, 1], a[:, 2+column], label=f'MPI {n}', **styles[n])
        ax.set(xlabel='x', ylabel=labels[column], xlim=(.6, .9) if zoom else (0, 1))
        ax.grid(alpha=.2)
        if zoom:
            window = (baseline[:, 1] >= .6) & (baseline[:, 1] <= .9)
            lo, hi = np.min(baseline[window, 2+column]), np.max(baseline[window, 2+column])
            pad = max((hi-lo)*.08, .005)
            ax.set_ylim(lo-pad, hi+pad)

    for solver in ('hlle', 'hllc'):
        fig, axes = plt.subplots(3, 3, figsize=(13.5, 9), constrained_layout=True,
                                 gridspec_kw={'height_ratios': [2.2, 2.2, 1.2]})
        for col in range(3):
            panel(axes[0, col], 'weno5', solver, col)
            panel(axes[1, col], 'weno5', solver, col, zoom=True)
            baseline = profiles['opt', 'weno5', solver, 400, 1]
            maximum = 0.
            for n in (2, 4, 8):
                a = profiles['opt', 'weno5', solver, 400, n]
                difference = a[:, 2+col]-baseline[:, 2+col]
                maximum = max(maximum, float(np.max(np.abs(difference))))
                axes[2, col].plot(a[:, 1], difference, label=f'MPI {n} − MPI 1', **styles[n])
            axes[2, col].set(xlabel='x', ylabel=r'$\Delta$'+labels[col], xlim=(0, 1))
            axes[2, col].grid(alpha=.2)
            axes[2, col].text(.02, .85, f'máx |diferencia| = {maximum:.2e}',
                              transform=axes[2, col].transAxes, fontsize=9)
            if maximum == 0:
                axes[2, col].set(ylim=(-1, 1), yticks=[0])
            else:
                axes[2, col].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axes[0, 0].legend(fontsize=8, ncol=3)
        axes[2, 0].legend(fontsize=8, loc='lower left', ncol=3)
        axes[1, 0].set_title('Ampliación de la región detrás del choque')
        fig.suptitle(f'Shu–Osher relativista de GRIAL · WENO5–{solver.upper()} · N=400 · t=0.4\n'
                     'Perfiles, ampliación y diferencias entre 1, 2, 4 y 8 procesos MPI', fontsize=13)
        fig.savefig(out/f'shu_weno5_{solver}_mpi.png', dpi=200)
        fig.savefig(out/f'shu_weno5_{solver}_mpi.pdf')
        plt.close(fig)
        fig, axes = plt.subplots(5, 2, figsize=(12, 13), constrained_layout=True)
        for row, method in enumerate(METHODS):
            panel(axes[row, 0], method, solver, 0)
            panel(axes[row, 1], method, solver, 0, zoom=True)
            axes[row, 0].set_title(f'{method.upper()} · dominio completo')
            axes[row, 1].set_title(f'{method.upper()} · ampliación')
        axes[0, 0].legend(fontsize=8, ncol=3)
        fig.suptitle(f'Shu–Osher · {solver.upper()} · sensor activo · N=400 · t=0.4\n'
                     'Densidad con cinco reconstructores y MPI 1/2/4/8', fontsize=13)
        fig.savefig(out/f'shu_{solver}_reconstructores.png', dpi=180)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--serial-exe', type=Path, required=True)
    parser.add_argument('--mpi-exe', type=Path, required=True)
    parser.add_argument('--debug-serial-exe', type=Path, required=True)
    parser.add_argument('--debug-mpi-exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    profiles = {}; checks = []; variables = []; binary_hashes = {}
    for build in ('opt', 'debug'):
        serial = (args.serial_exe if build == 'opt' else args.debug_serial_exe).resolve()
        mpi = (args.mpi_exe if build == 'opt' else args.debug_mpi_exe).resolve()
        binary_hashes[build] = dict(serial=sha(serial), mpi=sha(mpi))
        configs = [(m, s, 400) for m in (METHODS if build == 'opt' else ('mp5', 'weno5'))
                   for s in (('hlle', 'hllc') if build == 'opt' else ('hlle',))]
        configs += [(m, 'hlle', 401) for m in ('mp5', 'weno5')]
        for method, solver, nx in configs:
            folder = out/'runs'/build/f'shu_{method}_{solver}_n{nx}'
            baseline, trace, _, _, _ = execute(serial, folder/'serial', 'shu', method, solver, nx, 'serial', 1, build == 'debug')
            baseline_initial = initial_profile(folder/'serial', 1, nx)
            profiles[build, method, solver, nx, 0] = baseline
            for n in PROCESSES:
                work = folder/f'mpi{n}'
                final, current_trace, counts, repair_delta, seconds = execute(mpi, work, 'shu', method, solver, nx, 'mpi', n, build == 'debug')
                initial = initial_profile(work, n, nx)
                assert current_trace.tobytes() == trace.tobytes(), 'different time sequence'
                assert np.array_equal(final[:, :2], baseline[:, :2]), 'different global coordinates'
                initial_equal = initial.tobytes() == baseline_initial.tobytes()
                final_equal = final.tobytes() == baseline.tobytes()
                if initial_equal:
                    assert final_equal, (build, method, solver, nx, n, 'equal inputs produced different outputs')
                initial_error = np.max(np.abs(initial[:, 2:]-baseline_initial[:, 2:]), axis=0)
                final_error = np.max(np.abs(final[:, 2:]-baseline[:, 2:]), axis=0)
                scaled_error = final_error/np.maximum(1., np.max(np.abs(baseline[:, 2:]), axis=0))
                assert np.sum(counts[8:14]) == 0 and np.all(repair_delta == 0), 'repaired Shu state'
                profiles[build, method, solver, nx, n] = final
                checks.append(dict(build=build, method=method, solver=solver, nx=nx, processes=n,
                    steps=len(trace), time=.4, initial_byte_identical=initial_equal, final_byte_identical=final_equal,
                    initial_max_abs_error=float(np.max(initial_error)), final_max_abs_error=float(np.max(final_error)),
                    final_max_scaled_error=float(np.max(scaled_error)), time_sequence_identical=True,
                    c2p_calls=int(counts[7]), c2p_repairs=int(np.sum(counts[8:14])),
                    events=counts.tolist(), seconds=seconds))
                for field, e0, e1, scaled in zip(FIELDS, initial_error, final_error, scaled_error):
                    variables.append(dict(build=build, method=method, solver=solver, nx=nx, processes=n,
                        variable=field, initial_Linf=float(e0), final_Linf=float(e1), final_scaled_Linf=float(scaled)))
            print(f'{build} Shu {method} {solver} N={nx}: MPI 1/2/4/8 compared; steps={len(trace)}', flush=True)
    write_table(out/'differences_by_field.tsv', variables)
    write_table(out/'comparisons.tsv', [{k: v for k, v in row.items() if not isinstance(v, list)} for row in checks])
    (out/'PARAMETERS.json').write_text(json.dumps(dict(domain=[0., 1.], interface=.5, time=.4, gamma=5/3,
        rho_left=5., p_left=50., rho_right='2+0.3*sin(50*x)', p_right=5., initial_velocity=[0., 0., 0.],
        cfl=.3, sensor=True, tvd_limiter='MC', halo=3, omp_threads=1,
        initialization='canonical production set_initial_conditions, centre samples for both grids',
        comparison='measure native local initial/final differences; exact evolution required for identical input',
        exact_solution_used=False, numerical_high_resolution_reference_used=False), indent=2)+'\n')
    make_plots(out, profiles)
    assert len(checks) == 64
    report = dict(passed=True, scope='Shu native initial/final MPI comparison; equal input must give exact output',
        trajectories=80, serial_references=16, mpi_comparisons=len(checks), checks=checks,
        processes=list(PROCESSES), nx=[400, 401], time=.4, binary_sha256=binary_hashes,
        initial_nonidentical_comparisons=sum(not r['initial_byte_identical'] for r in checks),
        final_nonidentical_comparisons=sum(not r['final_byte_identical'] for r in checks),
        initial_max_abs_error=max(r['initial_max_abs_error'] for r in checks),
        final_max_abs_error=max(r['final_max_abs_error'] for r in checks),
        final_max_scaled_error=max(r['final_max_scaled_error'] for r in checks),
        c2p_repairs=sum(r['c2p_repairs'] for r in checks), production_mpi_io=False,
        distributed_restart_tested=False, scalability_benchmark=False)
    (out/'VERIFIED.json').write_text(json.dumps(report, indent=2)+'\n')
    print('MPI SHU COMPARISON PASSED:', out)


if __name__ == '__main__':
    main()
