#!/usr/bin/env python3
"""Reject corrupt/incomplete MPI block collections and real comparison failures."""
import importlib.util
import json
import os
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

spec = importlib.util.spec_from_file_location('long_campaign', Path(__file__).resolve().parents[1]/'tools/srhd_long_campaign.py')
campaign = importlib.util.module_from_spec(spec);spec.loader.exec_module(campaign)


def block(path, offset=0, nx=4, step=0, current=0., value=1.):
    data = np.full((2, 5, 3, nx), value, dtype='=f8')
    path.write_bytes(campaign.HEADER.pack(b'GRIAL_LONG_V1   ', 4, 3, offset, nx, step, 5, current, .05 if step else 0.)+data.tobytes())


class TestLongCampaign(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory();self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def make_run(self, label, ranks):
        folder = self.root/label;folder.mkdir()
        for index in [0, 1]:
            for rank in range(ranks):
                block(folder/f'frame_{index:06d}_rank_{rank:04d}.bin', rank*(4//ranks), 4//ranks,
                      index*2, index*.1)
        np.savetxt(folder/'trace.dat', [[1, .05, .05, .01], [2, .1, .05, .02]])
        diagnostics = folder/'data/weno5_hlle';diagnostics.mkdir(parents=True)
        budget = np.zeros((2, 57));budget[:, 2:12] = 1
        np.savetxt(diagnostics/'mpi_conservation_budget.dat', budget)
        budget[:, 0] = [0, .1];budget[:, 1] = [0, 2]
        np.savetxt(diagnostics/'mpi_conservation_budget.dat', budget)
        state = np.zeros((2, 28));state[:, :2] = budget[:, :2];state[:, 5:7] = 1;state[:, 9] = .1
        np.savetxt(diagnostics/'mpi_global_diagnostics.dat', state)
        counters = np.zeros((2, 26));counters[:, :2] = budget[:, :2];counters[:, 9] = [0, 72]
        counters[:, 2] = ranks  # Face work may depend on the decomposition.
        np.savetxt(diagnostics/'mpi_numerical_diagnostics.dat', counters)
        return folder

    def test_complete_partition(self):
        run = self.make_run('run', 2)
        state = campaign.frame(run, 1, 2)
        self.assertEqual(state['u'].shape, (4, 3, 5))
        self.assertTrue(np.array_equal(state['u'], np.ones((4, 3, 5))))

    def test_missing_or_overlapping_rank(self):
        run = self.make_run('run', 2)
        (run/'frame_000001_rank_0001.bin').unlink()
        with self.assertRaises(ValueError):campaign.frame(run, 1, 2)
        block(run/'frame_000001_rank_0001.bin', 1, 2, 2, .1)
        with self.assertRaises(ValueError):campaign.frame(run, 1, 2)

    def test_corrupt_blocks(self):
        path = self.root/'block.bin';block(path)
        data = path.read_bytes()
        for broken in [data[:20], data[:-1], data+b'x', b'X'+data[1:]]:
            path.write_bytes(broken)
            with self.assertRaises(ValueError):campaign.read_block(path)
        block(path, value=float('nan'))
        with self.assertRaises(ValueError):campaign.read_block(path)

    def test_comparison_and_physical_counters(self):
        ref, run = self.make_run('ref', 1), self.make_run('run', 2)
        result = campaign.compare(ref, run, 2, 'weno5_hlle')
        self.assertTrue(result['passed']);self.assertTrue(result['fields_identical'])
        block(run/'frame_000001_rank_0001.bin', 2, 2, 2, .1, 1.01)
        self.assertFalse(campaign.compare(ref, run, 2, 'weno5_hlle')['passed'])

    def test_budget_or_recovery_failure(self):
        ref, run = self.make_run('ref', 1), self.make_run('run', 2)
        path = run/'data/weno5_hlle/mpi_conservation_budget.dat'
        budget = np.loadtxt(path);budget[-1, -1] = 1e-5;np.savetxt(path, budget)
        self.assertFalse(campaign.compare(ref, run, 2, 'weno5_hlle')['passed'])
        budget[-1, -1] = 0;np.savetxt(path, budget)
        path = run/'data/weno5_hlle/mpi_numerical_diagnostics.dat'
        counters = np.loadtxt(path);counters[-1, 10] = 1;np.savetxt(path, counters)
        self.assertFalse(campaign.compare(ref, run, 2, 'weno5_hlle')['passed'])

    def test_affinity_overlap(self):
        for rank, cpus in enumerate([[0, 1], [2, 3]]):
            (self.root/f'affinity_rank_{rank:04d}.json').write_text(json.dumps(dict(rank=rank, cpus=cpus, threads=2)))
        campaign.check_affinity(self.root, 2, 2, [0, 1, 2, 3])
        (self.root/'affinity_rank_0001.json').write_text(json.dumps(dict(rank=1, cpus=[1, 2], threads=2)))
        with self.assertRaises(ValueError):campaign.check_affinity(self.root, 2, 2, [0, 1, 2, 3])

    def test_nonfinite_or_invalid_diagnostic(self):
        ref, run = self.make_run('ref', 1), self.make_run('run', 2)
        path = run/'data/weno5_hlle/mpi_global_diagnostics.dat'
        values = np.loadtxt(path);values[-1, -1] = 1;np.savetxt(path, values)
        with self.assertRaises(ValueError):campaign.compare(ref, run, 2, 'weno5_hlle')
        values[-1, -1] = 0;values[-1, 2] = float('nan');np.savetxt(path, values)
        with self.assertRaises(ValueError):campaign.compare(ref, run, 2, 'weno5_hlle')

    def test_singleton_environment_preserves_allocation_and_libraries(self):
        affinity = os.sched_getaffinity(0)
        with patch.dict(os.environ, {'SLURM_JOB_ID':'123', 'PMI_RANK':'0',
                        'PMIX_NAMESPACE':'job', 'OMPI_COMM_WORLD_SIZE':'1',
                        'LD_LIBRARY_PATH':'/test/mpi/lib', 'OMP_NUM_THREADS':'18'}):
            env = campaign.singleton_test_environment()
            self.assertFalse(any(key.startswith(('SLURM_', 'PMI_', 'PMIX_', 'OMPI_COMM_WORLD_')) for key in env))
            self.assertEqual(env['LD_LIBRARY_PATH'], '/test/mpi/lib')
            self.assertEqual(env['OMP_NUM_THREADS'], '1')
            self.assertEqual(os.environ['SLURM_JOB_ID'], '123')
        self.assertEqual(affinity, os.sched_getaffinity(0))


if __name__ == '__main__':
    unittest.main()
