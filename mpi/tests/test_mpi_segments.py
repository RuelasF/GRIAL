"""Small controls for the MPI chain launcher, without MPI or numerical simulations."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
import sys

TOOLS=Path(__file__).resolve().parents[1]/'tools'
sys.path.insert(0,str(TOOLS))
spec=importlib.util.spec_from_file_location('mpi_launcher',TOOLS/'run_mpi_segments.py')
launcher=importlib.util.module_from_spec(spec);spec.loader.exec_module(launcher)

class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory(prefix='grial-c03-launcher-')
        self.root=Path(self.temporary.name)
        self.exe=self.root/'fake';self.exe.write_text('#!/usr/bin/env python3\nraise SystemExit(7)\n');self.exe.chmod(0o700)
        self.par=self.root/'case.par';self.par.write_text('problem=convergence\n')
        self.args=SimpleNamespace(exe=self.exe,par=self.par,run_dir=self.root/'run',grid=[1,1],threads=1,
            max_segments=1,test_segment_steps=0,segment_seconds=0,margin_seconds=0,serial=True,
            mpirun='mpirun',launcher_arg=[],resume=False)
    def tearDown(self):self.temporary.cleanup()
    def test_crc_oracle(self):self.assertEqual(launcher.crc64(b'123456789'),0x6C40DF5F0B497347)
    def test_valid_controls(self):launcher.validate(self.args)
    def test_nonfinite(self):
        for value in [float('nan'),float('inf'),-float('inf')]:
            self.args.segment_seconds=value
            with self.assertRaises(ValueError):launcher.validate(self.args)
    def test_bad_margin(self):
        self.args.segment_seconds=10
        for margin in [-1,0,10,11,float('nan')]:
            self.args.margin_seconds=margin
            with self.assertRaises(ValueError):launcher.validate(self.args)
    def test_bad_counts(self):
        for key,value in [('grid',[0,1]),('threads',0),('max_segments',0),('test_segment_steps',-1)]:
            old=getattr(self.args,key);setattr(self.args,key,value)
            with self.assertRaises(ValueError):launcher.validate(self.args)
            setattr(self.args,key,old)
    def test_serial_grid(self):
        self.args.grid=[2,1]
        with self.assertRaises(ValueError):launcher.validate(self.args)
    def test_path_escape(self):
        with self.assertRaises(ValueError):launcher.checkpoint_path(self.root,'../escape')
    def test_missing_generation(self):
        with self.assertRaises(OSError):launcher.generation_hashes(self.root)
    def test_cold_restart_rejected(self):
        self.par.write_text('problem=convergence\nrestart=true\n')
        with self.assertRaises(ValueError):launcher.run(self.args)
        self.assertFalse(self.args.run_dir.exists())
    def test_child_failure_is_persistent(self):
        with self.assertRaisesRegex(ValueError,'failed with exit 7'):launcher.run(self.args)
        manifest=json.loads((self.args.run_dir/'mpi_segments.json').read_text())
        self.assertEqual(manifest['state'],'failed');self.assertEqual(len(manifest['segments']),1)
        self.assertFalse((self.args.run_dir/'mpi_segments.lock').exists())
        self.args.resume=True
        with self.assertRaisesRegex(ValueError,'did not stop cleanly'):launcher.run(self.args)
        self.assertFalse((self.args.run_dir/'session_000002.log').exists())
    def test_false_completion_rejected(self):
        self.exe.write_text('''#!/usr/bin/env python3
from pathlib import Path
Path('session_000001.dat').write_text('state=completed\\ntime=1\\ntarget=2\\nstep=1\\nranks=1\\ngrid_r=1\\ngrid_theta=1\\n')
''')
        with self.assertRaisesRegex(ValueError,'False physical completion'):launcher.run(self.args)
        self.assertEqual(json.loads((self.args.run_dir/'mpi_segments.json').read_text())['state'],'failed')

if __name__=='__main__':unittest.main()
