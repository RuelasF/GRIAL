#!/usr/bin/env python3
"""Launcher state/failure checks plus native partial-write and history controls."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('segments',ROOT/'tools/run_segments.py')
launcher=importlib.util.module_from_spec(spec);spec.loader.exec_module(launcher)
OUT=None


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.folder=Path(tempfile.mkdtemp(prefix='mock_',dir=OUT))
        self.exe=self.folder/'solver';self.exe.write_text('test executable\n');self.exe.chmod(0o755)
        self.par=self.folder/'case.par';self.par.write_text('problem=sod\nfinal_time=30\n')
        self.args=SimpleNamespace(exe=self.exe,par=self.par,run_dir=self.folder/'run',resume=False,
                                 segment_seconds=75.,margin_seconds=5.,checkpoint_seconds=25.,threads=4,max_segments=3)
        self.evolutions=0;self.mode='normal'

    def fake_run(self,cmd,*,cwd,env,stdout,stderr):
        cwd=Path(cwd)
        par=Path(cmd[-1]);text=par.read_text()
        if '--check' in cmd:
            cp=Path(next(line.split('=',1)[1] for line in text.splitlines() if line.startswith('restart_file=')))
            if self.mode=='bad_checkpoint':return SimpleNamespace(returncode=8)
            stored=json.loads(cp.read_text());stdout.write(f"CHECKPOINT_STATE step={stored['step']} time={stored['time']:.16e}\n".encode())
            return SimpleNamespace(returncode=0)
        self.evolutions+=1
        if self.mode=='bad_exit':return SimpleNamespace(returncode=7)
        if self.mode=='stale_status':return SimpleNamespace(returncode=0)
        data=cwd/'data';data.mkdir(exist_ok=True)
        step=10*self.evolutions;time_value=float(step)
        if self.mode=='no_progress':step=0;time_value=0.
        state='completed' if self.evolutions==3 else 'stopped_for_walltime'
        if self.mode=='false_completed':state='completed'
        cp=data/f'checkpoint_{step:09d}.rst';cp.write_text(json.dumps({'step':step,'time':time_value}))
        (data/'run_status.dat').write_text(f'status={state}\nstep={step}\nphysical_time={time_value}\ntarget_time=30\ncheckpoint={cp}\n')
        return SimpleNamespace(returncode=0)

    def test_three_segments_stop_at_completion(self):
        self.args.max_segments=6
        with patch.object(launcher.subprocess,'run',self.fake_run):
            self.assertEqual(launcher.run(self.args),0)
        self.assertEqual(self.evolutions,3)
        self.assertFalse((self.args.run_dir/'segment_0004').exists())
        manifest=json.loads((self.args.run_dir/'segments.json').read_text())
        self.assertEqual([s['outcome'] for s in manifest['segments']],['stopped_for_walltime']*2+['completed'])

    def test_quota_resume_and_completed_noop(self):
        self.args.max_segments=1
        with patch.object(launcher.subprocess,'run',self.fake_run):
            self.assertEqual(launcher.run(self.args),2)
            self.args.resume=True;self.args.max_segments=3
            self.assertEqual(launcher.run(self.args),0)
            self.assertEqual(launcher.run(self.args),0)
        self.assertEqual(self.evolutions,3)

    def test_failures_do_not_spawn_continuations(self):
        for mode in ['bad_exit','bad_checkpoint','stale_status','no_progress','false_completed']:
            with self.subTest(mode=mode):
                self.mode=mode;self.evolutions=0;self.args.run_dir=self.folder/mode
                with patch.object(launcher.subprocess,'run',self.fake_run), self.assertRaises((OSError,ValueError)):
                    launcher.run(self.args)
                self.assertEqual(self.evolutions,1)
                self.assertEqual(json.loads((self.args.run_dir/'segments.json').read_text())['state'],'failed')
                self.args.resume=True
                with self.assertRaises(ValueError):launcher.run(self.args)
                self.args.resume=False

    def test_changed_executable_is_rejected(self):
        self.args.max_segments=1
        with patch.object(launcher.subprocess,'run',self.fake_run):self.assertEqual(launcher.run(self.args),2)
        self.args.resume=True;self.exe.write_text('changed\n')
        with self.assertRaises(ValueError):launcher.run(self.args)

    def test_existing_directory_is_preserved(self):
        self.args.run_dir.mkdir();marker=self.args.run_dir/'keep';marker.write_text('keep\n')
        with self.assertRaises(FileExistsError):launcher.run(self.args)
        self.assertEqual(marker.read_text(),'keep\n')


def native_controls(exe,probe,out):
    out.mkdir()
    par=out/'case.par'
    par.write_text('''problem=sod
nx=800
final_time=0.02
save_interval=1
output_folder=data
numerical_diagnostics=true
diagnostic_stride=5
reconstruction=weno5
shock_sensor=true
''')
    env=dict(os.environ,OMP_NUM_THREADS='2',OMP_DYNAMIC='FALSE')
    def invoke(command,name,**kwargs):
        p=subprocess.run(command,cwd=out,env=env,capture_output=True,**kwargs)
        (out/(name+'.log')).write_bytes(p.stdout+p.stderr)
        return p
    assert invoke([str(exe),str(par)],'cold').returncode==0
    status=dict(line.split('=',1) for line in (out/'data/run_status.dat').read_text().splitlines())
    cp=out/status['checkpoint'];before=hashlib.sha256(cp.read_bytes()).hexdigest()
    def limit_file_size(size):
        signal.signal(signal.SIGXFSZ,signal.SIG_IGN)
        resource.setrlimit(resource.RLIMIT_CORE,(0,0))
        resource.setrlimit(resource.RLIMIT_FSIZE,(size,size))
    partial=Path(str(cp)+'.tmp')
    for limit in (4096,cp.stat().st_size-1):
        failed=invoke([str(probe),'--rewrite-failure',str(cp)],f'limited_write_{limit}',
                      preexec_fn=lambda: limit_file_size(limit))
        assert failed.returncode==0,(failed.stdout,failed.stderr)
        assert b'CHECKPOINT WRITE FAILURE HANDLED' in failed.stdout
        assert hashlib.sha256(cp.read_bytes()).hexdigest()==before
        assert partial.exists() and 0<partial.stat().st_size<=limit
    restart=out/'restart.par';restart.write_text(f'restart=true\nrestart_file={cp}\n')
    assert invoke([str(exe),'--check',str(restart)],'valid_generation').returncode==0
    bad=out/'partial.par';bad.write_text(f'restart=true\nrestart_file={partial}\n')
    assert invoke([str(exe),'--check',str(bad)],'partial_rejected').returncode!=0
    # Later samples are evidence: reject rollback rather than silently trim/append.
    history=out/'data/weno5_hlle/numerical_diagnostics.dat'
    with history.open('a') as f:f.write('1.0 1000 0\n')
    history_hash=hashlib.sha256(history.read_bytes()).hexdigest()
    rollback=invoke([str(exe),str(restart)],'history_ahead')
    assert rollback.returncode!=0 and b'history is ahead of checkpoint' in rollback.stdout
    assert hashlib.sha256(history.read_bytes()).hexdigest()==history_hash
    (out/'VERIFIED.json').write_text(json.dumps(dict(passed=True,original_generation_sha256=before,
        retained_generation_sha256=hashlib.sha256(cp.read_bytes()).hexdigest(),partial_bytes=partial.stat().st_size,
        history_unchanged=True,write_failure_returncode=failed.returncode),indent=2)+'\n')


def main():
    global OUT
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--exe',type=Path,default=ROOT/'grial')
    ap.add_argument('--probe',type=Path,default=ROOT/'tests/test_checkpoint_v5')
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    OUT=args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-segments-'))
    OUT.mkdir(parents=True,exist_ok=True)
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(LauncherTests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():return 1
    native_controls(args.exe.resolve(),args.probe.resolve(),OUT/'native')
    print('SEGMENT LAUNCHER TESTS PASSED:',OUT)
    return 0


if __name__=='__main__':sys.exit(main())
