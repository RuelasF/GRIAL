#!/usr/bin/env python3
"""Run GRIAL sequentially across wall-clock segments, preserving logs and receipts.

No scheduler submissions. Exit 0: physical target reached; 2: segment quota
reached with a valid restart; 1: failure or no progress. Resume explicitly.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def publish(path, data):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2)+'\n')
    temporary.replace(path)


def parameters(text, overrides):
    lines = []
    for line in text.splitlines():
        key = line.split('#', 1)[0].split('=', 1)[0].strip().lower()
        if key not in overrides:
            lines.append(line)
    return '\n'.join(lines)+'\n'+'\n'.join(f'{key}={value}' for key,value in overrides.items())+'\n'


def run(args):
    exe = args.exe.resolve()
    case = args.par.resolve()
    directory = args.run_dir.resolve()
    if not exe.is_file() or not os.access(exe, os.X_OK):
        raise ValueError('Executable is missing or not executable')
    for value in [args.segment_seconds, args.margin_seconds, args.checkpoint_seconds]:
        if not math.isfinite(value):
            raise ValueError('Wall-clock controls must be finite')
    if not 0 < args.margin_seconds < args.segment_seconds or args.checkpoint_seconds < 0:
        raise ValueError('Require 0 < margin < segment seconds; checkpoint seconds >= 0')
    if args.threads < 1 or args.max_segments < 1:
        raise ValueError('Threads and max-segments must be positive')
    manifest_path = directory/'segments.json'
    if args.resume:
        manifest = json.loads(manifest_path.read_text())
        if manifest['executable_sha256'] != sha(exe) or manifest['parameter_sha256'] != sha(case):
            raise ValueError('Executable or original parameters changed; start a separate run')
        if manifest['state'] == 'completed':
            print('Physical target already completed:', directory)
            return 0
        if manifest['state'] != 'stopped_for_walltime':
            raise ValueError('Previous run did not stop cleanly; inspect its logs before restarting')
        checkpoint = Path(manifest['segments'][-1]['checkpoint'])
        previous_step = manifest['segments'][-1]['step']
        previous_time = manifest['segments'][-1]['physical_time']
    else:
        directory.mkdir(parents=True, exist_ok=False)
        text = case.read_text()
        # A cold configuration avoids ambiguous original relative checkpoint paths.
        if re.search(r'^\s*restart\s*=\s*(true|\.true\.|t)\s*(#.*)?$',text,re.M|re.I):
            raise ValueError('Use a cold-start case; --resume continues this launcher run')
        manifest = dict(state='new', created=stamp(), executable=str(exe), executable_sha256=sha(exe),
                        parameters=str(case), parameter_sha256=sha(case), segments=[])
        (directory/'original.par').write_text(text)
        publish(manifest_path, manifest)
        checkpoint = None
        previous_step, previous_time = 0, 0.0
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC='FALSE')
    controls = dict(output_folder='data', walltime_limit_seconds=args.segment_seconds,
                    walltime_margin_seconds=args.margin_seconds, walltime_checkpoint_seconds=args.checkpoint_seconds)
    for _ in range(args.max_segments):
        number = len(manifest['segments'])+1
        folder = directory/f'segment_{number:04d}'
        folder.mkdir()
        entry = dict(number=number, started=stamp(), threads=args.threads, controls=controls,
                     restart_from=str(checkpoint) if checkpoint else None, outcome='running')
        manifest['segments'].append(entry)
        manifest['state'] = 'running'
        publish(manifest_path, manifest)
        try:
            overrides = dict(controls)
            if checkpoint:
                overrides.update(restart='true', restart_file=str(checkpoint))
                text = parameters('', overrides)
                check_file = folder/'check.par'
                check_file.write_text(text)
                with (folder/'check.log').open('wb') as log:
                    checked = subprocess.run([str(exe),'--check',str(check_file)],cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT)
                if checked.returncode:
                    raise ValueError('Restart checkpoint failed --check')
                found = re.search(r'CHECKPOINT_STATE step=(\d+) time=\s*([\d.Ee+\-]+)',(folder/'check.log').read_text())
                if not found or (int(found[1]),float(found[2])) != (previous_step, previous_time):
                    raise ValueError('Checkpoint disagrees with recorded step/time')
            else:
                text = parameters((directory/'original.par').read_text(),overrides)
            par = folder/'run.par'
            par.write_text(text)
            print(f'Segment {number}: start, OMP={args.threads}, restart={checkpoint}',flush=True)
            start = time.monotonic()
            with (folder/'run.log').open('wb') as log:
                p = subprocess.run([str(exe),str(par)],cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT)
            entry.update(returncode=p.returncode, elapsed_seconds=time.monotonic()-start, ended=stamp())
            if p.returncode:
                raise ValueError(f'Solver exited with code {p.returncode}; see {folder / "run.log"}')
            status = dict(line.split('=',1) for line in (directory/'data/run_status.dat').read_text().splitlines())
            state = status['status']
            step, physical_time, target = int(status['step']),float(status['physical_time']),float(status['target_time'])
            checkpoint = (directory/status['checkpoint']).resolve()
            if state not in ['completed','stopped_for_walltime'] or not all(map(math.isfinite,[physical_time,target])):
                raise ValueError('Invalid terminal status')
            if checkpoint.suffix != '.rst' or not checkpoint.is_file() or (directory/'data') not in checkpoint.parents:
                raise ValueError('Missing or unexpected checkpoint path')
            if step < previous_step or physical_time < previous_time:
                raise ValueError('Evolution went backwards')
            if state == 'completed' and physical_time < target:
                raise ValueError('Completed status precedes physical target')
            if state == 'stopped_for_walltime' and (step <= previous_step or physical_time <= previous_time or physical_time >= target):
                raise ValueError('No physical progress, or inconsistent stop; increase segment budget')
            # Validate the new generation before recording any continuation.
            verify_file=folder/'verify.par'
            verify_file.write_text(parameters('',dict(controls,restart='true',restart_file=str(checkpoint))))
            with (folder/'verify.log').open('wb') as log:
                checked=subprocess.run([str(exe),'--check',str(verify_file)],cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT)
            found=re.search(r'CHECKPOINT_STATE step=(\d+) time=\s*([\d.Ee+\-]+)',(folder/'verify.log').read_text())
            if checked.returncode or not found or (int(found[1]),float(found[2]))!=(step,physical_time):
                raise ValueError('New checkpoint failed validation or disagrees with terminal status')
            entry.update(outcome=state,step=step,physical_time=physical_time,target_time=target,checkpoint=str(checkpoint))
            manifest['state']=state
            publish(folder/'receipt.json',entry)
            publish(manifest_path,manifest)
            print(f'Segment {number}: {state}, step={step}, t={physical_time:.16g}',flush=True)
            previous_step,previous_time=step,physical_time
            if state=='completed':
                return 0
        except (OSError,ValueError,KeyError) as error:
            entry.update(outcome='failed',error=str(error),ended=stamp())
            manifest['state']='failed'
            publish(folder/'receipt.json',entry)
            publish(manifest_path,manifest)
            raise
    print('Segment quota reached. Checkpoint retained; use --resume:',directory,flush=True)
    return 2


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('par',type=Path)
    ap.add_argument('--exe',type=Path,default=Path(__file__).resolve().parents[1]/'grial')
    ap.add_argument('--run-dir',type=Path,required=True)
    ap.add_argument('--segment-seconds',type=float,default=86400)
    ap.add_argument('--margin-seconds',type=float,default=900)
    ap.add_argument('--checkpoint-seconds',type=float,default=1800)
    ap.add_argument('--threads',type=int,default=4)
    ap.add_argument('--max-segments',type=int,default=3)
    ap.add_argument('--resume',action='store_true')
    args=ap.parse_args()
    try:
        return run(args)
    except (OSError,ValueError,KeyError) as error:
        print('ERROR:',error,file=sys.stderr)
        return 1


if __name__=='__main__':
    sys.exit(main())
