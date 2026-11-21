#!/usr/bin/env python3
"""Local C03 chain. Exit 0: target completed; 2: clean segment quota; 1: failure.

No scheduler submissions. MPI executable/grid and OpenMP threads are explicit.
Changing the grid on --resume is supported; changing code/physics is not.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import time
from run_segments import parameters, publish, stamp, sha

POLY = 0x42F0E1EBA9EA3693
MASK = (1 << 64)-1
TABLE = []
for index in range(256):
    value = index << 56
    for bit in range(8):
        value = ((value << 1) ^ (POLY if value & (1 << 63) else 0)) & MASK
    TABLE.append(value)


def crc64(data):
    value = 0
    for byte in data:
        value = TABLE[(value >> 56) ^ byte] ^ ((value << 8) & MASK)
    return value if value < (1 << 63) else value-(1 << 64)


def generation_hashes(directory):
    """Metadata CRC links plus streaming SHA256 of every immutable field block."""
    physical = (directory/'COMMITTED').read_bytes()
    observer = (directory/'OBSERVED').read_bytes()
    sidecar = (directory/'observer.bin').read_bytes()
    body = (directory/'manifest.bin').read_bytes()
    if len(body) < 132 or body[:16].rstrip() != b'GRIAL_MPI_CP_V1':
        raise ValueError('Invalid physical manifest')
    if not (len(physical) == 48 and physical[:16].rstrip() == b'GRIAL_CP_COMMIT1'
            and physical[-16:].rstrip() == b'GRIAL_CP_END1'
            and len(observer) == 64 and observer[:16].rstrip() == b'GRIAL_OBS_COM1'
            and observer[-16:].rstrip() == b'GRIAL_OBS_END1'):
        raise ValueError('Incomplete physical/observer commit')
    if (struct.unpack_from('=2q', physical, 16) != (len(body), crc64(body))
            or struct.unpack_from('=2q', observer, 16) != (len(physical), crc64(physical))
            or struct.unpack_from('=2q', observer, 32) != (len(sidecar), crc64(sidecar))):
        raise ValueError('Checkpoint and observer are not a coherent generation')
    ranks = struct.unpack_from('=i', body, 60)[0]
    if ranks < 1 or len(body) != 132+72*ranks:
        raise ValueError('Invalid checkpoint manifest')
    configuration = (directory/'configuration.bin').read_bytes()
    if struct.unpack_from('=2q', body, 100) != (len(configuration), crc64(configuration)):
        raise ValueError('Checkpoint configuration changed')
    expected = {'COMMITTED', 'OBSERVED', 'configuration.bin', 'manifest.bin', 'observer.bin'}
    for rank in range(ranks):
        descriptor = struct.unpack_from('=9q', body, 116+rank*72)
        name = f'block_{rank:09}.bin'
        if descriptor[0] != rank or (directory/name).stat().st_size != descriptor[7]:
            raise ValueError('Missing/truncated checkpoint block')
        expected.add(name)
    actual = {p.name for p in directory.iterdir()}
    if actual != expected:
        raise ValueError('Unexpected or incomplete generation files')
    result = {}
    for name in sorted(expected):
        digest = hashlib.sha256()
        with (directory/name).open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024*1024), b''):
                digest.update(chunk)
        result[name] = digest.hexdigest()
    return result


def validate(args):
    if any(not math.isfinite(x) for x in [args.segment_seconds, args.margin_seconds]):
        raise ValueError('Wall-clock controls must be finite')
    if args.segment_seconds < 0 or args.margin_seconds < 0:
        raise ValueError('Wall-clock controls cannot be negative')
    if args.segment_seconds > 0 and not 0 < args.margin_seconds < args.segment_seconds:
        raise ValueError('Require 0 < margin < segment seconds')
    if min(*args.grid, args.threads, args.max_segments) < 1 or args.test_segment_steps < 0:
        raise ValueError('Grid, threads and quota must be positive; test step quota nonnegative')
    if args.serial and args.grid != [1, 1]:
        raise ValueError('Non-MPI build requires --grid 1 1')
    if not args.exe.is_file() or not os.access(args.exe, os.X_OK) or not args.par.is_file():
        raise ValueError('Missing executable/parameter file')
    if not args.serial and shutil.which(args.mpirun) is None:
        raise ValueError('MPI launcher is unavailable')


def checkpoint_path(directory, value):
    path = (directory/value).resolve()
    path.relative_to(directory)
    return path


def run(args):
    args.exe = args.exe.resolve(); args.par = args.par.resolve(); directory = args.run_dir.resolve()
    validate(args)
    manifest_file = directory/'mpi_segments.json'
    if args.resume:
        manifest = json.loads(manifest_file.read_text())
        if manifest['executable_sha256'] != sha(args.exe) or manifest['parameter_sha256'] != sha(args.par):
            raise ValueError('Executable/parameters changed; use a separate run')
        if manifest['serial'] != args.serial:
            raise ValueError('Cannot change MPI build mode within a chain')
        if sha(directory/'case.par') != manifest['effective_parameter_sha256']:
            raise ValueError('Effective parameters were modified')
        if manifest['state'] not in ['completed', 'stopped_for_walltime']:
            raise ValueError('Previous session did not stop cleanly; inspect the logs')
        last = manifest['segments'][-1]
        checkpoint = checkpoint_path(directory, last['checkpoint'])
        if generation_hashes(checkpoint) != last['generation_sha256']:
            raise ValueError('Immutable checkpoint files changed')
        if manifest['state'] == 'completed':
            print('Physical target already completed:', directory)
            return 0
        mode = 'resume'; source = str(checkpoint)
        previous_step = last['step']; previous_time = last['time']
    else:
        text = args.par.read_text()
        for line in text.splitlines():
            item = line.split('#', 1)[0].strip().split('=', 1)
            if len(item) == 2 and item[0].strip().lower() == 'restart':
                if item[1].strip().lower() not in ['false', 'f', 'no', '0', '.false.']:
                    raise ValueError('Cold chain requires restart=false; use --resume for an existing chain')
        directory.mkdir(parents=True, exist_ok=False)
        text = parameters(text, {'output_folder': 'data', 'restart': 'false'})
        (directory/'case.par').write_text(text)
        manifest = dict(format='GRIAL_MPI_SEGMENTS_1', created_utc=stamp(), state='prepared',
                        executable_sha256=sha(args.exe), parameter_sha256=sha(args.par),
                        effective_parameter_sha256=sha(directory/'case.par'), serial=args.serial, segments=[])
        mode = 'cold'; source = str(directory/'case.par'); previous_step = 0; previous_time = 0.0
    lock = directory/'mpi_segments.lock'
    with lock.open('x') as stream:
        stream.write(str(os.getpid())+'\n')
    try:
        env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC='FALSE', OMP_WAIT_POLICY='PASSIVE')
        for unused in range(args.max_segments):
            number = len(manifest['segments'])+1
            log = directory/f'session_{number:06}.log'
            receipt = directory/f'session_{number:06}.dat'
            if log.exists() or receipt.exists() or receipt.with_suffix('.dat.tmp').exists():
                raise ValueError('Session files already exist; inspect the previous attempt')
            command = [] if args.serial else [args.mpirun, *args.launcher_arg, '-n', str(math.prod(args.grid))]
            command += [str(args.exe), mode, source, *map(str, args.grid), str(args.segment_seconds),
                        str(args.margin_seconds), str(args.test_segment_steps), str(number)]
            manifest.update(state='running', updated_utc=stamp())
            publish(manifest_file, manifest)
            start = time.monotonic()
            with log.open('xb') as stream:
                process = subprocess.run(command, cwd=directory, env=env, stdout=stream, stderr=subprocess.STDOUT)
            row = dict(number=number, command=command, grid=args.grid, threads=args.threads,
                       returncode=process.returncode, elapsed_launcher_seconds=time.monotonic()-start, log=log.name)
            try:
                if process.returncode != 0:
                    raise ValueError(f'MPI session failed with exit {process.returncode}')
                values = dict(line.split('=', 1) for line in receipt.read_text().splitlines())
                state = values['state']; step = int(values['step']); physical = float(values['time']); target = float(values['target'])
                if state not in ['completed', 'stopped_for_walltime'] or not all(map(math.isfinite, [physical, target])):
                    raise ValueError('Invalid segment receipt')
                if (int(values['ranks']) != math.prod(args.grid) or int(values['grid_r']) != args.grid[0]
                        or int(values['grid_theta']) != args.grid[1]):
                    raise ValueError('Receipt resource grid does not match the launch')
                if state == 'completed' and physical != target:
                    raise ValueError('False physical completion')
                if state == 'stopped_for_walltime' and not (previous_step < step and previous_time < physical < target):
                    raise ValueError('Stopped session made no scientific progress')
                checkpoint = checkpoint_path(directory, values['checkpoint'])
                checkpoint_epoch = struct.unpack_from('=2di', (checkpoint/'manifest.bin').read_bytes(), 76)
                if checkpoint_epoch[0] != physical or checkpoint_epoch[2] != step:
                    raise ValueError('Receipt and checkpoint epoch differ')
                row.update(state=state, step=step, time=physical, target=target, checkpoint=str(checkpoint.relative_to(directory)),
                           generation_sha256=generation_hashes(checkpoint), receipt=values)
            except (ValueError, OSError, KeyError) as error:
                row.update(state='failed', error=str(error));manifest['segments'].append(row)
                manifest.update(state='failed', updated_utc=stamp());publish(manifest_file, manifest)
                raise
            manifest['segments'].append(row);manifest.update(state=state, updated_utc=stamp())
            publish(manifest_file, manifest)
            print(f'Session {number}: {state}; step={step}; t={physical:.9g}/{target:.9g}', flush=True)
            if state == 'completed':
                return 0
            mode = 'resume';source = str(checkpoint);previous_step = step;previous_time = physical
        print('Segment quota reached with a valid checkpoint; resume explicitly.', flush=True)
        return 2
    finally:
        lock.unlink()


def parser():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--exe', type=Path, required=True)
    ap.add_argument('--par', type=Path, required=True)
    ap.add_argument('--run-dir', type=Path, required=True)
    ap.add_argument('--grid', nargs=2, type=int, default=[1, 1])
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--segment-seconds', type=float, default=0)
    ap.add_argument('--margin-seconds', type=float, default=0)
    ap.add_argument('--max-segments', type=int, default=1)
    ap.add_argument('--mpirun', default='mpirun')
    ap.add_argument('--launcher-arg', action='append', default=[])
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--serial', action='store_true', help='single-process build without MPI')
    ap.add_argument('--test-segment-steps', type=int, default=0, help='deterministic validation stops; zero disables')
    return ap


if __name__ == '__main__':
    try:
        sys.exit(run(parser().parse_args()))
    except (ValueError, OSError, KeyError) as error:
        print('MPI SEGMENT ERROR:', error, file=sys.stderr)
        sys.exit(1)
