#!/usr/bin/env python3
"""Distributed checkpoint round trips, independent checksums and failed publication."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import struct
import subprocess
import tempfile

POLYNOMIAL = 0x42F0E1EBA9EA3693
MASK = (1 << 64) - 1
TABLE = []
for byte in range(256):
    value = byte << 56
    for _ in range(8):
        value = ((value << 1) ^ (POLYNOMIAL if value & (1 << 63) else 0)) & MASK
    TABLE.append(value)

def crc64(data):
    value = 0
    for byte in data:
        value = TABLE[(value >> 56) ^ byte] ^ ((value << 8) & MASK)
    return value if value < 1 << 63 else value - (1 << 64)

assert crc64(b'123456789') == 0x6C40DF5F0B497347

def hashes(path):
    return {str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(path.rglob('*')) if p.is_file()}

def inspect(work, count):
    ranks = []
    for rank in range(count):
        row = (work/f'checkpoint_rank_{rank}.txt').read_text().split()
        assert list(map(int, row[:2])) == [rank, count] and int(row[3]) == 34
        errors = list(map(float, row[4:]))
        assert errors[0] == 0 and all(math.isfinite(x) and x < 1e-11 for x in errors)
        ranks.append(dict(rank=rank, cases=34, field_error=errors[0], budget_error=errors[1]))
    generations = []
    for case in sorted(work.glob('case_*')):
        for generation in ['generation_a', 'generation_b']:
            directory = case/generation
            commit = (directory/'COMMITTED').read_bytes()
            assert len(commit) == 48 and commit[:16].rstrip() == b'GRIAL_CP_COMMIT1'
            assert commit[-16:].rstrip() == b'GRIAL_CP_END1'
            body = (directory/'manifest.bin').read_bytes()
            assert struct.unpack_from('=2q', commit, 16) == (len(body), crc64(body))
            assert body[:16].rstrip() == b'GRIAL_MPI_CP_V1'
            assert struct.unpack_from('=5i', body, 16) == (32, 64, 64, 32, 0x01020304)
            assert struct.unpack_from('=i', body, 60)[0] == count
            assert len(body) == 132 + 72*count
            configuration = (directory/'configuration.bin').read_bytes()
            assert configuration[:12].rstrip() == b'GRHDCP_V9'
            assert struct.unpack_from('=2q', body, 100) == (len(configuration), crc64(configuration))
            descriptors = []
            for rank in range(count):
                descriptor = struct.unpack_from('=9q', body, 116+rank*72)
                assert descriptor[0] == rank and all(x > 0 for x in descriptor[1:4])
                data = (directory/f'block_{rank:09}.bin').read_bytes()
                assert descriptor[7:] == (len(data), crc64(data))
                assert data[:16].rstrip() == b'GRIAL_CP_BLOCK1' and data[-16:].rstrip() == b'GRIAL_CP_END1'
                assert struct.unpack_from('=2i', data, 36) == (rank, count)
                assert struct.unpack_from('=3i', data, 56) == descriptor[1:4]
                assert struct.unpack_from('=3i', data, 68) == descriptor[4:7]
                assert struct.unpack_from('=2di', data, 100) == struct.unpack_from('=2di', body, 76)
                cells = math.prod(descriptor[1:4])
                assert len(data) == 292 + 80*cells + 16
                counters = struct.unpack_from('=16q', data, 124)
                assert all(x > 2**32 for x in counters)
                descriptors.append(list(descriptor))
            assert len(list(directory.iterdir())) == count+3
            generations.append(dict(path=str(directory.relative_to(work)), descriptors=descriptors))
    assert len(generations) == 68
    return ranks, generations

def repair_manifest(directory, rank=None):
    body = bytearray((directory/'manifest.bin').read_bytes())
    if rank is not None:
        data = (directory/f'block_{rank:09}.bin').read_bytes()
        struct.pack_into('=2q', body, 116+rank*72+56, len(data), crc64(data))
    (directory/'manifest.bin').write_bytes(body)
    commit = bytearray((directory/'COMMITTED').read_bytes())
    struct.pack_into('=2q', commit, 16, len(body), crc64(body))
    (directory/'COMMITTED').write_bytes(commit)

def corrupt(source, destination, kind, rank):
    shutil.copytree(source, destination)
    block = destination/f'block_{rank:09}.bin'
    if kind == 'missing-commit':
        (destination/'COMMITTED').rename(destination/'COMMITTED.hidden')
    elif kind == 'missing-block':
        block.rename(block.with_suffix('.hidden'))
    elif kind == 'manifest':
        data = bytearray((destination/'manifest.bin').read_bytes());data[40] ^= 1
        (destination/'manifest.bin').write_bytes(data)
    elif kind == 'configuration':
        data = bytearray((destination/'configuration.bin').read_bytes());data[80] ^= 1
        (destination/'configuration.bin').write_bytes(data)
    elif kind == 'abi':
        data = bytearray((destination/'manifest.bin').read_bytes());struct.pack_into('=i', data, 16, 16)
        (destination/'manifest.bin').write_bytes(data);repair_manifest(destination)
    else:
        data = bytearray(block.read_bytes())
        if kind == 'truncate':
            data = data[:4096]
        elif kind == 'nan-with-valid-crc':
            shape = struct.unpack_from('=3i', data, 56)
            struct.pack_into('=d', data, 292+40*math.prod(shape), math.nan)
        elif kind == 'wrong-rank-with-valid-crc':
            struct.pack_into('=i', data, 36, rank+1)
        else:
            data[300] ^= 1
        block.write_bytes(data)
        if kind.endswith('valid-crc'):
            repair_manifest(destination, rank)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--grid', nargs=2, type=int)
    args = parser.parse_args()
    probe = args.probe.resolve()
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-c01-'))
    out.mkdir(parents=True, exist_ok=True)
    assert not (out/'PROGRESS.json').exists(), 'Use a new output directory'
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE', OMP_WAIT_POLICY='PASSIVE')
    env.pop('OMPI_MPI_THREAD_LEVEL', None)
    results = []
    def run(work, grid, mode, target='', success=True):
        command = ['mpirun', '--oversubscribe', '--bind-to', 'none', '-np', str(math.prod(grid))]
        # Fault modes set the last rank's limit inside the probe, after MPI initialization.
        command += [str(probe), *map(str, grid), mode, str(target)]
        log_path = work/(mode+'_'+(Path(target).name if target else 'main')+'.log')
        with log_path.open('wb') as log:
            process = subprocess.Popen(command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = process.wait(timeout=600 if mode=='report' else 60)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL);process.wait()
                raise RuntimeError('C01 probe timed out: '+str(log_path))
        output = log_path.read_text(errors='replace')
        assert (code == 0) == success, (log_path, code, output[-4000:])
        return code, output

    grids = [tuple(args.grid)] if args.grid else [(1,1),(2,1),(1,2),(4,1),(1,4),(2,2),(4,2),(2,4)]
    fault_grids = [(1,1),(2,2),(2,4)] if not args.grid else []
    corruption_kinds = ['missing-commit','missing-block','manifest','configuration','truncate','block',
                        'abi','nan-with-valid-crc','wrong-rank-with-valid-crc']
    for grid in grids:
        work = out/f'checkpoint_{grid[0]}x{grid[1]}';work.mkdir()
        run(work, grid, 'report')
        ranks, generations = inspect(work, math.prod(grid))
        result = dict(grid=grid, passed=True, ranks=ranks, generations=generations, corruption_checks=[])
        if grid in fault_grids:
            older = work/'case_1_1/generation_a';previous = hashes(older)
            for kind in corruption_kinds:
                broken = work/('broken_'+kind)
                corrupt(work/'case_1_1/generation_b', broken, kind, math.prod(grid)-1)
                _, output = run(work, grid, 'read-failure', broken)
                assert 'CHECKPOINT RECOVERY PROBE PASSED' in output
                result['corruption_checks'].append(kind)
            for fault, mode in [('returned','write-failure'),('crash','write-crash')]:
                failed = work/('failed_'+fault)
                code, output = run(work, grid, mode, failed, success=fault=='returned')
                if fault=='returned':
                    assert 'CHECKPOINT IO FAILURE PROBE PASSED' in output
                    rank = math.prod(grid)-1
                    block = failed/f'block_{rank:09}.bin'
                    assert re.search(
                        rf'CHECKPOINT I/O rank={rank} (write|close|verify) {re.escape(str(block))} iostat=-?\d+',
                        output,
                    ), output
                else:
                    assert 'SIGXFSZ' in output or 'signal 25' in output
                assert failed.exists() and not (failed/'COMMITTED').exists()
                partial = failed/f'block_{math.prod(grid)-1:09}.bin'
                assert partial.exists() and partial.stat().st_size == 4096
                assert hashes(older) == previous
                result.setdefault('write_failures', []).append(dict(mode=fault, returncode=code, bytes=4096, old_generation_unchanged=True))
            run(work, grid, 'read-success', older)
        results.append(result)
        (out/'PROGRESS.json').write_text(json.dumps(results,indent=2)+'\n')
        print('checkpoint',grid,'PASS',flush=True)

    invalid = {'bad-layout':'inconsistent domain layout','bad-path':'checkpoint path differs',
               'bad-clock':'inconsistent diagnostic time epoch','bad-seed':'checkpoint configuration or epoch differs',
               'bad-flag':'checkpoint configuration or epoch differs','bad-state':'invalid physical checkpoint state',
               'bad-counters':'invalid checkpoint counters','bad-eos':'invalid checkpoint epoch or equation of state',
               'bad-midstep':'cannot checkpoint inside an RK step',
               'bad-load-active':'checkpoint observer restoration requires C03','bad-no-halos':'checkpoint recovery requires prepared halos'}
    if not args.grid:
        for mode, marker in invalid.items():
            for grid in [(1,2),(2,2)]:
                work = out/f'{mode}_{grid[0]}x{grid[1]}';work.mkdir()
                code, output = run(work, grid, mode, success=False)
                assert marker in output, (mode, output[-4000:])
                assert not list(work.glob('invalid_generation/*')), 'Invalid request wrote checkpoint data'
                results.append(dict(grid=grid, mode=mode, rejected=True, returncode=code))
    (out/'VERIFIED.json').write_text(json.dumps(dict(passed=True, results=results, threshold=1e-11,
        scope='Same prepared partition; physical state/counters only. No C02 redistribution, C03 observer continuity, production, multi-node or scaling.'),indent=2)+'\n')
    print('CHECKPOINT MPI MATRIX PASSED:',out)

if __name__ == '__main__':
    main()
