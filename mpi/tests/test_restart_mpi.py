#!/usr/bin/env python3
"""Cross-launch restart chains, polar redistribution and V7/V8/V9 import."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import struct
import subprocess
import tempfile
from test_checkpoint_mpi import corrupt,crc64,repair_manifest
from test_walltime_restart import records

CASES=[(fixture,config) for fixture in range(1,14) for config in range(1,4 if fixture>=6 else 3)]

def history(path):
    data=path.read_text().splitlines()
    return [list(map(int,data[0].split())),list(map(float,data[1].split())),
            list(map(int,data[2].split())),list(map(float,data[3].split()))]

def convert_legacy(source,destination,version):
    r=records(source)
    if version==7:
        r=[b'GRHDCP_V7   ']+r[1:24]+r[-4:-1]
    elif version==8:
        r=[b'GRHDCP_V8   ']+r[1:27]+r[-4:-1]+[r[-1].replace(b'GRIAL_END_V9',b'GRIAL_END_V8',1)]
    destination.write_bytes(b''.join(struct.pack('=i',len(x))+x+struct.pack('=i',len(x)) for x in r))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--probe',type=Path,required=True)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--serial',action='store_true')
    ap.add_argument('--smoke',action='store_true')
    args=ap.parse_args()
    probe=args.probe.resolve()
    out=args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-c02-'))
    out.mkdir(parents=True,exist_ok=True)
    assert not (out/'PROGRESS.json').exists(),'Use a new output directory'
    env=dict(os.environ,OMP_NUM_THREADS='2',OMP_DYNAMIC='FALSE',OMP_WAIT_POLICY='PASSIVE')
    env.pop('OMPI_MPI_THREAD_LEVEL',None)
    results=[]
    def hashes(paths):
        return {str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest()
                for path in paths for p in ([path] if path.is_file() else sorted(path.rglob('*'))) if p.is_file()}
    def invoke(work,grid,mode,source='',destination='',target=8,reference='ref',success=True):
        command=[] if args.serial else ['mpirun','--oversubscribe','--bind-to','none','-n',str(math.prod(grid))]
        command += [str(probe),*map(str,grid),mode,str(source),str(destination),str(target),reference]
        log=work/(mode+'_'+str(destination or source).replace('/','_')+'.log')
        with log.open('wb') as stream:
            p=subprocess.Popen(command,cwd=work,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            try:code=p.wait(timeout=600)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid,signal.SIGKILL);p.wait();raise RuntimeError('C02 timeout: '+str(log))
        output=log.read_text(errors='replace')
        assert (code==0)==success,(log,code,output[-4500:])
        return code,output
    def prepare_outputs(label):
        for f,c in CASES:
            for scheme in ['c01','weno5_hlle','mp5_hllc']:
                (out/f'case_{f}_{c}'/('events_'+label)/scheme).mkdir(parents=True,exist_ok=True)
    def run(grid,mode,source='',destination='',target=8,reference='ref'):
        label=reference if mode in ['produce','serial-reference'] else destination
        prepare_outputs(label)
        inputs=[out/f'case_{f}_{c}'/source for f,c in CASES] if source else []
        before_hashes=hashes(inputs)
        _,output=invoke(out,grid,mode,source,destination,target,reference)
        assert hashes(inputs)==before_hashes,'Restart edited its input generation/file'
        assert 'PARALLEL RESTART TESTS PASSED' in output
        ranks=[]
        for rank in range(math.prod(grid)):
            path=out/f'restart_{mode}_rank_{rank}.txt'
            row=path.read_text().split()
            assert list(map(int,row[:2]))==[rank,math.prod(grid)] and int(row[3])==34
            field,budget=map(float,row[4:]);assert field==0 and math.isfinite(budget) and budget<1e-11
            ranks.append(dict(rank=rank,cases=34,field_error=field,budget_error=budget))
            path.rename(out/f'{destination or reference}_{mode}_rank_{rank}.txt')
        if mode=='resume':
            for f,c in CASES:
                case=out/f'case_{f}_{c}'
                before=history(case/('history_'+source+'.txt'))[2:]
                actual=history(case/('history_'+destination+'.txt'))[:2]
                assert before[0]==actual[0],(f,c,source,destination,'historical counts changed')
                assert max(abs(a-b) for a,b in zip(before[1],actual[1]))<1e-11
        results.append(dict(grid=grid,mode=mode,source=source,destination=destination,target=target,reference=reference,
                            ranks=ranks,input_files_unchanged=True))
        (out/'PROGRESS.json').write_text(json.dumps(results,indent=2)+'\n')
        print(mode,grid,source,'->',destination or reference,'PASS',flush=True)

    run((1,1),'produce')
    if args.serial or args.smoke:
        run((1,1),'resume','ref_2','single',8)
    else:
        for grid,source,destination,target in [((2,1),'ref_2','chain2',4),((4,1),'chain2','chain4',6),
            ((2,1),'chain4','chain2_final',8),((2,2),'chain2','polar',8),((1,4),'chain4','theta',8),
            ((4,2),'ref_2','eight',8),((2,4),'ref_2','eight_early',4),((1,1),'eight_early','back_one',8)]:
            run(grid,'resume',source,destination,target)
    versions=[9] if args.smoke else [7,8,9]
    for version in versions:
        for f,c in CASES:
            case=out/f'case_{f}_{c}'
            convert_legacy(case/'events_ref/checkpoint_c01_step_000000002.rst',case/f'legacy_v{version}.rst',version)
        prefix=f'legacy_ref{version}'
        run((1,1),'serial-reference',f'legacy_v{version}.rst',reference=prefix)
        profiles=[(1,1)] if args.serial or args.smoke else [(1,1),(2,1),(2,2),(2,4)]
        for grid in profiles:
            target=2 if math.prod(grid)==8 else 8
            dest=f'import{version}_{grid[0]}x{grid[1]}'
            run(grid,'import',f'legacy_v{version}.rst',dest,target,prefix)
        if not args.serial and not args.smoke:
            run((2,1),'resume',f'import{version}_2x4',f'import{version}_repartitioned',8,prefix)
    rejection_results=[]
    if not args.serial and not args.smoke:
        original=out/'case_1_1/chain4'
        original_hash=hashes([original])
        for kind in ['missing-commit','manifest','configuration','missing-block','truncate','block',
                     'nan-with-valid-crc','wrong-rank-with-valid-crc','counter-overflow','settings','metadata-short','bad-gamma']:
            work=out/('reject_'+kind);work.mkdir();bad=work/'source'
            if kind in ['counter-overflow','settings','metadata-short','bad-gamma']:
                shutil.copytree(original,bad)
                if kind=='counter-overflow':
                    for rank in range(4):
                        path=bad/f'block_{rank:09}.bin';data=bytearray(path.read_bytes())
                        struct.pack_into('=q',data,124,(1<<63)-1);path.write_bytes(data);repair_manifest(bad,rank)
                elif kind in ['metadata-short','bad-gamma']:
                    config=bytearray((bad/'configuration.bin').read_bytes())
                    if kind=='metadata-short':config=config[:20]
                    else:struct.pack_into('=d',config,84,-1.0)
                    (bad/'configuration.bin').write_bytes(config)
                    body=bytearray((bad/'manifest.bin').read_bytes())
                    struct.pack_into('=2q',body,100,len(config),crc64(config))
                    (bad/'manifest.bin').write_bytes(body);repair_manifest(bad)
            else:corrupt(original,bad,kind,3)
            flavor='config-failure' if kind in ['missing-commit','manifest','configuration','metadata-short'] else 'state-failure'
            if kind=='settings':flavor='settings-failure'
            code,text=invoke(work,(2,2),'probe','source',flavor,success=kind!='bad-gamma')
            if kind=='bad-gamma':assert 'adiabatic_index debe ser mayor que 1' in text
            else:assert 'RESTART REJECTION PROBE PASSED' in text
            rejection_results.append(dict(kind=kind,returncode=code,distributed_rejected=True))
        assert hashes([original])==original_hash
        for version in [7,9]:
            for kind in ['truncate','version','trailing','mutated-after-configuration']:
                work=out/f'reject_serial_{version}_{kind}';work.mkdir()
                data=bytearray((out/f'case_1_1/legacy_v{version}.rst').read_bytes())
                if kind=='truncate':data=data[:-5]
                elif kind=='version':data[4:16]=b'GRHDCP_V6   '
                elif kind=='trailing':data+=struct.pack('=2i',0,0)
                (work/'source.rst').write_bytes(data)
                flavor='state-failure' if kind=='mutated-after-configuration' else 'config-failure'
                code,text=invoke(work,(2,2),'serial-probe','source.rst',flavor)
                assert 'RESTART REJECTION PROBE PASSED' in text
                rejection_results.append(dict(version=version,kind=kind,returncode=code,serial_rejected=True))
    # Every successful launch compares hashes of all input files before and after.
    source_hashes={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in out.glob('case_*/legacy_v*.rst')}
    (out/'VERIFIED.json').write_text(json.dumps(dict(passed=True,results=results,threshold=1e-11,
        serial_inputs=source_hashes,rejections=rejection_results,
        scope='Physical fields/epoch/injection/counter totals; new diagnostic baselines; no C03 histories or production MPI'),indent=2)+'\n')
    print('PARALLEL RESTART MATRIX PASSED:',out)

if __name__=='__main__':main()
