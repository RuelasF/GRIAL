#!/usr/bin/env python3
"""C03 restart equivalence, histories and corruption, using real independent launches."""
import argparse
from array import array
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import time
from test_checkpoint_mpi import crc64,repair_manifest


def receipt(directory, session):
    return dict(line.split('=',1) for line in (directory/f'session_{session:06}.dat').read_text().splitlines())


def field_state(directory):
    body=(directory/'manifest.bin').read_bytes();nx,ny,nz=struct.unpack_from('=3i',body,36)
    ranks=struct.unpack_from('=i',body,60)[0];cells=nx*ny*nz
    state=array('d',[0.0])*(10*cells);events=[0]*16
    for rank in range(ranks):
        descriptor=struct.unpack_from('=9q',body,116+72*rank)
        _,lx,ly,lz,ox,oy,oz,size,digest=descriptor
        data=(directory/f'block_{rank:09}.bin').read_bytes()
        assert (len(data),crc64(data))==(size,digest)
        local=array('d');local.frombytes(data[292:-16]);assert len(local)==10*lx*ly*lz
        for event,value in enumerate(struct.unpack_from('=16q',data,124)):events[event]+=value
        for v in range(10):
            for k in range(lz):
                for j in range(ly):
                    source=((v*lz+k)*ly+j)*lx
                    dest=((v*nz+k+oz)*ny+j+oy)*nx+ox
                    state[dest:dest+lx]=local[source:source+lx]
    return state,events,(nx,ny,nz)


def rows(path):
    if not path.exists():return []
    return [list(map(float,line.split())) for line in path.read_text().splitlines() if line and not line.startswith('#')]


def compare(reference,actual,segments):
    a=receipt(reference,1);b=receipt(actual,len(segments));assert a['time']==b['time'] and a['step']==b['step']
    x,events,shape=field_state(reference/a['checkpoint']);y,new_events,_=field_state(actual/b['checkpoint'])
    assert x.tobytes()==y.tobytes(),'physical fields differ after segmented continuation'
    assert events[7:13]==new_events[7:13],'physical C2P counters differ'
    root_a=next((reference/'data').iterdir());root_b=next((actual/'data').iterdir())
    maximum=0.0;common=0
    for path in root_a.glob('mpi_*.dat'):
        ra=rows(path);rb=rows(root_b/path.name)
        assert rb and all(rb[i][0]<rb[i+1][0] for i in range(len(rb)-1)),path.name
        with_step=path.name in ['mpi_global_diagnostics.dat','mpi_conservation_budget.dat','mpi_numerical_diagnostics.dat']
        keyed={tuple(row[:2] if with_step else row[:1]):row for row in rb}
        for row in ra:
            key=tuple(row[:2] if with_step else row[:1])
            if key not in keyed:continue
            other=keyed[key];assert len(row)==len(other)
            if path.name=='mpi_numerical_diagnostics.dat':continue # face work legitimately depends on the partition
            scale=max(1.0,*map(abs,row),*map(abs,other))
            error=max(abs(u-v) for u,v in zip(row,other))/scale
            if path.name=='mpi_ppi_modes.dat':
                # A phase is ill-conditioned when C_m approaches zero. Compare
                # the conserved Fourier coefficients and their normalized chord.
                scale=max(abs(row[1]),abs(other[1]),1e-200)
                error=abs(row[1]-other[1])/scale
                for m in range(1,5):
                    q=2+4*(m-1)
                    chord=2*min(row[q+2],other[q+2])*abs(math.sin(0.5*m*(row[q+3]-other[q+3])))
                    error=max(error,abs(row[q]-other[q])/scale,abs(row[q+1]-other[q+1])/scale,
                              abs(row[q+2]-other[q+2]),chord)
                    for sample in [row,other]:
                        phase=-math.atan2(sample[q+1],sample[q])/m
                        assert 2*abs(math.sin(0.5*m*(sample[q+3]-phase)))<1e-11,'phase does not match its own coefficients'
            assert math.isfinite(error) and error<1e-11,(path.name,key,error)
            maximum=max(maximum,error);common+=1
    budget=rows(root_b/'mpi_conservation_budget.dat');assert len(budget)>=2
    scale=max(1e-200,*map(abs,budget[-1][2:12]))
    balance=max(map(abs,budget[-1][-5:]))/scale
    assert balance<1e-11,('persistent residual',balance)
    assert max(abs(u-v) for u,v in zip(budget[0][2:7],budget[-1][2:7]))/scale<1e-11,'initial baseline changed'
    expected_extra=0;previous=0
    for entry in segments:
        step=int(entry['step']);grid=[int(entry['grid_r']),int(entry['grid_theta'])]
        duplicates=(grid[0]-1)*shape[1]*shape[2]+(grid[1]-1)*shape[0]*shape[2]
        expected_extra+=3*(step-previous)*duplicates;previous=step
    numeric=rows(root_b/'mpi_numerical_diagnostics.dat')[-1]
    assert int(numeric[-1])==expected_extra,(numeric[-1],expected_extra)
    event=rows(root_b/'mpi_perturbation_event.dat')
    assert len(event)<=1
    return dict(field_error=0,series_error=maximum,budget_error=balance,common_samples=common,
                steps=int(b['step']),perturbation_events=len(event),extra_faces=expected_extra)


def cases():
    common='''nx=73
ny=1
nz=1
cfl=0.3
final_time=0.12
save_interval=0.05
checkpoint_interval=10
output_folder=data
output_prefix=c03
numerical_diagnostics=true
diagnostic_stride=3
extraction_stride=4
walltime_checkpoint_seconds=0
extract_gw=false
extract_mdot=false
ppi_diagnostics=false
'''
    fixtures=[('wave','problem=convergence\n'),('sod','problem=sod\n'),
              ('khi','problem=khi\nnx=37\nnz=19\n'),
              ('jet','problem=jet\nnx=37\nnz=19\nfinal_time=0.8\n'),
              ('michel','problem=michel\n')]
    for model in ['fd','fm']:
        for spin in [0,0.9]:
            fixtures.append((f'{model}_3d_a{spin}',f'''problem={model}_3d
spin={spin}
metric=kerr_schild
geometry=spheroidal
{model}_inner_radius=8
{model}_pressure_max_radius=12
nx=37
ny=19
nz=8
r_min=1.2
r_max=40
extract_gw=true
extract_mdot=true
ppi_diagnostics=true
perturbation_time=0.03
perturbation_amplitude=0.01
perturbation_seed=3435
'''))
        fixtures.append((model+'_sag',f'''problem={model}_sagittal
spin=0.9
metric=kerr_schild
geometry=spheroidal
{model}_inner_radius=8
{model}_pressure_max_radius=12
nx=37
ny=19
nz=1
r_min=1.2
r_max=40
extract_gw=true
extract_mdot=true
perturbation_time=0.03
perturbation_amplitude=0.01
perturbation_seed=3435
'''))
    for name,extra in fixtures:
        for config in range(3 if name.startswith(('fd','fm')) else 2):
            method='weno5' if config!=1 else 'mp5';solver='hlle' if config!=1 else 'hllc'
            perturbation=['pressure_noise','density_noise','density_mode'][config] if config<3 else 'none'
            if not name.startswith(('fd','fm')):perturbation='none'
            fields=dict(line.split('=',1) for line in (common+extra).splitlines() if line)
            fields.update(reconstruction=method,riemann_solver=solver,shock_sensor='true',perturbation=perturbation)
            text=''.join(key+'='+value+'\n' for key,value in fields.items())
            yield name+'_'+str(config),text


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--driver',type=Path,required=True);ap.add_argument('--probe',type=Path,required=True)
    ap.add_argument('--output',type=Path);ap.add_argument('--serial',action='store_true');ap.add_argument('--smoke',action='store_true')
    args=ap.parse_args();driver=args.driver.resolve();probe=args.probe.resolve()
    out=args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-c03-'))
    out.mkdir(parents=True,exist_ok=True);assert not (out/'PROGRESS.json').exists(),'Use a new output directory'
    env=dict(os.environ,OMP_NUM_THREADS='2',OMP_DYNAMIC='FALSE',OMP_WAIT_POLICY='PASSIVE');results=[]
    def launch(exe,grid,arguments,work,label,success=True):
        command=[] if args.serial else ['mpirun','--oversubscribe','--bind-to','none','-n',str(math.prod(grid))]
        command += [str(exe),*map(str,arguments)]
        log=work/(label+'.log');assert not log.exists()
        start=time.monotonic()
        with log.open('wb') as stream:
            p=subprocess.run(command,cwd=work,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=240)
        assert (p.returncode==0)==success,(command,p.returncode,log.read_text()[-2500:])
        return dict(command=command,returncode=p.returncode,seconds=time.monotonic()-start,log=str(log))
    profiles=[(1,1)] if args.serial or args.smoke else [(1,1),(2,1),(4,1),(1,2),(2,2),(1,4),(4,2),(2,4)]
    for grid in profiles:
        work=out/('unit_'+'_'.join(map(str,grid)));work.mkdir()
        row=launch(probe,grid,[*grid,'roundtrip'],work,'roundtrip')
        assert (work/'probe_result.txt').read_text().strip()=='roundtrip_exact'
        row.update(kind='unit',grid=grid);results.append(row)
    selected=list(cases())
    if args.smoke:selected=[selected[0],selected[2],selected[8],selected[10],selected[-3]]
    for name,text in selected:
        case=out/name;case.mkdir();par=case/'case.par';par.write_text(text)
        reference=case/'reference';reference.mkdir()
        launch(driver,(1,1),['cold',par,1,1,0,0,0,1],reference,'run')
        case_profiles=profiles
        if 'ny=19' not in text and not (args.serial or args.smoke):case_profiles=[(1,1),(2,1),(4,1),(8,1)]
        for grid in case_profiles:
            work=case/('chain_'+'_'.join(map(str,grid)));work.mkdir();entries=[]
            launch(driver,(1,1),['cold',par,1,1,0,0,2,1],work,'cold')
            first=receipt(work,1);entries.append(first)
            if first['state']=='stopped_for_walltime':
                launch(driver,grid,['resume',first['checkpoint'],*grid,0,0,0,2],work,'resume')
                entries.append(receipt(work,2))
            data=compare(reference,work,entries)
            data.update(case=name,grid=grid,kind='continuation');results.append(data)
        (out/'PROGRESS.json').write_text(json.dumps(results,indent=2)+'\n')
        print(name,'PASS',flush=True)
    if not args.serial:
        for mode in ['partial','active','calendar-disagree','path-disagree','status-disagree']:
            work=out/('invalid_'+mode);work.mkdir()
            results.append(launch(probe,(2,1),[2,1,mode],work,'probe',False))
    # An actual root file-size limit truncates only the new observer sidecar.
    io_work=out/'file_limit';io_work.mkdir();io_grid=(1,1) if args.serial else (2,1)
    launch(probe,io_grid,[*io_grid,'empty-produce','previous'],io_work,'produce_previous')
    (io_work/'probe_result.txt').rename(io_work/'previous_result.txt')
    previous={str(p.relative_to(io_work/'previous')):p.read_bytes() for p in (io_work/'previous').iterdir()}
    row=launch(probe,io_grid,[*io_grid,'io','broken'],io_work,'file_limit')
    assert (io_work/'broken/COMMITTED').exists(),'fault did not reach the C03 writer'
    assert not (io_work/'broken/OBSERVED').exists(),'truncated observer was published'
    assert all((io_work/'previous'/name).read_bytes()==data for name,data in previous.items())
    (io_work/'probe_result.txt').rename(io_work/'io_result.txt')
    launch(probe,io_grid,[*io_grid,'accept','previous'],io_work,'previous_recovered')
    assert (io_work/'probe_result.txt').read_text().strip()=='recovered'
    row.update(kind='real_file_limit',limit=2000,previous_recovered=True);results.append(row)
    # The real launcher must distinguish clean quota, completion, and rejected resume.
    launcher=Path(__file__).resolve().parents[1]/'tools/run_mpi_segments.py'
    par=out/selected[0][0]/'case.par';chain=out/'launcher_chain'
    def control(extra,label,expected):
        command=[sys.executable,str(launcher),'--exe',str(driver),'--par',str(par),'--run-dir',str(chain),
                 '--threads','2','--launcher-arg=--oversubscribe','--launcher-arg=--bind-to=none']
        if args.serial:command+=['--serial']
        command+=list(map(str,extra))
        with (out/(label+'.log')).open('wb') as stream:
            process=subprocess.run(command,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=240)
        assert process.returncode==expected,(label,process.returncode,(out/(label+'.log')).read_text()[-2000:])
    sequence=[(1,1)] if args.serial else [(1,1),(2,1),(4,1)]
    for number,grid in enumerate(sequence):
        control((['--resume'] if number else [])+['--grid',*grid,'--max-segments',1,'--test-segment-steps',2],
                'launcher_quota_'+str(number),2)
    final_grid=(1,1) if args.serial else (2,1)
    control(['--resume','--grid',*final_grid,'--max-segments',1],'launcher_finish',0)
    manifest=json.loads((chain/'mpi_segments.json').read_text());assert manifest['state']=='completed'
    entries=[receipt(chain,i+1) for i in range(len(manifest['segments']))]
    results.append(dict(kind='launcher_chain',**compare(out/selected[0][0]/'reference',chain,entries)))
    control(['--resume','--grid',*final_grid],'launcher_already_completed',0)
    assert len(json.loads((chain/'mpi_segments.json').read_text())['segments'])==len(entries)
    # A block modification is rejected before launching another MPI session.
    checkpoint=chain/manifest['segments'][-1]['checkpoint'];block=next(checkpoint.glob('block_*'))
    original_block=block.read_bytes();changed=bytearray(original_block);changed[300]^=1;block.write_bytes(changed)
    control(['--resume','--grid',*final_grid],'launcher_tamper',1)
    block.write_bytes(original_block)
    assert not (chain/f'session_{len(entries)+1:06}.log').exists()
    # True elapsed-time stop on one allocation, then explicit continuation by five steps.
    clock=out/'wall_clock';clock.mkdir();clock_par=clock/'case.par'
    clock_fields=dict(line.split('=',1) for line in list(cases())[4][1].splitlines() if line)
    clock_fields['final_time']='100';clock_par.write_text(''.join(k+'='+v+'\n' for k,v in clock_fields.items()))
    clock_grid=(1,1) if args.serial else (2,1)
    launch(driver,clock_grid,['cold',clock_par,*clock_grid,1.5,0.3,0,1],clock,'clock_stop')
    stopped=receipt(clock,1);assert stopped['state']=='stopped_for_walltime' and int(stopped['step'])>0
    assert float(stopped['elapsed'])<1.5,'wall-clock margin was exhausted'
    launch(driver,clock_grid,['resume',stopped['checkpoint'],*clock_grid,0,0,5,2],clock,'clock_resume')
    continued=receipt(clock,2);reference=out/'wall_clock_reference';reference.mkdir()
    launch(driver,(1,1),['cold',clock_par,1,1,0,0,int(continued['step']),1],reference,'clock_reference')
    results.append(dict(kind='wall_clock',**compare(reference,clock,[stopped,continued])))
    # Corrupt the first cold generation in independent copies; histories are at that epoch.
    original=out/selected[0][0]/'chain_1_1'
    # The finished chain has histories ahead of its first checkpoint, deliberately rejected.
    mutations=['ahead','marker_missing','sidecar_missing','sidecar_truncated','sidecar_flip','physical_missing',
               'physical_flip','physical_truncated','manifest_flip','history_missing','history_changed','calendar','epoch',
               'nonfinite','history_steps','cursor','hist_sig','packet_shape','observer_link']
    first_checkpoint=receipt(original,1)['checkpoint']
    # A fresh cold session at step two gives coherent histories for each corruption test.
    pristine=out/'pristine';pristine.mkdir();par=out/selected[0][0]/'case.par'
    launch(driver,(1,1),['cold',par,1,1,0,0,2,1],pristine,'cold')
    for mutation in mutations:
        source=original if mutation=='ahead' else pristine
        work=out/('corrupt_'+mutation);shutil.copytree(source,work)
        generation=work/first_checkpoint;body=bytearray((generation/'observer.bin').read_bytes())
        if mutation=='marker_missing':(generation/'OBSERVED').unlink()
        elif mutation=='sidecar_missing':(generation/'observer.bin').unlink()
        elif mutation=='sidecar_truncated':(generation/'observer.bin').write_bytes(body[:-7])
        elif mutation=='sidecar_flip':body[-24]^=1;(generation/'observer.bin').write_bytes(body)
        elif mutation=='physical_missing':next(generation.glob('block_*')).unlink()
        elif mutation=='physical_flip':
            path=next(generation.glob('block_*'));data=bytearray(path.read_bytes());data[300]^=1;path.write_bytes(data)
        elif mutation=='physical_truncated':
            path=next(generation.glob('block_*'));path.write_bytes(path.read_bytes()[:-11])
        elif mutation=='manifest_flip':
            path=generation/'manifest.bin';data=bytearray(path.read_bytes());data[44]^=1;path.write_bytes(data)
        elif mutation=='history_missing':(work/'data/weno5_hlle/mpi_conservation_budget.dat').unlink()
        elif mutation=='history_changed':
            path=work/'data/weno5_hlle/mpi_conservation_budget.dat';path.write_bytes(path.read_bytes()+b'1 999\n')
        elif mutation not in ['ahead']:
            # These corruptions keep the sidecar checksum valid to exercise semantic checks.
            offsets={'calendar':(92,'d',-1.0),'epoch':(60,'i',999),'nonfinite':(84,'d',float('nan')),
                     'history_steps':(2200,'q',-1),'cursor':(116,'c',b'X'),
                     'hist_sig':(2232,'q',-5),'packet_shape':(56,'i',1000000000),'observer_link':(2216,'q',42)}
            offset,fmt,value=offsets[mutation];struct.pack_into('='+fmt,body,offset,value)
            (generation/'observer.bin').write_bytes(body)
            marker=bytearray((generation/'OBSERVED').read_bytes());struct.pack_into('=2q',marker,32,len(body),crc64(body))
            (generation/'OBSERVED').write_bytes(marker)
        before={str(p.relative_to(work)):p.read_bytes() for p in work.rglob('*.dat')}
        grid=(1,1) if args.serial else (2,1)
        row=launch(probe,grid,[*grid,'reject',first_checkpoint],work,'reject')
        assert (work/'probe_result.txt').read_text().strip() in ['rejected_without_mutation','configuration_rejected']
        assert all((work/name).read_bytes()==data for name,data in before.items()),'rejected recovery changed histories'
        row.update(kind='corruption',mutation=mutation);results.append(row)
    summary=dict(passed=True,results=results,continuations=sum(x.get('kind')=='continuation' for x in results),
                 corruption_cases=len(mutations),purely_local=True,threshold=1e-11)
    (out/'VERIFIED.json').write_text(json.dumps(summary,indent=2)+'\n');print('C03 MATRIX PASSED',flush=True)

if __name__=='__main__':main()
