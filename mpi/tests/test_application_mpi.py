#!/usr/bin/env python3
"""Exercise the actual grial distributed entry, including restart and rejection."""
import argparse
from array import array
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from test_observer_mpi import cases, compare, field_state
from test_walltime_restart import state as serial_state,records

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from parallel_vtk import collect_frames, read_frame


def status(work):
    paths = list((work/'data').glob('*/run_status.dat'))
    assert len(paths) == 1
    return dict(line.split('=', 1) for line in paths[0].read_text().splitlines())


def adapt_receipt(work, entry, session=1):
    data = dict(entry, time=entry['physical_time'])
    (work/f'session_{session:06}.dat').write_text(''.join(f'{k}={v}\n' for k,v in data.items()))
    return data


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--exe', type=Path, required=True)
    ap.add_argument('--driver', type=Path)
    ap.add_argument('--baseline', type=Path)
    ap.add_argument('--output', type=Path)
    ap.add_argument('--serial', action='store_true')
    ap.add_argument('--smoke', action='store_true')
    args = ap.parse_args()
    exe = args.exe.resolve()
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-control-c-'))
    out.mkdir(parents=True, exist_ok=True)
    assert not (out/'PROGRESS.json').exists(), 'Use a new output directory'
    env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE', OMP_WAIT_POLICY='PASSIVE', OPENBLAS_NUM_THREADS='1')
    results, rejections = [], []

    def launch(work, label, arguments, grid=(1,1), success=True, program=exe, threads=2):
        prefix = [] if args.serial else ['mpirun','--oversubscribe','--bind-to','none','-n',str(math.prod(grid))]
        log = work/(label+'.log')
        assert not log.exists()
        with log.open('wb') as stream:
            p = subprocess.run([*prefix,str(program),*map(str,arguments)],cwd=work,
                               env=dict(env,OMP_NUM_THREADS=str(threads)),
                               stdout=stream,stderr=subprocess.STDOUT,timeout=240)
        assert (p.returncode==0)==success, (label,p.returncode,log.read_text()[-3500:])
        if not success: rejections.append(label)

    selected = {'wave_0','sod_0','khi_0','jet_1','michel_0','fd_3d_a0_0',
                'fd_3d_a0.9_1','fm_3d_a0_2','fm_3d_a0.9_0','fm_sag_1'}
    if args.smoke: selected = {'sod_0','fd_3d_a0.9_1'}
    assert selected <= set(dict(cases())), 'Unknown fixture name'
    for name, text in cases():
        if name not in selected: continue
        reference = out/(name+'_reference');reference.mkdir()
        par = reference/'case.par';par.write_text(text)
        launch(reference,'check',['--distributed','--check',par])
        assert not (reference/'data').exists()
        launch(reference,'cold',['--distributed',par])
        ref = status(reference);assert ref['status']=='completed'
        adapt_receipt(reference,ref)
        native,_,shape = field_state(reference/ref['checkpoint'])
        grids = [(1,1)] if args.serial else ([(2,2),(8,1)] if shape[1]>1 else [(2,1),(4,1),(8,1)])
        if args.smoke: grids = grids[:1]
        for grid in grids:
            cold=out/(name+'_cold_'+str(grid[0])+'x'+str(grid[1]));cold.mkdir()
            cold_par=cold/'case.par';cold_par.write_text(text)
            launch(cold,'cold',['--distributed',cold_par,'--grid',*grid],grid)
            cold_status=status(cold)
            assert cold_status['status']=='completed' and cold_status['physical_time']==ref['physical_time']
            cold_values,_,cold_shape=field_state(cold/cold_status['checkpoint'])
            assert shape==cold_shape and len(native)==len(cold_values)
            cold_error=max(abs(u-v)/max(1.,abs(u),abs(v)) for u,v in zip(native,cold_values))
            assert cold_error < 1e-11, (name,grid,'cold field error',cold_error)
            work = out/(name+'_resume_'+str(grid[0])+'x'+str(grid[1]));work.mkdir()
            par = work/'case.par';par.write_text(text)
            launch(work,'stop_initial',['--distributed',par,'--walltime','1e-9','--margin','5e-10'])
            first = status(work);assert first['status']=='stopped_for_walltime' and int(first['step'])==0
            old = work/first['checkpoint']
            assert not list((work/'data').glob('*/.mpi_run_lock'))
            before = {str(p.relative_to(work)):hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (work/'data').rglob('*') if p.is_file()}
            launch(work,'check_restart',['--distributed','--resume','--check',old,'--grid',*grid],grid)
            assert before == {str(p.relative_to(work)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in (work/'data').rglob('*') if p.is_file()}
            launch(work,'resume',['--distributed','--resume',old,'--grid',*grid],grid)
            end = status(work);assert end['status']=='completed'
            adapt_receipt(work,end,2)
            comparison = compare(reference,work,[adapt_receipt(work,first),end])
            fields = next((work/'data').glob('*/mpi_fields'))
            frames = collect_frames(fields,fields/'frames.pvd')
            actual = read_frame(fields/f'step_{int(end["step"]):010}',reconstruct=True)
            assert tuple(actual['shape']) == shape
            assert not list((work/'data').glob('*/.mpi_run_lock'))
            results.append(dict(case=name,grid=grid,cold_field_error=cold_error,
                                cold_fields_exact=native.tobytes()==cold_values.tobytes(),**comparison))
            if args.smoke: break
        if args.driver:
            work = out/(name+'_driver');work.mkdir();par=work/'case.par';par.write_text(text)
            launch(work,'driver',['cold',par,1,1,0,0,0,1],program=args.driver.resolve())
            # The existing driver exposes all native UP/P values, not a raster oracle.
            driver_receipt = dict(line.split('=',1) for line in (work/'session_000001.dat').read_text().splitlines())
            old,_,_ = field_state(work/driver_receipt['checkpoint'])
            assert native.tobytes()==old.tobytes(), (name,'entry differs from C03/C04 driver')
        if args.baseline:
            a = out/(name+'_legacy_before');b = out/(name+'_legacy_after')
            for work, program in [(a,args.baseline.resolve()),(b,exe)]:
                work.mkdir();par=work/'case.par';par.write_text(text)
                # Strict legacy bytes use one thread; atomics may reorder FP sums with two.
                launch(work,'legacy',[par],program=program,threads=1)
            def legacy(work):
                info=dict(line.split('=',1) for line in (work/'data/run_status.dat').read_text().splitlines())
                return serial_state(work/info['checkpoint'])
            assert legacy(a)==legacy(b), (name,'legacy entry changed')
            def checkpoint(work):
                info=dict(line.split('=',1) for line in (work/'data/run_status.dat').read_text().splitlines())
                return work/info['checkpoint']
            assert records(checkpoint(a))==records(checkpoint(b)), (name,'legacy V9 records changed')

    # Real wall-clock stop with progress, followed by a different process grid.
    if not args.smoke:
        text = dict(cases())['khi_0']
        text = text.replace('nx=37','nx=193').replace('nz=19','nz=97').replace('final_time=0.12','final_time=0.4')
        reference=out/'wall_reference';reference.mkdir();par=reference/'case.par';par.write_text(text)
        launch(reference,'continuous',['--distributed',par]);ref=status(reference);adapt_receipt(reference,ref)
        work=out/'wall_segmented';work.mkdir();par=work/'case.par';par.write_text(text)
        launch(work,'stop',['--distributed',par,'--walltime','1.5','--margin','0.3'])
        first=status(work);assert first['status']=='stopped_for_walltime' and 0<int(first['step'])<int(ref['step'])
        grid=(1,1) if args.serial else (2,1)
        launch(work,'finish',['--distributed','--resume',work/first['checkpoint'],'--grid',*grid],grid)
        end=status(work);adapt_receipt(work,end,2)
        results.append(dict(case='real_walltime',grid=grid,**compare(reference,work,[adapt_receipt(work,first),end])))

    controls=out/'rejections';controls.mkdir();par=controls/'case.par'
    par.write_text(dict(cases())['sod_0'])
    grid=(1,1) if args.serial else (2,1)
    for label, options in [
        ('missing_source',[]),('unknown_option',['--unknown']),('missing_grid',['--grid',2]),
        ('bad_grid',[par,'--grid',-1,2]),('wrong_grid',[par,'--grid',3,3]),
        ('missing_file',['missing.par']),('negative_clock',[par,'--walltime',-1]),
        ('nan_clock',[par,'--walltime','NaN']),('invalid_margin',[par,'--walltime',1,'--margin',2]),
        ('duplicate_option',[par,'--check','--check']),('missing_checkpoint',['--resume','missing']),
        ('unexpected_source',[par,par])]:
        launch(controls,label,['--distributed',*options],grid,False)
    bad=controls/'bad.par';bad.write_text('problem=sod\nnot_a_parameter=3\n')
    launch(controls,'bad_parameter',['--distributed',bad],grid,False)
    bad.write_text(dict(cases())['sod_0'].replace('numerical_diagnostics=true','numerical_diagnostics=false'))
    launch(controls,'missing_d01',['--distributed',bad],grid,False)
    original=next(p for p in out.glob('*_reference') if (p/'data').exists())
    launch(original,'refuse_overwrite',['--distributed',original/'case.par'],success=False)
    source=original/status(original)['checkpoint'];corrupt=original/'corrupt';shutil.copytree(source,corrupt)
    piece=next(corrupt.glob('block_*.bin'));piece.write_bytes(piece.read_bytes()[:100])
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (original/'data').rglob('*') if p.is_file()}
    launch(original,'corrupt_restart',['--distributed','--resume','--check',corrupt,'--grid',*grid],grid,False)
    assert hashes=={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (original/'data').rglob('*') if p.is_file()}
    locked=controls/'locked';locked.mkdir();locked_par=locked/'case.par';locked_par.write_text(dict(cases())['sod_0'])
    (locked/'data/weno5_hlle/.mpi_run_lock').mkdir(parents=True)
    launch(locked,'locked_output',['--distributed',locked_par],grid,False)
    assert not list((locked/'data').glob('*/mpi_fields'))
    blocked=controls/'blocked';blocked.mkdir();blocked_par=blocked/'case.par';blocked_par.write_text(dict(cases())['sod_0'])
    (blocked/'data').write_text('preserve this file\n')
    launch(blocked,'unwritable_output',['--distributed',blocked_par],grid,False)
    assert (blocked/'data').read_text()=='preserve this file\n'
    if not args.serial:
        for label, second in [('different_arguments',['--parallel-check']),('different_paths',['--distributed','other.par'])]:
            command=['mpirun','--oversubscribe','--bind-to','none','-n','1',str(exe),'--distributed',str(par),
                     ':','-n','1',str(exe),*map(str,second)]
            with (controls/(label+'.log')).open('wb') as stream:
                p=subprocess.run(command,cwd=controls,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=30)
            assert p.returncode!=0;rejections.append(label)
        directories=[controls/'rank0',controls/'rank1']
        for k, directory in enumerate(directories):
            directory.mkdir();(directory/'case.par').write_text(dict(cases())['sod_0']+f'# rank {k}\n')
        command=['mpirun','--oversubscribe','--bind-to','none']
        for directory in directories:
            if directory!=directories[0]:command+=[':']
            command+=['-n','1','--wdir',str(directory),str(exe),'--distributed','--check','case.par']
        with (controls/'different_contents.log').open('wb') as stream:
            p=subprocess.run(command,cwd=controls,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=30)
        assert p.returncode!=0;rejections.append('different_contents')
    report=dict(passed=True,comparisons=results,rejections=rejections,legacy_checked=bool(args.baseline),
                driver_checked=bool(args.driver),real_walltime_checked=not args.smoke)
    (out/'PROGRESS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=True,comparisons=len(results),rejections=len(rejections))))


if __name__=='__main__': main()
