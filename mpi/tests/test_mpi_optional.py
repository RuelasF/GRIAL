#!/usr/bin/env python3
"""M01 MPI lifecycle/error smoke tests. No distributed fluid evolution."""
import argparse
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import tempfile
import time


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--exe',type=Path,required=True)
    ap.add_argument('--probe',type=Path,required=True)
    ap.add_argument('--serial-exe',type=Path)
    ap.add_argument('--mpiexec',default='mpirun')
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    exe,probe=args.exe.resolve(),args.probe.resolve()
    out=args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='grial-mpi-m01-'))
    out.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,OMP_NUM_THREADS='2',OMP_DYNAMIC='FALSE')
    env.pop('OMPI_MPI_THREAD_LEVEL',None)
    results=[]
    def run(name,n,command,success=True,extra_env=None):
        work=out/name;work.mkdir()
        cmd=[args.mpiexec,'--oversubscribe','-n',str(n),*map(str,command)]
        started=time.monotonic()
        process=subprocess.Popen(cmd,cwd=work,env=dict(env,**(extra_env or {})),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,start_new_session=True)
        try:output,_=process.communicate(timeout=20)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGKILL);output,_=process.communicate()
            (work/'run.log').write_bytes(output)
            raise AssertionError('MPI timeout/hang: '+name)
        (work/'run.log').write_bytes(output)
        assert (process.returncode==0)==success,(name,process.returncode,output.decode(errors='replace'))
        assert not (work/'data').exists(),'Runtime checks must not create simulation output'
        results.append(dict(name=name,processes=n,returncode=process.returncode,seconds=time.monotonic()-started))
        return output.decode(errors='replace')
    for n in [1,2,4]:
        output=run(f'probe_{n}',n,[probe,'ok'])
        found=re.findall(r'PARALLEL_TEST rank=(\d+) size=(\d+) workers=(\d+)',output)
        assert sorted(int(rank) for rank,_,_ in found)==list(range(n)),found
        assert all(int(size)==n and int(workers)==2 for _,size,workers in found),found
        assert output.count('PARALLEL RUNTIME TESTS PASSED')==1
        assert 'floating-point exceptions' not in output
        output=run(f'grial_check_{n}',n,[exe,'--parallel-check'])
        assert f'PARALLEL_CHECK mpi=T size={n}' in output and output.count('PARALLEL_CHECK')==1
    for n in [2,4]:
        for mode in ['collective-fail','abort','mpi-error','thread-fail']:
            output=run(f'{mode}_{n}',n,[probe,mode],success=False)
            assert 'PARALLEL ERROR' in output,(mode,output)
        par=out/f'case_{n}.par';par.write_text('problem=sod\nnx=16\nfinal_time=0\noutput_folder=data\n')
        output=run(f'evolution_guard_{n}',n,[exe,par],success=False)
        assert 'Distributed evolution is not enabled yet' in output
        output=run(f'check_guard_{n}',n,[exe,'--check',par],success=False)
        assert 'Distributed evolution is not enabled yet' in output
    if args.serial_exe:
        serial=args.serial_exe.resolve()
        for n in [2,4]:
            output=run(f'nonmpi_executable_guard_{n}',n,[serial,out/f'case_{n}.par'],success=False)
            assert 'This executable has no MPI' in output
    output=run('invalid_usage',1,[exe],success=False)
    assert 'invalid command line' in output
    # Open MPI 5.0.7 rejects overrides different from the requested level.
    version=subprocess.run([args.mpiexec,'--version'],env=env,capture_output=True,text=True,check=True).stdout
    if 'Open MPI' in version:
        output=run('funneled_denied',2,[probe,'ok'],success=False,extra_env={'OMPI_MPI_THREAD_LEVEL':'0'})
        assert 'MPI_THREAD_FUNNELED not provided' in output or ('MPI_Init_thread' in output and 'MPI_ERR_ARG' in output)
        output=run('same_level_accepted',2,[probe,'ok'],extra_env={'OMPI_MPI_THREAD_LEVEL':'1'})
        assert 'PARALLEL RUNTIME TESTS PASSED' in output
    (out/'VERIFIED.json').write_text(json.dumps(dict(passed=True,mpi_version=version.strip(),checks=results,
                                                   distributed_evolution=False),indent=2)+'\n')
    print('OPTIONAL MPI TESTS PASSED:',out)


if __name__=='__main__':main()
