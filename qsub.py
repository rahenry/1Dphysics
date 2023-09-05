#!/usr/bin/env python3

from rah_utility import proc, newliner, mkdir
import sys
from args import Q
import os
c = ['nci_ssh']

def send_qsub():
    c = ['sh', 'git_add']
    p = proc(c)
    print(p)
    c = 'gitc'
    p = proc(c)
    print(p)
    cbase = f'''
    module load intel-mkl intel-compiler python3/3.9.2
    pip install pathos sympy
    cd t3
    git pull
    '''
    import subprocess
    command = []
    for line in cbase.split('\n'):
        line = line.strip()
        if line == '': continue
        line = line + ';'
        line = line.split()
        command += line

    ssh = ['ssh', 'rh8579@gadi.nci.org.au']
    p = subprocess.Popen(ssh, stdin=subprocess.PIPE)
    command += ['python3', 'run.py'] + sys.argv[1:]
    command.remove('-s')
    command.append('-q')
    command = bytes(' '.join(command), 'ascii')
    p.stdin.write(command)
    a, b = p.communicate()

def make_qsub():
    q = newliner()
    q += "#!/bin/bash"
    q += "#PBS -P my43"
    q += "#PBS -q normal"
    q += f"#PBS -l walltime={Q.walltime}:00:00"
    q += "#PBS -l mem=50GB"
    ncpus = Q.nprocs
    q += f"#PBS -l ncpus={ncpus}"
    q += "#PBS -l jobfs=40GB"
    q += "#PBS -l storage=scratch/my43"
    q += "#PBS -l wd"
    qsub_folder = os.path.join('qsub', Q.target)
    mkdir(qsub_folder)
    args = sys.argv[1:]
    if '-q' in args:
        args.remove('-q')
    args = ' '.join(args)
    q += f"#PBS -o {os.path.join(qsub_folder, 'qsub_output.o')}"
    q += f"#PBS -e {os.path.join(qsub_folder, 'qsub_output.e')}"
    q += "module load intel-mkl intel-compiler python3/3.9.2"
    q += f'python3 run.py {args}'
    qfile = os.path.join(qsub_folder, 'qsub.qsub')
    open(qfile, 'w+').write(str(q))
    c = ['qsub', qfile]
    return c

if __name__ == "__main__":
    c = ['ssh', 'rh8579@gadi.nci.org.au']

    command = b'''
    echo "123";
    echo "420";
    echo "----------------";
    cd t3
    ls
    git pull
    '''


    c = ['gc']
    p = proc(c)
    import subprocess
    p = subprocess.Popen(c, stdin=subprocess.PIPE)
    p.stdin.write(command)
    a, b = p.communicate()
