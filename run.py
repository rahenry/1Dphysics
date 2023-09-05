#!/usr/bin/env python3


from args import Q
from SystemSet import SystemSet
from Slicer import SystemSlicer
from exact_fp import exact_parafermions_matrix, make_pf_spectrum
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from figures import graph1, graph_all
from qsub import send_qsub, make_qsub
from rah_utility import proc
from output_figs import *
from polytest import *
import subprocess
import numpy
import warnings
warnings.filterwarnings(action='ignore', category=UserWarning)
warnings.filterwarnings(action='ignore', category=FutureWarning)
warnings.filterwarnings(action='ignore', category=numpy.ComplexWarning)
# warnings.filterwarnings(action='ignore', category=ComplexWarning)
from output_new import process_output
from output_tables import process_tables

def retriever(t):
    print(f"Attempting to retrieve data for {t}")
    remote = 'rh8579@gadi.nci.org.au:~/t3'
    target = os.path.join('results', t)
    target = os.path.join(remote, target)
    silent_remove(target)
    # print(target)
    command = ['rsync', '-uvdr', target, 'results']
    command = ' '.join(command)
    p = proc(command, shell=True)
    print(p)

if __name__ == "__main__":
    if Q.retrieve:
        retriever(Q.target)

        input_file = os.path.join("data", Q.target)
        for l in open(input_file).readlines():
            if l.strip() == '': continue
            sp = l.split()
            if 'add' in sp[0]:
                for a in sp[1:]:
                    print(a)
                    retriever(a)
    if Q.send_qsub:
        send_qsub()
        exit()
    if Q.qsub:
        q = make_qsub()
        p = proc(q)
        print(p)
        exit()
    S = SystemSet(Q.target)
    S.prepare()
    S.run()
    S.purge_graphs()
    S.extras()
    process_output(S)
    process_tables(S)



    modes = S['output']
    # if modes == []:
    #     modes = ['standard_output']
    for m in modes:
        x = f'{m}(S)'
        print(f'Output mode {m}')
        try:
            exec(x)
        except NameError:
            print(x)
            # print(f'Unknown output mode {m}')
            raise
        # if m == 'fp_test':
        #     fp_test1(S)
        #     # fp_test2(S)
        # elif m == 'exact_test':
        #     exact_test1(S)
        # elif m == 'XYW_test':
        #     XYW_test(S)
        # elif m == 'XYWp_test':
        #     XYWP_test(S)
        # elif m == 'standard':
        #     standard_output(S)
        # elif m == 'pf_tester':
        #     pf_tester(S)
        # else:
        #     print(f'Unknown output mode {m}')
