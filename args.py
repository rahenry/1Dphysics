#!/usr/bin/env python3


import argparse
import os

parser = argparse.ArgumentParser()
A = parser.add_argument
A("target", nargs="?", default=None)
parser.add_argument("--rerun_simulations", "-Rs", default=False, action="store_true")
parser.add_argument("--rerun_measurements", "-Rm", default=False, action="store_true")
parser.add_argument("--rerun_exact_values", "-Rev", default=False, action="store_true")
parser.add_argument("--rerun_fidelity", "-Rf", default=False, action="store_true")
parser.add_argument("--rerun_eps", "-Reps", default=False, action="store_true")
A("--storage_root", default="~/.var/tenpy_data", type=str)
A("--qsub", "-q", default=False, action="store_true")
A("--send_qsub", "-s", default=False, action="store_true")
A("--nprocs", "-n", default=1, type=int)
A("--skip_simulations", "-Ss", default=False, action="store_true")
A("--skip_measurements", "-Sm", default=False, action="store_true")
A("--walltime", "-W", type=str, default="5")
A("--retrieve", "-r", default=False, action="store_true")
A("--skip_pgf", "-sp", default=False, action="store_true")
A("--load_only", "-lo", default=False, action="store_true")
A("--profile_run", "-pr", default=False, action="store_true")
A("--log_print", "-lp", default=0, type=int)
A("--storage_prefix", default="~/.var/t4_data")
Q = parser.parse_args()

D = Q.__dict__
for key, val in D.items():
    if isinstance(val, str) and "~" in val:
        D[key] = os.path.expanduser(val)

if Q.nprocs == 1 and D["log_print"] == 0:
    D["log_print"] = 1
