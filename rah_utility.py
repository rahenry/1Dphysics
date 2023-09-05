# Miscellaneous utility functions

import subprocess
import sys
import itertools
import numpy
import os
import shutil
import time
from hashlib import md5
import timeit
import json
import select
import struct
import re
print_output_default = False
pdflatex_command = ['pdflatex', '-interaction',
                    'nonstopmode', '-halt-on-error', '-file-line-error']


def reduce_string(s):
    s = s.replace('_', '')
    return re.sub(r'\W+', '', s.lower())


def now():
    return timeit.default_timer()


def simple_hash(s):
    h = md5()
    h.update(s.encode('ascii'))
    return h.hexdigest()[:10]


def mkdir(d):
    d = os.path.expanduser(d)
    if not d:
        return
    if not os.path.exists(d):
        os.makedirs(d)


bcs = {
    'open': [
        0, '0'],
    'periodic': [
        1, '1'],
}


def make_bc(bc):
    if bc in bcs:
        return bc
    for x, aliases in bcs.items():
        for a in aliases:
            if bc == a:
                return x


def get_values_present(systems, key, allow_duplicates=True):
    vals = set()
    for s in systems:
        val = getattr(s, key)
        if not allow_duplicates and val in vals:
            print("Error: duplicate entry for " + key + str(val))
            return
        vals.add(val)
    return sorted(vals)


def dict_from_eqs(s):
    s = s.split()
    res = {}
    for x in s:
        if len(x.split('=')) != 2:
            return ''
        y, z = x.split('=')
        try:
            z = int(z)
        except ValueError:
            try:
                z = float(z)
            except ValueError:
                return ''
        res[y] = z
    return res


def make_correlation_nice_name(key):
    res = ''
    key = key.split('_')
    res += key[0]
    for x in key[1:]:
        res += '(' + x + ')'
    return res


def update_command_history(root_dir, argv):
    # s = ' '.join(sys.argv)
    run_file = os.path.join(root_dir, 'command_history')
    current_history = open(run_file).read()
    new_entry = 'python3 ' + ' '.join(sys.argv) + '\n'
    if new_entry in current_history:
        current_history = current_history.replace(new_entry, '')
    current_history += new_entry
    open(run_file, 'w+').write(current_history)


def make_N_L_pairs(args):
    Lmaxs = args.Lmax
    if len(args.Lmax) == 1:
        Lmaxs = Lmaxs * len(args.N)

    res = []
    if args.L_specific:
        for (N, L) in list(itertools.product(args.N, args.Lmax)):
            res.append((N, L))
    else:
        for (N, Lmax) in zip(args.N, Lmaxs):
            for L in range(args.Lmin, Lmax+1):
                res.append((N, L))
    return res


def silent_remove(filename):
    filename = os.path.expanduser(filename)
    try:
        os.remove(filename)
        # print("deleted " + filename)
    except OSError:
        try:
            shutil.rmtree(filename)
        except OSError:
            pass


def approximately_integer(x):
    if rel_diff(x, int(x)) < 1E-14:
        return True
    return False


def rel_diff(x, y, tol=1E-15):
    if isinstance(x, str):
        x = float(x)
    if isinstance(y, str):
        y = float(y)
    # x = abs(x)
    # y = abs(y)
    if abs(x-y) < tol:
        return 0
    if x == 0 and y == 0:
        return 0
    return abs(x-y)/(abs(x)+abs(y))


def c_numeric_output_convert(z):
    if ',' not in z:
        return(float(z))
    z = z.strip('()')
    imag_neg = '-' in z[1:].replace('e-', '>????')
    if imag_neg:
        z = z.replace(',', '')
    else:
        z = z.replace(',', '+')
    z = z + 'j'
    return complex(z)


def get_value_range(runs, key):
    examples = set([])
    for r in runs:
        x = r.get(key)
        if x is not None:
            examples.add(x)
    return sorted(examples)


def update_general(task):
    out = task.proc.stdout
    stdout_new = task.stdout_pending
    time.sleep(0.2)
    while (True):
        rlist, wlist, elist = select.select([out], [], [], 0)
        if rlist == []:
            break
        line = task.proc.stdout.readline().decode()
        stdout_new += line
        if out.print_output:
            print('"'+line.strip('\n'))
        if line == '':
            break

    if "JSON_START" in stdout_new:
        if "JSON_END" not in stdout_new:
            task.stdout_pending = stdout_new
        jsons = []
        for x in stdout_new.split("JSON_START"):
            if "JSON_END" not in x:
                continue
            x = x.split("JSON_END")[0]
            jsons.append(x.strip('\n'))
        for j in jsons:
            if j.strip('\n') == '':
                print("Bad json?")
                continue
            j = json.loads(j)
            task.latest_json = j
            if "write_to" not in j:
                continue
            for key, val in j.items():
                if key == "write_to":
                    continue
                setattr(task, key, val)
                getattr(task.ob, j["write_to"])[key] = val
            if "just_saved" in j and j["just_saved"]:
                # print("saving state")
                task.ob.save()
    task.stdout += stdout_new


def make_complex_str(z):
    return str(z).replace('j', '').replace('+', ',').replace('-', ',-')


def apply_jsons(stdout_new, task=None, ob=None):
    if "JSON_START" in stdout_new:
        jsons = []
        for x in stdout_new.split("JSON_START"):
            if not "JSON_END" in x:
                continue
            x = x.split("JSON_END")[0]
            jsons.append(x.strip('\n'))
        for j in jsons:
            if j.strip('\n') == '':
                print("Bad json?")
                continue
            j = json.loads(j)
            if task:
                task.latest_json = j
            if not "write_to" in j:
                continue
            for key, val in j.items():
                if key == "write_to":
                    continue
                if task:
                    setattr(task, key, val)
                if ob:
                    getattr(ob, j["write_to"])[key] = val
            if "just_saved" in j and j["just_saved"]:
                # print("saving state")
                if ob:
                    ob.save()
    if task:
        task.stdout += stdout_new


def normalised(s):
    norm = numpy.linalg.norm(s)
    if norm < 1E-12:
        return s
    return s / norm


inventory_types = ['MISC', 'WEAP', 'ARMO', 'CLOT', 'REPA',
                   'APPA', 'LOCK', 'PROB', 'INGR', 'BOOK', 'ALCH', 'LEVI']
naming_types = ['FNAM', 'NAME']
NULL = '\x00'
NEWLINE = '\r\n'

attributes = ['strength', 'speed', 'agility', 'willpower',
              'endurance', 'intelligence', 'luck', 'personality']


def mkdir(d):
    d = os.path.expanduser(d)
    if not d:
        return
    if not os.path.exists(d):
        os.makedirs(d)


required_directories = [
    'report',
]
for r in required_directories:
    mkdir(r)


def fix_file_name(file_name):
    if not file_name:
        return None
    return os.path.abspath(os.path.expanduser(file_name))


def read_newline_separated(file_name):
    data = open(file_name).readlines()
    data = [x.strip() for x in data]
    pieces = []
    piece = []
    for d in data:
        if len(d) > 0 and d[0] == '#':
            continue
        if d == '' and piece != []:
            pieces.append(piece)
            piece = []
            continue
        if d != '':
            piece.append(d)
    if piece not in pieces:
        pieces.append(piece)

    return pieces


repeating_data = {
    'ARMO': {'INDX': ['BNAM', 'CNAM']},
    'CLOT': {'INDX': ['BNAM', 'CNAM']},
    'FACT': {'ANAM': ['INTV']},
    'NPC_': {'NPCO': [], 'NPCS': [], 'DODT': ['DNAM'], 'AI_T': ['AI_W', 'AI_F'], 'AI_W': ['AI_T', 'AI_F']},
    'CREA': {'NPCO': []},
    'CONT': {'NPCO': []},
}


def get_repeating_data(t):
    if t not in repeating_data:
        repeating_data[t] = {}
    return repeating_data[t]


def proc(c, cwd='.', print_output=False, env=os.environ.copy(), **kwargs):
    mkdir(cwd)
    if "shell" not in kwargs and isinstance(c, str):
        c = c.split()
    p = subprocess.Popen(c, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, cwd=cwd, env=env, **kwargs)
    res = ''
    while True:
        line = p.stdout.readline()
        line = line.decode()
        if print_output:
            sys.stdout.write(line)
        if line == '' and p.poll() != None:
            break
        if line:
            res += line
    return_code = p.wait()
    return res


def reduce_name(x):
    remove = ['_', ' ', "'", '"', '-']
    x = x.lower()
    for r in remove:
        x = x.replace(r, '')
    return x

def redeq(x, y):
    return reduce_name(x) == reduce_name(y)


def unpack_binary(data, t, names):
    if t == 'c':
        t = '32c'
    return struct.unpack('<'+t, data)


def make_binary_record(data, t):
    return t.upper().encode('cp1252') + struct.pack('<iii', len(data), 0, 0) + data


def make_binary_subrecord(data, t):
    return t.upper().encode('cp1252') + struct.pack('<i', len(data)) + data


def find_previous(l, x):
    if l == []:
        return None
    if l[0] == x:
        return None
    for i in range(1, len(l)-1):
        if l[i+1] == x:
            return l[i]


def insert_after(l, x, y):
    if x not in l:
        return
    j = None
    for i in range(len(l)):
        if l[i] == x:
            j = i
            break

    return l[:j+1] + [y] + l[j+1:]


def is_int(s):
    try:
        x = int(s)
        return True
    except ValueError:
        return False


def is_numeric(s):
    try:
        x = float(s)
        return True
    except ValueError:
        return False


def convert_to_number_if_possible(x):
    try:
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return x


def flatten(l):
    return [item for sublist in l for item in sublist]


def printvars(x):
    for i, j in vars(x).items():
        j = str(j)
        j = j.split('\n')[0]
        print(f'{i:<20} {j:<10} ')
        # print(f"{i} {j}")
        # print(f'{i:<10} {j:<30}')


def cmnd(c, D=None, shell=False):
    if D:
        mkdir(D)
    if shell:
        c = ' '.join(c)
    P = subprocess.Popen(c, cwd=D, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, shell=shell)
    stdout, stderr = P.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    return stdout, stderr


def listify(x):
    if not isinstance(x, list):
        x = [x]
    return x


def restore_from_backup(target_dir, backup_dir, delete_extras=True):
    base_b, dirs_b, files_b = next(os.walk(backup_dir))
    base_t, dirs_t, files_t = next(os.walk(target_dir))
    for f in files_t:
        f_t = os.path.join(base_t, f)
        f_b = os.path.join(base_b, f)
        if f in files_b:
            s_t = os.stat(f_t)
            s_b = os.stat(f_b)
            if s_t.st_size != s_b.st_size:
                print(f'Reverting file {f_t} to {f_b}')
                c = ['cp', f_b, f_t]
                cmnd(c)
        else:
            if delete_extras:
                print(f'Deleting file {f_t}')
                silent_remove(f_t)
    for f in files_b:
        f_t = os.path.join(base_t, f)
        f_b = os.path.join(base_b, f)
        if not os.path.exists(f_t):
            print(f'Copying file {f_b} to {f_t}')
            c = ['cp', f_b, f_t]
            cmnd(c)

    for d in dirs_b:
        d_t = os.path.join(base_t, d)
        d_b = os.path.join(base_b, d)
        if not os.path.exists(d_t):
            print(f'Restoring dir {d_t}')
            c = ['cp', '-dr', d_b, d_t]
            cmnd(c)
    for d in dirs_t:
        d_t = os.path.join(base_t, d)
        d_b = os.path.join(base_b, d)
        if d in dirs_b:
            # print(f'Descending into {d_b}')
            restore_from_backup(d_t, d_b, delete_extras=delete_extras)
        else:
            if delete_extras:
                print(d, d_t, d_b)
                print(f'Deleting dir {d_t}')
                silent_remove(d_t)


def lowercase_rename(dir):
    def rename_all(root, items):
        for name in items:
            try:
                os.rename(os.path.join(root, name),
                          os.path.join(root, name.lower()))
            except OSError:
                pass  # can't rename it, so what

    for root, dirs, files in os.walk(dir, topdown=False):
        rename_all(root, dirs)
        rename_all(root, files)


def copy_tree(src, dest, overwrite=False):
    mkdir(dest)
    for (dirpath, dirname, filenames) in os.walk(src):
        r = os.path.relpath(dirpath, src)
        y = os.path.join(dest, r)
        mkdir(y)
        for f in filenames:
            x = os.path.join(y, f)
            if os.path.exists(x):
                if overwrite:
                    os.remove(x)
                    shutil.copyfile(os.path.join(dirpath, f), x)
                else:
                    pass
            else:
                shutil.copyfile(os.path.join(dirpath, f), x)


def inheritors(klass):
    subclasses = set()
    work = [klass]
    while work:
        parent = work.pop()
        for child in parent.__subclasses__():
            if child not in subclasses:
                subclasses.add(child)
                work.append(child)
    return subclasses

class newliner:
    def __init__(self):
        self.lines = []

    def __str__(self):
        return '\n'.join(self.lines)+'\n'

    def __add__(self, other):
        self.lines.append(str(other))
        return self


def test_list_type(list1, t):
    for x in list1:
        try:
            t(x)
        except ValueError:
            return False
    return True

def interpret_float_range(data):
    if len(data) == 3:
        n = data[-1]
        if (rel_diff(n, int(n)) < 1E-15) and n >= 3:
            n = int(n)
            return numpy.linspace(data[0], data[1], n)
        else:
            return data
    else:
        return data

def interpret_int_range(data):
    if len(data) == 3:
        if data[2] < data[1] and data[1] > data[0]:
            data = numpy.arange(data[0], data[1], data[2])
            data = [int(x) for x in data]
    return data

def interpret_parameter_line(line, use_ranges=True):
    line = line.strip('\n')
    line_name = line.split()[0]
    data = line.split()[1:]
    test_types = [int, float, str]
    T = None
    for t in test_types:
        if test_list_type(data, t):
            T = t
            break
    if T is None:
        print(f"Couldn't interpret {line}")
    data = [t(x) for x in data]

    if use_ranges:
        if t is float:
            data = interpret_float_range(data)
        if t is int:
            data = interpret_int_range(data)

    return [line_name, list(data)]

def read_parameter_lines_str(f):
    res = {}
    for l in f:
        if l.strip() == '':
            continue
        z = interpret_parameter_line(l)
        res[z[0]] = z[1]

    return res

def read_parameter_lines(f):
    res = {}
    for l in open(f):
        if l.strip() == '':
            continue
        z = interpret_parameter_line(l)
        res[z[0]] = z[1]

    return res

def dict_product(d):
    res = [{}]
    for key, vals in d.items():
        res_new = []
        for r in res:
            for v in vals:
                x_new = dict(r)
                x_new[key] = v
                res_new.append(x_new)
        res = res_new
    return res

def matplotlib_set_aspect(ax, ratio):
    # ratio = 1.0
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)

def extract_archive(f, into='.'):
    mkdir(into)
    c = None
    if '.zip' in f:
        c = ['unzip', f]
    elif '.rar' in f:
        c = ['unrar', 'x', f]
    else:
        print(f'Uknown archive {f}')

    if c:
        p = proc(c, cwd=into)
        return True

def extrema(x):
    return numpy.array([numpy.min(x), numpy.max(x)])

def deduce_raster_params(X):
    """
    Computes raster dimensions based on min/max coordinates in X
    sample step computed from 2nd - smallest coordinate values
    """
    things = [numpy.unique(v) for v in X.T]
    unique_sorted = numpy.vstack(things).T
    d_min = unique_sorted[0] # x min, y min
    d_max = unique_sorted[-1] # x max, y max
    d_step = unique_sorted[1]-unique_sorted[0] # x, y step
    nsamples = (numpy.round((d_max - d_min) / d_step) + 1).astype(int)
    return d_min, d_max, d_step, nsamples

def rastify(X, y):
    X = numpy.array(X)
    d_min, d_max, d_step, nsamples = deduce_raster_params(X)
    ind = numpy.round((X - d_min) / d_step).T.astype(int)
    ind = ind.T
    z = numpy.zeros(nsamples)
    for i in range(len(X)):
        a, b = ind[i]
        z[b][a] = y[i]


    extent = numpy.vstack((d_min, d_max)).T.ravel()
    z = numpy.array(z)
    return z, extent

def cartesian_product(*arrays):
    la = len(arrays)
    dtype = numpy.result_type(*arrays)
    arr = numpy.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(numpy.ix_(*arrays)):
        arr[...,i] = a
    return arr.reshape(-1, la)


def find_minima(z):
    X, Y = z.shape
    res = []
    for i in range(X-2):
        x = i+1
        for j in range(Y-2):
            y = j+1
            u = z[x][y+1]
            d = z[x][y-1]
            l = z[x-1][y]
            r = z[x+1][y]
            t = z[x][y]
            things = [-1, 0, 1]
            accepted = 1
            for ox in things:
                for oy in things:
                    if oy == ox == 0:
                        continue
                    if t > z[x+ox][y+oy]:
                        accepted = 0

            if accepted:
                res.append((x, y))
            # if t < u and t < d and t < l and t < r:
            #     res.append((x, y))

    return res

def extract_archive(f, target):
    extractors = {
        'rar': ['unrar', 'x'],
        '7z': ['7za', 'x'],
        'zip': ['unzip'],
        'tar': ['tar', '-xf'],
    }
    data = None
    for ex, d in extractors.items():
        if ex in reduce_name(f):
            data = d
            break
    if data is None:
        print(f"Could not extract {f}")

    c = [*data, f]
    print(c)
    p = proc(c, cwd=target)
    print(p)
