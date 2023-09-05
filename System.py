#!/usr/bin/env python3

from rah_utility import simple_hash, mkdir, silent_remove, rel_diff
from args import Q as cli_args
from State import State
import os
import json
from tenpy.tools.params import Config
from ZnModel import ZnModel
from sxxModel import sxxModel
from solve import solve
from exact_values import make_exact_values, gamma_new, fix_parameters, make_exact_values_post_model
import numpy
import measure
import logging
import dill
import warnings
from tenpy.networks.mps import MPS
from tenpy.linalg import np_conserved as npc
from log import nuHandler
warnings.filterwarnings("ignore", message="unused options for config")

class System:
    def __str__(self):
        return f'System: {self.identifier}'
    def __reprt__(self):
        return f'System: {self.identifier}'

    def generate_hash(self):
        ignore = ['n_states']
        H = dict(self.config.options)
        for x in ignore:
            if x in H:
                del H[x]
        self.hash = simple_hash(json.dumps(H))

    def init_model(self, force=0):
        if not force:
            if self.parent and self.Q.load_only: return
            if self['skip_model']: return
        # print(f"Initisalising model for {self}")
        bc = str(self['bc'])
        if bc == 'infinite':
            self.config['bc_MPS'] = 'infinite'
        else:
            self.config['bc_MPS'] = 'finite'
        if bc == '1':
            self.config['bc_x'] = 'periodic'

        if self['model_type'] == 'sxx':
            self.model = sxxModel(self.config)
        else:
            self.model = ZnModel(self.config)
        if hasattr(self.model, 'export'):
            for x in self.model.export:
                if hasattr(self.model, x):
                    self[x] = getattr(self.model, x)

    def init_states(self):
        if not self['n_states']: return
        if self['skip_model']: return
        # self.log.info(f'Initialising states for {self}')
        self.states = []
        self.states_left = []
        for i in range(self['n_states']):
            s = State(i, self)
            self.states.append(s)

    def get(self, key, default=None):
        return self.config.get(key, default)

    def get_base(self, key):
        if key in self.config_base:
            return self.config_base[key]

    def evaluate_params(self):
        for x, y in self.config_base.items():
            if isinstance(y, str) and '{' in y:
                y = y.strip('{}')
                z = eval(y)
                self.config_base[x] = z

    def __init__(self, data, parent=None, auxillary=0, just_go=0):
        self.config_base = dict(data)
        self.evaluate_params()
        self.identifier = "Unnamed"
        # self.config_base = Config(dict(data), "system_config")
        self.config = Config(dict(data), "system_config")
        self.parent = parent
        self.generate_hash()
        self.fidelity_system = None
        self.states_loaded = False
        self.Q = cli_args
        self.measurement_data = {}
        self.complex_data = {}
        self.other_data = {}
        self.will_run = 0
        self.will_measure = 0
        self.run_all_finished = 0
        self.next = None
        self.previous = None
        self.submitted = 0
        self.prog_tested = 0
        self.chargeless_states = None
        self.chargeless_states_left = None
        self.auxillary = auxillary
        self.model = None
        # self.will_measure = 0

        # output and other file names
        self.sset_name = self.get("sset_name", "unknown")
        self.storage_dir = os.path.join(self.Q.storage_prefix, f'storage/{self.sset_name}')
        self.result_dir = os.path.join(self.Q.storage_prefix, f'results/{self.sset_name}')
        self.log_dir = os.path.join(self.Q.storage_prefix, f'logs/{self.sset_name}')
        mkdir(self.storage_dir)
        mkdir(self.result_dir)
        mkdir(self.log_dir)
        self.measurement_file = os.path.join(self.result_dir, f"{self.hash}_measurements.dill")
        self.other_file = os.path.join(self.result_dir, f"{self.hash}.other")
        if just_go:
            self.prepare()
            self.run()
            self.measure()

    def retrieve_measurements(self):
        if self.auxillary: return
        if not self.parent: return
        if self.parent and not self.parent['retrieve_from']:
            return
        print(f"Attempting retrieval for {self}")
        if not self.model:
            self.init_model()
        d = os.path.join('results', self.parent['retrieve_from'])
        dirpath, dirnames, filenames = next(os.walk(d))

        for f in filenames:
            f = os.path.join(dirpath, f)
            if not 'measure' in f or not '.dill' in f:
                continue
            with open(f, 'rb') as f1:
                data = dill.load(f1)
            self.parent.retrieval_data[f] = data

        keys = ['M', 'cX', 'cZ']
        keys2 = ['e_0']
        found = 0
        for f, data in self.parent.retrieval_data.items():
            if f in self.parent.used_retrievals:
                continue

            success = 1
            tol = 1E-3
            for k in keys2:
                if k not in data:
                    success = 0
            if not success:
                continue
            for k in keys:
                delta = rel_diff(self[k], data[k])
                # if data['M'] == self['M']:
                #     print(k, self[k], data[k], delta)
                if delta < tol:
                    pass
                else:
                    success = 0
            if success:
                found = 1
                self.measurement_data = data
                self.save_measurements()
                self.parent.used_retrievals.add(f)
            else:
                pass
                # print(0)

        if found:
            pass
            # print(self, found)
        else:
            print(self, self.auxillary)
            exit()

    def prepare(self):
        # print(f'Preparing {self}')
        fix_parameters(self)
        self.retrieve_measurements()
        self.load_measurements()
        # print(self)
        # print(self["min_E_distance"])
        # exit()
        self.load_other()
        self.run_test()

        if self.will_run or self.Q.rerun_measurements or self.Q.rerun_exact_values or self['method'] == 'exact':
            if not self.auxillary:
                make_exact_values(self)
            self.init_model()
            if not self.auxillary:
                make_exact_values_post_model(self)
            # self.save_measurements()
        self.log_file = os.path.join(self.log_dir, f"{self.identifier}.log")
        with open(self.log_file, 'w+') as f:
            f.write('')
        self.log = logging.getLogger(self.identifier)
        self.handler = nuHandler()
        # self.log.addHandler(self.handler)
        # self.log.info("Logging started")

        self.auxillary_systems = []
        fp = self['fidelity_parameter']
        if not self.auxillary and self.parent and (fp or self.parent["fidelity_parameters"]) and not self.Q.load_only:
            if not self['fidelity_delta']:
                print('Error: need to specify fidelity_delta')
                exit()
            fps = self.parent["fidelity_parameters"]
            config_fidelity = dict(self.config_base)
            if fps:
                for p in fps:
                    config_fidelity[p] = self[p] + self['fidelity_delta']
            else:
                config_fidelity[fp] = self[fp] + self['fidelity_delta']
            self.fidelity_system = System(config_fidelity, parent=self.parent, auxillary=1)
            self.fidelity_system.prepare()
            self.auxillary_systems.append(self.fidelity_system)
        raf = 1
        if self.will_run or (self.will_measure and not self.auxillary):
            raf = 0
        for s in self.auxillary_systems:
            if s.will_run or not s.run_all_finished:
                raf = 0
        self.run_all_finished = raf


    def rename(self, x):
        self.identifier = x
        if self.fidelity_system:
            self.fidelity_system.identifier = f'{self.identifier}_fidelity'

    def save_other(self):
        with open(self.other_file, 'wb+') as f:
            dill.dump(self.other_data, f)
        # for s in self.auxillary_systems:
        #     s.save_other()

    def save_measurements(self):
        if not self.parent: return
        if self.Q.load_only:
            return
        # for x, y in self.measurement_data.items():
        #     if numpy.iscomplex(y):
        #         self.log.info(f'Tried to save a complex measurement! {x} {y}')
        #         self.measurement_data[x] = y.real
        # with open(self.measurement_file, "wb+") as f:
            # dill.dump(self.measurement_data, f)
        m = {}
        for key, val in self.measurement_data.items():
            try:
                json.dumps(val)
                m[key] = val
            except:
                pass
        with open(self.measurement_file, "w+") as f:
            try:
                json.dump(m, f)
            except TypeError:
                self.log.info(f'Big save fail!!!')
                print(f'Big save fail!!!')

                # for x, y in self.measurement_data.items():
        # for s in self.auxillary_systems:
        #     s.save_measurements()

    def load_measurements(self, force=0):
        if self.Q.rerun_measurements: return

        if not os.path.exists(self.measurement_file):
            self.save_measurements()
        else:
            # try:
            #     with open(self.measurement_file, 'rb') as f1:
            #         self.measurement_data = dill.load(f1)
            # except:
            #     print(f'Loading {self.measurement_file} failed')
            #     self.save_measurements()
            with open(self.measurement_file) as f:
                try:
                    self.measurement_data = json.load(f)
                except json.decoder.JSONDecodeError:
                    # self.log.info(f"Json decode failed {self}; deleting measurement file")
                    os.remove(self.measurement_file)

        for x, y in self.measurement_data.items():
            if y is None:
                self.measurement_data[x] = []

    def load_other(self):
        if not os.path.exists(self.other_file):
            self.save_other()
        else:
            try:
                with open(self.other_file, 'rb') as f1:
                    self.other_data = dill.load(f1)
            except:
                pass
            # with open(self.other_file) as f:
            #     self.other_data = json.load(f)

    def load_states(self):
        if self['skip_model']: return
        if self.parent and self.Q.load_only: return
        if not self.model:
            self.init_model()
        if not hasattr(self, "states"):
            self.init_states()
        else:
            # self.log.info(f'Loading states for {self}')
            for s in self.states:
                s.load()
        for s in self.auxillary_systems:
            s.load_states()

    def save_states(self):
        if not self.parent: return
        if self['skip_model']: return
        if self.parent and self.Q.load_only: return
        for s in self.states:
            s.save()
        # for s in self.auxillary_systems:
        #     s.save_states()

    def unload_states(self):
        for s in self.states:
            s.unload()
        for s in self.auxillary_systems:
            s.unload_states()

    def run_test(self):
        self.will_run = 1
        if not self.parent: return

        if not self.Q.rerun_simulations and self['n_states']:
            finished = True
            for i in range(self['n_states']):
                if not self[f'converged_{i}']:
                    finished = False
            if finished:
                self.will_run = 0
                # self.log.info(f"{self} already converged for all states")
        # print(self.hash, self['e_0'], self['converged_0'], self.will_run)
        if self.Q.rerun_simulations:
            self.will_run = 1
        if self.Q.skip_simulations:
            self.will_run = 0
        # self.will_measure = 1
        # if self['measured'] and not self.Q.rerun_measurements and not self.Q.rerun_simulations:
            self.log.info(f"Already measured {self}")
        #     self.will_measure = 0
        # if self.Q.skip_measurements:
        #     self.will_measure = False
        # if self.will_run:
        #     self.log.info(f"{self} will run")
        # if self.will_measure:
        #     self.log.info(f"{self} will measure")
        if self['e_0'] == [] and not self.auxillary:
            # print(f'{self}, {self["e_0"]}')
            # print(self.measurement_data)
            self.will_measure = 1
        if self.Q.rerun_measurements == 1:
            self.will_measure = 1
        if self.Q.skip_measurements == 1:
            self.will_measure = 0
        if self.Q.load_only:
            self.will_run = 0
            self.will_measure = 0

    def prep_force_run(self):
        self.Q.rerun_simulations = 1
        self.will_run = 1
        self.init_model()
    def force_run(self):
        self.prep_force_run()
        for s in self.auxillary_systems:
            s.prep_force_run()
        self.run()

    def copy(self, changes={}):
        c = self.config_base
        c.update(changes)
        return System(c)

    def run(self):
        if self['method'] == 'none':
            return
        if self.will_run:
            self.load_states()
            self.measurement_data['solve_energies'] = [0]*self['n_states']
            if self['n_states']:
                for i in range(self['n_states']):
                    if self[f'converged_{i}'] and not self.Q.rerun_simulations:
                        continue
                    if self.Q.rerun_simulations and not self['skip_model'] and not self['method'] == 'exact':
                        self.states[i] = State(i, self)
                    if i == 0:
                        if self['detect_phase']:
                            detect_phase2(self)
                    solve(self, i)
                self.save_states()
            else:
                solve(self, 0)
            z = numpy.argsort(self['solve_energies'])
            self.other_data['eig_indices'] = z
            self.save_other()

            # if self['fidelity_parameter'] and not self['method'] == 'exact':
            #     sf = self.auxillary_systems[0]
            #     for i in range(self['n_states']):
            #         psi = self.states[i].psi.copy()
            #         sf.states[i].psi = psi

            for s in self.auxillary_systems:
                s.run()
        self.run_all_finished = 1
        self.dump_log()

    def dump_log(self):
        self.handler.dump(self.log_file)

    def make_left_state(self, psi):
        N = self['N']
        if N == 2:
            return psi
        z = numpy.zeros([N, N])
        z[0][0] = 1
        z[2][1] = 1
        z[1][2] = 1
        o = npc.Array.from_ndarray_trivial(z, labels=['p', 'p*'])
        ops = [o for i in range(self["L"])]
        psi = psi.copy()
        psi.apply_product_op(ops)
        # psi = measure.invert(psi)
        return psi


    def measure(self):
        if not self.will_measure and not self.Q.rerun_measurements:
        # if not self.will_run and not self.Q.rerun_measurements:
            return
        self.load_states()
        measure.measure(self)
        self.save_measurements()
        self.dump_log()

    def __getitem__(self, key, default=[]):
        if key in self.config.options:
            return self.config.get(key, [])
        if key in self.measurement_data:
            return self.measurement_data[key]
        if key in self.complex_data:
            return self.complex_data[key]
        if key in self.other_data:
            return self.other_data[key]
        if hasattr(self, key):
            return getattr(self, key)
        # if self.model and hasattr(self.model, key):
        #     return getattr(self.model, key)
        return default

    def __setitem__(self, key, val):
        try:
            iter(val)
            self.measurement_data[key] = val
            return
        except:
            pass

        data_new = {}
        # try:
        #     iter(val)
        #     data_new[key] = val
        # except:
        #     if numpy.iscomplex(val):
        #         data_new[key] = val
        #         data_new[f'{key}_real'] = val.real
        #         data_new[f'{key}_imag'] = val.imag
        #     else:
        #         if numpy.iscomplexobj(val):
        #             val = val.real
        #         data_new[key] = val
        #         data_new[f'{key}_real'] = val
        #         data_new[f'{key}_imag'] = 0
        if numpy.iscomplex(val):
            data_new[f'{key}_abs'] = abs(val)
            data_new[f'{key}_real'] = val.real
            data_new[f'{key}_imag'] = val.imag
            self.complex_data[key] = val
        else:
            if numpy.iscomplexobj(val):
                val = val.real
            data_new[key] = val
            data_new[f'{key}_real'] = val
            data_new[f'{key}_imag'] = 0
        self.measurement_data.update(data_new)

    def make_label(self, keys):
        res = ''
        for k in keys:
            res += f'{k}={self[k]}, '
        res = res.strip()
        res = res.strip(',')
        return res

    def drop_charges(self):
        if self.chargeless_states:
            return
        self.load_states()
        self.chargeless_states = []
        # c = dict(self.config_base)
        c = dict(self.config)
        c['conserve'] = None
        M = ZnModel(c)
        self.chargeless_model = M
        prod = ["0"] * self['L']
        psi0 = MPS.from_product_state(M.lat.mps_sites(), prod, 'finite')
        for s in self.states:
            psi = s.psi.copy()
            psi._B = [B.drop_charge() for B in psi._B]
            psi.sites = M.lat.mps_sites()
            psi.chinfo = psi0.chinfo
            self.chargeless_states.append(psi)
        for s in self.auxillary_systems:
            s.drop_charges()

        self.chargeless_states_left = []
        for s in self.states_left:
            psi = s.psi.copy()
            psi._B = [B.drop_charge() for B in psi._B]
            psi.sites = M.lat.mps_sites()
            psi.chinfo = psi0.chinfo
            self.chargeless_states_left.append(psi)

if __name__ == "__main__":
    S = System({})
    S.save_states()
    print(S)

def detect_phase(s):
    print('detecting phase')
    1
    c = dict(s.config_base)
    h = c['h']
    l = c['lambda']
    L = 31
    c['L'] = L
    c['bc'] = 0
    c['n_states'] = 2
    c['detect_phase'] = 0
    test_sys = System(c, just_go=1)
    test_sys.run()
    test_sys.measure()
    best = 1E6
    best_state = None
    for i in range(test_sys['n_states']):
        e = test_sys[f'e_{i}']
        # print(i, e)
        if e < best:
            best = e
            best_state = test_sys.states[i].psi
    k = L // 2 - (L//2)%2
    seg = best_state.extract_segment(k, k+1)
    seg.bc = 'infinite'
    s.states[0].psi = seg


def detect_phase2(s):
    return
    print('detecting phase')
    c = dict(s.config_base)
    h = c['lambda']
    l = c['lambda']
    c['h'] = l
    c['detect_phase'] = 0
    test_sys = System(c, just_go=1)
    test_sys.run()
    test_sys.measure()
    # s.states[0].psi = test_sys.states[0].psi
