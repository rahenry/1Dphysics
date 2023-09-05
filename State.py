#!/usr/bin/env python3
import tenpy.tools.hdf5_io as hdf5
import os
from tenpy.networks.mps import MPS

class State:
    def __init__(self, index, system):
        self.index = index
        self.system = system
        self.psi = None
        # self.log = self.system.log
        self.init_method = system.get('init_method', 'product')
        self.alt_state = system.get('alt_state', '0')
        self.file = os.path.join(self.system.storage_dir, f"{self.system.hash}_state_{self.index}.h5")
        if not system.Q.rerun_simulations:
            self.load()
        if not self.psi:
            self.init()

    def unload(self):
        del self.psi
        self.psi = None

    def load(self):
        # self.log.info(f'Attempting to load {self.index} from {self.file}')
        if self.psi:
            # self.log.info(f'Already have a state')
            return
        if os.path.exists(self.file):
            # self.log.info(f'Loaded state {self.index} from {self.file}')
            try:
                self.psi = hdf5.load(self.file)
                return
            except OSError:
                os.remove(self.file)
                # self.log.info('Bad file')
            except hdf5.Hdf5ImportError:
                os.remove(self.file)
                # self.log.info('Bad file')
        # self.log.info('Not found')

    def save(self):
        if self.psi:
            hdf5.save(self.psi, self.file)

    def init(self):
        # self.log.info(f'Initialising state {self.index}')
        M = self.system.model
        bc = self.system["bc_MPS"]
        if self.init_method == 'product':
            prod = ["0"] * self.system['L']
            prod[0] = self.alt_state
            self.psi = MPS.from_product_state(M.lat.mps_sites(), prod, bc=bc)
        elif self.init_method == 'cycle':
            prod = ["0"] * self.system['L']
            for i in range(self.system['L']):
                N = self.system['N']
                prod[i] = str(i % N)
            self.psi = MPS.from_product_state(M.lat.mps_sites(), prod, bc=bc)
        else:
            print("Unknown init_method")
            exit()
