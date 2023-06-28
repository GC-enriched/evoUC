#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
	Author: Angel Garcia-Casillas del Riego (angel.garcia-casillasdelriego@uzh.ch)
"""
import numpy as np
import time

mode = 0
delta_t = 0.05

def readConfigFile(filename):
    import configparser
    import sys
    import os

    global env_params
    global delta_t
    global mode

    if os.path.isfile(filename):
        config = configparser.ConfigParser()
        config.read(filename)
        env_params = {'Kdiff': 1.0, 'C0': 100.0, 'Kprod': 1.0}
        env_options = config.options("ENVIRONMENT")
        for o in env_options:
            found = 0
            if o == "diff":
                env_params['Kdiff'] = config.getfloat("ENVIRONMENT","diff")
                found = 1
            elif o == "init":
                env_params['C0'] = config.getfloat("ENVIRONMENT","init")
                found = 1
            elif o == "kprod":
                env_params['Kprod'] = config.getfloat("ENVIRONMENT","kprod")
                found = 1
            elif o == "ks":
                env_params['Ks'] = config.getint("ENVIRONMENT","ks")
                found = 1
            elif o == "beta":
                env_params['beta'] = config.getfloat("ENVIRONMENT","beta")
                found = 1
            elif o == "env_size":
                env_params['env_size'] = config.getint("ENVIRONMENT","env_size")
                found = 1
            if found == 0:
                print("\033[1mOption %s in section [ENVIRONMENT] of %s file is unknown\033[0m" % (o, filename))
                sys.exit()
        solver_options = config.options("SOLVER")
        for o in solver_options:
            found = 0
            if o == "delta_t":
                delta_t = config.getfloat("SOLVER","delta_t")
                found = 1
            elif o == "mode":
                ## Should in theory be a cell parameter but this is more practical
                mode = config.getint("SOLVER","mode")
                found = 1
            if found == 0:
                print("\033[1mOption %s in section [SOLVER] of %s file is unknown\033[0m" % (o, filename))
                sys.exit()
        env_params['nparticles'] = int(env_params['env_size'] *
            env_params['beta'] / env_params['Ks'])

class Simulation(object):
    """
    A population-based simulation
    """
    def __init__(self, gen_configs, env_config, **kwargs):
        readConfigFile(env_config)
        self.ngens = len(gen_configs)
        self.delta_t = delta_t
        self.init_pop(gen_configs)
        self.init_env()
        self.init_pops_in_env()

    def init_pop(self, gen_configs):
        from . import Population
        self.pops = [Population(gen_configs[i],
            env_params['nparticles'], mode=mode)
            for i in range(self.ngens)]

    def init_env(self):
        from . import Environment
        self.env = Environment(env_params)

    def init_pops_in_env(self):
        for pop in self.pops:
            pop.delta_t = self.delta_t
            pop.beta = self.env.beta
            for subpop in pop.subpops:
                subpop.capacity = self.env.Ks
                if mode == 2:
                    subpop.crit = lambda conc: pop.Kdet * (
                        self.env.C0/conc + self.env.C0/(conc + pop.KMC))

    def run(self, final_time, time_inc, display=False):
        nparticles = self.env.nparticles
        nits = final_time // time_inc
        abiotic_freq = nparticles // self.delta_t
        popsizes = [np.empty(nits)] * self.ngens
        freesizes = [np.empty(nits)] * self.ngens
        partpops = np.empty(nits * env_params['nparticles']
            * self.ngens).reshape(self.ngens,
            env_params['nparticles'], nits)
        #resources = np.empty(nits)
        for t in range(nits):
            particle_cells = []
            for j in range(nparticles):
                particle_cells.append(np.sum(
                    [pop.subpops[j].size
                    for pop in self.pops]))
            ct, Ct = self.env.biotic_update(particle_cells)
            self.pop_update(ct, Ct)
            for i in range(self.ngens):
                popsizes[i][t] = self.pops[i].natt
                freesizes[i][t] = self.pops[i].nfree
                for j in range(env_params['nparticles']):
                    partpops[i,j,t] = self.pops[i].subpops[j].size
            #resources[t] = np.sum([part.conc for part in self.env.particles])
            if not t % abiotic_freq:
                self.abiotic_update()
        if display:
            import matplotlib.pyplot as plt
            for i in range(self.ngens):
                plt.plot(np.arange(nits), popsizes[i], 'b-')
                plt.plot(np.arange(nits), freesizes[i], 'y--')
                for j in range(env_params['nparticles']):
                    plt.plot(np.arange(nits), partpops[i,j,:], '-')
            #plt.plot(np.arange(nits), resources, 'r-')
            plt.show()
        return popsizes

    def pop_update(self, ct, Ct):
        """
        Simulate growth for all populations with current nutrient concentrations
        """
        for pop in self.pops:
            pop.update(ct, Ct)

    def abiotic_update(self):
        self.env.abiotic_update()
        for pop in self.pops:
            pop.abiotic_update()