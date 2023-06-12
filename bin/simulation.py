#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
	Author: Angel Garcia-Casillas del Riego (angel.garcia-casillasdelriego@uzh.ch)
"""

def readConfigFile(filename):
    import configparser
    import sys
    import os

    global env_params

    if os.path.isfile(filename):
        config = configparser.ConfigParser()
        config.read(filename)
        env_params = {'D': 1.0, 'C0': 100.0, 'Kprod': 1.0}
        env_options = config.options("ENVIRONMENT")
        for o in env_options:
            found = 0
            if o == "diff":
                env_params['D'] = config.getfloat("ENVIRONMENT","diff")
                found = 1
            elif o == "init":
                env_params['C0'] = config.getfloat("ENVIRONMENT","init")
                found = 1
            elif o == "kprod":
                env_params['Kprod'] = config.getfloat("ENVIRONMENT","kprod")
                found = 1
            if found == 0:
                print("\033[1mOption %s in section [ENVIRONMENT] of %s file is unknown\033[0m" % (o, filename))
                sys.exit()


class Simulation(object):
    """
    A population-based simulation
    """
    def __init__(self, gen_configs, env_config, **kwargs):
        readConfigFile(env_config)
        self.init_pop(gen_configs)
        self.init_env()

    def init_pop(self, gen_configs):
        from . import Population
        self.pops = [Population(genfile, 100) for genfile in gen_configs]

    def init_env(self):
        from . import Environment
        self.env = Environment(100, env_params)

    def run(self, final_time, time_inc, display=False):
        nits = final_time // time_inc
        popsizes = np.empty(nits)
        for t in range(nits):
            Ct = self.env.update([pop.size for pop in self.pops])
            self.pop_update(Ct)
            popsizes[it] = [pop.size for pop in self.pops]
        if display:
            import matplotlib.pyplot as plt
            plt.plot(np.arange(nits), popsizes)
        return popsizes