#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
	Author: Angel Garcia-Casillas del Riego (angel.garcia-casillasdelriego@uzh.ch)
"""
import np
import Environment
import Population

delta_t = 0.01
nparticles = 10
C0 = 1.0
Kprod = 0.05
Gmax = 1.0
Kgrowth = 1.0
Kdet = 1.0
Kdiff = 1.0

def readConfigFile(filename):
	import configparser
    import sys

    global delta_t
    global nparticles
    global C0
    global Kprod
    global Gmax
    global Kgrowth
    global Kdet
    global Kdiff

    if os.path.isfile(filename):
        config = configparser.ConfigParser()
        config.read(filename)
        sim_options = config.options("SIMULATION")
        for o in sim_options:
            found = 0
            if o == "delta_t":
                delta_t = config.getfloat("SIMULATION","delta_t")
                found = 1
            if found == 0:
                print("\033[1mOption %s in section [SIMULATION] of %s file is unknown\033[0m" % (o, filename))
                sys.exit()

        model_args = {}

        pop_options = config.options("POPULATION")
        for o in pop_options:
            found = 0
            if o == "Kgrowth":
                kgrowth = config.getfloat('POPULATION','Kgrowth')
                found = 1
            elif o == "Gmax":
                Gmax = config.getfloat('POPULATION','Gmax')
                found = 1
            elif o == "Kdet":
                Kdet = config.getfloat('POPULATION','Kdet')
                found = 1
            elif o == "Kdiff":
                Kdiff = config.getfloat('POPULATION','Kdiff')
                found = 1
            if found == 0:
                print("\033[1mOption %s in section [POPULATION] of %s file is unknown\033[0m" % (o, filename))
                sys.exit()

        env_options = config.options("ENVIRONMENT")
        for o in env_options:
            found = 0
            if o == "nparticles":
                nparticles = config.getint("ENVIRONMENT","nparticles")
                found = 1
            elif o == "Kprod":
                kprod = config.getfloat("ENVIRONMENT","Kprod")
                found = 1
            elif o == "C0":
                C0 = config.getfloat("ENVIRONMENT","C0")
                found = 1
            if found == 0:
                print("\033[1mOption %s in section [ENVIRONMENT] of %s file is unknown\033[0m" % (o, filename))
                sys.exit()
    else:
        raise IOError(filename + ' was not found or is not a valid file')

class Ecosystem(object):
	"""
	An interface between microbial populations and their abiotic environment
	"""
	def __init__(self, configfile):
		readConfigFile(configfile)
		self.init_agents()
	
	def init_agents(self):
		self.pop = Population(nparticles)
		self.env = Environment()

	def run_dynamical_sim(self):
		change = np.inf
		while not np.isclose(change,0):
			concs = self.env.update()
			change = self.pop.update(concs)
		popstruct = {particle.norgs:particle.conc
			for particle in self.env.particles}
		popstruct['free'] = self.pop.free
		return popstruct