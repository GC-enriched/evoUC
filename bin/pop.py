#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
	Author: Angel Garcia-Casillas del Riego (angel.garcia-casillasdelriego@uzh.ch)
"""

class Population(object):
	"""
	A microbial population that grows on a limiting nutrient
	"""
	def __init__(self, genfile, nparticles, mode=0):
		gen_params = self.read_genfile(genfile)
		for param in gen_params:
			self.__setattr__(param, gen_params[param])
		self.subpops = nparticles
		self.subpops[0].size += 1
		self.freepop = FreePop()
		self.freepop.crit = lambda conc: self.density
		if not mode:
			for subpop in subpops:
				subpop.crit = lambda conc: self.Kdet
		elif mode == 1:
			for subpop in subpops:
				subpop.crit = lambda conc: self.Kdet * self.size
		else:
			for subpop in subpops:
				subpop.crit = lambda conc: self.Kdet * conc

	@property
	def subpops(self):
		return self._subpops
	
	@subpops.setter
	def subpops(self, nparticles):
		self._subpops = [ParticlePop() for n in range(nparticles)]

	@property
	def natt(self):
		return np.sum([subpop.size for subpop in subpops])

	def read_genfile(self, genfile):
		import os
		import configparser
		genotype = {'Kgrowth': 1.0, 'Gmax': 1.0,
			'Kdet': 1.0, 'KMC': 1.0}
		if os.path.isfile(genfile):
			config = configparser.ConfigParser()
			config.read(genfile)
			gen_options = config.options("GENOTYPE")
			for o in gen_options:
				found = 0
				if o == "kgrowth":
					genotype['Kgrowth'] = config.getfloat('CELL','kgrowth')
					found = 1
				elif o == "gmax":
					genotype['Gmax'] = config.getfloat('CELL','gmax')
					found = 1
				elif o == "kdet":
					genotype['Kdet'] = config.getfloat("CELL","kdet")
					found = 1
				elif o == "kmc":
					genotype['KMC'] = config.getfloat("CELL","kmc")
					found = 1
				if found == 0:
					print("\033[1mOption %s in section [CELL] of \
						%s file is unknown\033[0m" % (o, filename))
					sys.exit()
			return genotype
		else:
			raise Exception('File %s not found' % (genfile))

	def Monod(self, nutconc):
		"""
		Estimate biomass growth rate according to a Monod growth model
		"""
		return self.Gmax * self.Kgrowth/(self.Kgrowth + nutconc)

	def update(self, concs):
		change = np.sum([subpop.update(concs[subpop])
			for subpop in concs])
		attachment = self.freepop.flow()
		detachment = np.sum([self.subpops[i].flow()
			for i in range(1, len(self.subpops))])
		self.free.size += detachment
		for subpop in self.subpops:
			self.free.size += detachment/nparticles
		return change

class Subpopulation(Population):
	"""docstring for Subpopulation"""
	def __init__(self, arg):
		self.size = 0

	def flow(self):
		inc = self.crit * self.delta_t
		self.size -= inc

class ParticlePop(Subpopulation):
	"""docstring for ParticlePop"""
	def __init__(self):
		super().__init__()

	def update(self):
		"""
		Estimate population growth under current environmental conditions
		"""
		oldsize = self.size
		delta_biomass = self.Monod(nutconc)
		self.size += delta_biomass * self.delta_t
		change = np.abs(oldsize - self.size)
		return change		
		
class FreePop(Subpopulation):
	"""
	A population of free-living cells
	"""
	def __init__(self):
		super().__init__()

	def update(self):
		"""
		Free-living populations cannot grow. Dummy method returns 0 growth if its called
		"""
		return 0.0