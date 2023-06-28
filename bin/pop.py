#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
	Author: Angel Garcia-Casillas del Riego (angel.garcia-casillasdelriego@uzh.ch)
"""
import sys
import numpy as np
import scipy.stats as st

class Population(object):
	"""
	A microbial population that grows on a limiting nutrient
	"""
	def __init__(self, genfile, nparticles, mode=0):
		gen_params = self.read_genfile(genfile)
		for param in gen_params:
			self.__setattr__(param, gen_params[param])
		self.subpops = nparticles
		self.subpops.reverse()
		self.freepops = nparticles
		for j in range(nparticles):
			self.subpops[j].mother_pop = self
			if not mode:
				self.subpops[j].crit = lambda conc: self.Kdet
			elif mode == 1:
				self.subpops[j].crit = lambda subpop: self.Kdet * subpop.size
			else:
				self.subpops[j].crit = lambda conc: self.Kdet * conc
			self.freepops[j].crit = lambda conc: self.density
		## Assume the number of cells per particle after dispersal
		## is Zipf-distributed, with multicellularity parameter
		## ZMC determining the slope of the distribution
		self.cell_dist = st.zipf(self.ZMC)

	@property
	def subpops(self):
		return self._subpops
	
	@subpops.setter
	def subpops(self, nparticles):
		self._subpops = [ParticlePop(1) if n < self.ncells
			else ParticlePop(0) for n in range(nparticles)]

	@property
	def freepops(self):
		return self._freepops
	
	@freepops.setter
	def freepops(self, nparticles):
		"""
		Preassign free-living cells to the vicinity of
		the particle that they are going to occupy
		This helps assign free-living collectives to a single particle
		"""
		self._freepops = [FreePop()
			for n in range(nparticles)]

	@property
	def npops(self):
		return len(self.subpops)

	@property
	def natt(self):
		return np.sum([subpop.size for subpop in self.subpops])

	@property
	def nfree(self):
		return np.sum([freepop.size for freepop in self.freepops])

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
				if o == "ncells":
					genotype['ncells'] = config.getint("GENOTYPE",'ncells')
					found = 1
				elif o == "kgrowth":
					genotype['Kgrowth'] = config.getfloat("GENOTYPE",'kgrowth')
					found = 1
				elif o == "gmax":
					genotype['Gmax'] = config.getfloat("GENOTYPE",'gmax')
					found = 1
				elif o == "kdet":
					genotype['Kdet'] = config.getfloat("GENOTYPE","kdet")
					found = 1
				elif o == "kmc":
					genotype['KMC'] = config.getfloat("GENOTYPE","kmc")
					found = 1
				elif o == "zmc":
					genotype['ZMC'] = config.getfloat("GENOTYPE","zmc")
					found = 1
				if found == 0:
					print("\033[1mOption %s in section [GENOTYPE] of \
						%s file is unknown\033[0m" % (o, genfile))
					sys.exit()
			return genotype
		else:
			raise Exception('File %s not found' % (genfile))

	def update(self, nutrient, precursor):
		change = 0
		for i in range(self.npops):
			change += self.subpops[i].update(
				nutrient[i]) if self.subpops[i].size > 0 else 0.0
			print('A', [sp.size for sp in self.subpops])
			print('A', [fp.size for fp in self.freepops])
			self.migration(precursor, i)
			print('B', [sp.size for sp in self.subpops])
			print('B', [fp.size for fp in self.freepops])
		return change

	def abiotic_update(self):
		"""
		Simulate constant addition and removal of particles
		"""
		self.turnover(*[self.subpops,
			self.freepops])

	def turnover(self, *populations):
		"""
		Remove first population, then add an empty instance in the last position
		"""
		for pop in populations:
			subpop = pop.pop(0)
			subpop.size = 0
			pop.append(subpop)

	def migration(self, precursor, pid):
		attachment, detachment = self.subpops[pid].flow(
			precursor[pid], pid)
		## Add detaching cells to free populations
		if int(detachment) > 0 and pid < self.npops - 1:
			self.dispersal(int(detachment), pid)
		## Remove attaching cells from adjacent free population
		if int(attachment) > 0:
			incoming = self.freepops[pid].del_cells(int(attachment))
			## Add attaching cells to current particle population
			if int(incoming) > 0:
				self.subpops[pid].size += int(incoming)
				self.subpops[pid].incoming -= int(incoming)

	def dispersal(self, ncells, current_id=0):
		cell_dist = [self.cell_dist.pmf(i)
			for i in range(1, self.npops - current_id)]
		## Normalize so that sum of frequencies is 1
		cell_dist = list(cell_dist/np.sum(cell_dist))
		added = 0
		for i in range(len(cell_dist)):
			particle_cells = round(ncells * cell_dist[i])
			self.freepops[current_id + i + 1].size += particle_cells
			added += particle_cells
			if added > 0:
				print('det', current_id, 'to', current_id + i + 1)
		while added < ncells:
			for r in range(current_id + 1, self.npops):
				self.freepops[r].add_cells(1)
				print('det', current_id, 'to', r)
				added += 1
				if added >= ncells:
					break


class Subpopulation(Population):
	"""
	A group of cells within a population
	"""
	def __init__(self, ncells):
		## Default is detached cell
		self.size = ncells

	def del_cells(self, old_cells):
		"""
		Remove cells from subpopulation
		"""
		ndels = np.min([old_cells, self.size])
		self.size -= ndels
		return ndels

	def Monod(self, nutconc):
		"""
		Estimate biomass growth rate according to a Monod growth model
		"""
		return self.size * nutconc


class ParticlePop(Subpopulation):
	"""
	A subpopulation of cells the live associated to a particle
	"""
	def __init__(self, ncells):
		super().__init__(ncells)
		self.biomass = 0.0
		self.incoming = 0.0
		self.outgoing = 0.0

	@property
	def delta_t(self):
		return self.mother_pop.delta_t

	def flow(self, conc, pid):
		## Estimate the rate at which cells encounter particle
		self.incoming += self.mother_pop.freepops[pid].size * (
			self.estimate_encounter_rate() * self.delta_t)
		if self.size > 0 and pid < self.mother_pop.npops - 1:
			## Need at least one cell to detach
			## Assume cells in the ``surface`` cannot detach
			if conc > 0:
				## Estimate rate at which cells deatch from particle
				self.outgoing += self.crit(conc) * self.delta_t
			else:
				## If particle exhausted, all cells detach
				self.outgoing += self.size * self.delta_t
		incoming = self.incoming
		outgoing = self.outgoing
		if int(self.outgoing) > 0:
			## Remove detaching cells from particle population
			outgoing = self.del_cells(int(self.outgoing))
			self.outgoing -= int(outgoing)
		return incoming, outgoing

	def estimate_encounter_rate(self):
		"""
		Estimate particle encounter rate
		The larger fraction of space cells occupy,
		the more cells will encounter a particle
		The larger the particle, the more cells will find it
		Assume cells are less likely to stick to more
		highly-populated particles
		"""
		return self.mother_pop.beta * self.capacity/(
			self.size + 1)

	def update(self, nutconc):
		"""
		Estimate population growth under current environmental conditions
		"""
		oldsize = self.size
		delta_biomass = self.Monod(nutconc)
		new_biomass = delta_biomass * self.delta_t
		self.biomass += new_biomass
		if int(self.biomass) > 0:
			self.size += int(self.biomass)
			self.biomass -= int(self.biomass)
		change = np.abs(oldsize - self.size)
		return change
		
class FreePop(Subpopulation):
	"""
	A population of free-living cells that can find a particle
	"""
	def __init__(self):
		super().__init__(0)

	def update(self):
		"""
		Free-living populations cannot grow. Dummy method returns 0 growth if its called
		"""
		return 0.0