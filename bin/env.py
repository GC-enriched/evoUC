#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
	Author: Angel Garcia-Casillas del Riego (angel.garcia-casillasdelriego@uzh.ch)
"""

class Environment(object):
	"""
	A solution of precursor molecule, patchly distributed in highly concentrated particles
	"""
	def __init__(self, env_params):
		for param in env_params:
			self.__setattr__(param, env_params[param])
		## Number of particles constant
		self.init_env(env_params)

	@property
	def particles(self):
		return self._particles

	@particles.setter
	def particles(self, nparticles):
		self._particles = [Particle(self.C0,
			self.Kprod) for n in range(nparticles)]

	def init_env(self, env_params):
		self.particles = self.nparticles
		for particle in self.particles:
			for param in env_params:
				particle.__setattr__(param, env_params[param])
			particle.C0 = self.C0/self.Ks
			particle.locs = self.Ks
			particle.nlocs = len(particle.locs)

	def biotic_update(self, ncells):
		"""
		Estimate population growth rate under a Monod model
		given nutrient concentrations,
		then update the number of cells applying Euler method
		"""
		private_nut = [self.particles[i].estimate_public_good(
			ncells[i]) if ncells[i] > 0 else 0.0
			for i in range(self.nparticles)]
		precursor_nut = [self.particles[i].conc
			for i in range(self.nparticles)]
		return private_nut, precursor_nut

	def abiotic_update(self):
		"""
		Simulate constant addition and removal of particles
		"""
		particle = self.particles.pop(0)
		particle.conc = self.C0
		self.particles.append(particle)

class Particle(Environment):
	"""
	A highly concentrated particle of precursor molecule
	"""
	def __init__(self, C0, Kprod):
		self.conc = C0
		self.release = lambda x: x/(Kprod + x)

	@property
	def locs(self):
		return self._locs

	@locs.setter
	def locs(self, capacity):
		self._locs = [self.C0 for i in range(capacity)]

	def estimate_public_good(self, ncells):
		production = 0.0
		for i in range(self.nlocs):
			cell_mining = self.release(self.locs[i])
			self.locs[i] = max([self.locs[i] - cell_mining, 0.0])
			self.conc = max([self.conc - cell_mining, 0.0])
			if self.locs[i] < 0.0:
				self.locs[i] = 0.0
			production += cell_mining
		return self.calculate_diffusion_efficiency(
			production, ncells)

	def calculate_diffusion_efficiency(self, nutconc, ncells):
		diff_fraction = self.Kdiff/(ncells + self.Kdiff)
		return nutconc * diff_fraction