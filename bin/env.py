#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
	Author: Angel Garcia-Casillas del Riego (angel.garcia-casillasdelriego@uzh.ch)
"""

class Environment(object):
	"""
	A solution of precursor molecule, patchly distributed in highly concentrated particles
	"""
	def __init__(self, nparticles, env_params):
		for param in env_params:
			self.__setattr__(param, env_params[param])
		self.nparticles = nparticles
		self.init_env()

	@property
	def particles(self):
		return self._particles

	@particles.setter
	def particles(self, nparticles):
		self._particles = {n:Particle() for n in range(nparticles)}

	def init_env(self):
		self.particles = self.nparticles

class Particle(Environment):
	"""
	A highly concentrated particle of precursor molecule
	"""
	def __init__(self):
		self.conc = self.C0
