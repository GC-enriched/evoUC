#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
	Author: Angel Garcia-Casillas del Riego (angel.garcia-casillasdelriego@uzh.ch)
"""

class Environment(object):
	"""
	A population of cells that grows by consuming
	a limiting nutrient following a Monod model
	"""
	def __init__(self, C0):
		self.C0 = C0
		self.C = C0
		self.Kprod = Kprod

	def update(self, popsize):
		"""
		Estimate population growth rate under a Monod model
		given nutrient concentrations,
		then update the number of cells applying Euler method
		"""
		private_nut = self.estimate_public_good(popsize)

	def estimate_public_good(self):
		production = self.C * self.Kprod
		return self.calculate_diffusion_efficiency(production)

	def calculate_diffusion_efficiency(self, nutconc):
		