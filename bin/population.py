#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
	Author: Angel Garcia-Casillas del Riego (angel.garcia-casillasdelriego@uzh.ch)
"""

class Population(object):
	"""
	A population of cells that grows by consuming
	a limiting nutrient following a Monod model
	"""
	def __init__(self, Kgrowth, Gmax, Kdet, time_inc):
		self.size = 1
		self.Kgrowth = Kgrowth
		self.Gmax = Gmax
		self.Kdet = Kdet
		self.time_inc = time_inc

	def update(self, nutconc):
		"""
		Estimate population growth rate under a Monod model
		given nutrient concentrations,
		then update the number of cells applying Euler method
		"""
		growth_rate = self.get_Monod_growth(nutconc)
		self.size += int(self.size * growth_rate / self.time_inc)

	def get_Monod_growth(self, nutconc):
		"""
		Apply Monod's equation to estimate population growth under
		current environmental composition
		"""
		growth_rate = self.Gmax * self.Kgrowth/(self.Kgrowth + nutconc)
		return growth_rate