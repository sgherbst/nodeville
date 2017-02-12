#!/usr/bin/env python

# Steven Herbst
# sgherbst@gmail.com

# January 25, 2017

import numpy as np
import math
import matplotlib.pyplot as plt
import random
import numpy as np

def main():

	c = Circuit()

	n1 = c.add(Node('n1', cap=1e-6))
	n2 = c.add(Node('n2'))
	gnd = GndNode()

	V1 = c.add(VoltSource(n1, gnd, 2.0))
	R1 = c.add(Resistor(n1, n2, 1e3))
	R2 = c.add(Resistor(n2, gnd, 1e3))

	c.run(75e-9)
	c.plot(c.nodes)
	c.plot_step()

	# show all plots
	plt.show()

class Circuit(object):

	def __init__(self):
		self.nodes = []
		self.devices = []

		# simulation parameters

		self.dt_def_subdiv = 100
		self.dt_min = 1e-18

		self.dt_safety_factor = 0.9
		self.dt_grow_max = 2.0
		self.dt_shrink_min = 0.3

	def add(self, thing):
		if isinstance(thing, Device):
			self.devices.append(thing)
		elif isinstance(thing, Node):
			self.nodes.append(thing)
		else:
			raise Exception('Trying to add object of unknown type to circuit.')
		return thing

	def run(self, tstop):
		# time and voltage storage
		self.t_vec = []
		self.volt_vec = {}
		for node in self.nodes:
			self.volt_vec[node.name] = []

		# keep track of loop parameters for debugging purposes
		self.dt_vec = [] 

		# initialize loop variables
		t = 0
		dt = float(tstop)/self.dt_def_subdiv

		while t<tstop:

			# three time steps
			# n=0: full step
			# n=1: rewind, half-step
			# n=2: half-step
			for n in range(3):
				if n==0:
					t_iter = t
					dt_iter = dt
				elif n==1:
					t_iter = t
					dt_iter = 0.5*dt
				else:
					t_iter = t+0.5*dt
					dt_iter = 0.5*dt

				for node in self.nodes:
					node.clear(n)
				# update device currents
				for device in self.devices:
					device.update(t_iter, dt_iter)
				# update node voltages
				for node in self.nodes:
					node.update(n, dt_iter)
			
			# check if timestep is OK
			min_ratio = float('inf')
			for node in self.nodes:
				min_ratio = min(node.get_error_ratio(), min_ratio)
			if min_ratio < 1.0:
				dt = dt*self.dt_safety_factor*max(min_ratio, self.dt_shrink_min)
				if dt >= self.dt_min:
					continue
				else:
					raise Exception('Timestep too small.')

			# if it is then accept the timestep
			for device in self.devices:
				device.accept(t, dt)
			for node in self.nodes:
				node.accept()

			# update time
			t = t+dt
			self.t_vec.append(t)

			# update timestep
			self.dt_vec.append(dt)
			dt = dt*self.dt_safety_factor*min(min_ratio, self.dt_grow_max)

			# then save all the node voltages
			for node in self.nodes:
				self.volt_vec[node.name].append(node.volt)

		# convert time and voltage lists to arrays
		print len(self.t_vec)
		self.t_vec = np.asarray(self.t_vec)
		self.dt_vec = np.asarray(self.dt_vec)
		for node in self.nodes:
			self.volt_vec[node.name] = np.asarray(self.volt_vec[node.name])

	def plot(self, nodes):
		fig1 = plt.figure()
		ax1 = fig1.add_subplot(111)

		styles = ['r-', 'g-', 'b-', 'c-', 'm-', 'k-', 'y-']
		for node, style in zip(nodes,styles):
			ax1.plot(self.t_vec, self.volt_vec[node.name], style, label=node.name)

		ax1.set_xlabel('Time (s)')
		ax1.set_ylabel('Voltage (V)')
		ax1.set_title('Transient simulation')
		ax1.legend(loc='lower right')

	def plot_step(self):
		fig1 = plt.figure()
		ax1 = fig1.add_subplot(111)

		ax1.step(self.t_vec, self.dt_vec, '-r', label='dt')

		ax1.set_xlabel('Time (s)')
		ax1.set_ylabel('Timestep (s)')
		ax1.set_title('Transient simulation')
		ax1.legend(loc='lower right')

class GndNode(object):
	def __init__(self):
		self.volt = 0.0
		self.nextVolt = 0.0

	def add(self, value):
		pass

class Node(object):
	def __init__(self, name, cap=1e-12, vtol=1e-3):
		# store settings
		self.name = name
		self.cap = cap
		self.vtol = vtol

		# initialize voltages
		self.volt = 0
		self.nextVolt = 0

	def clear(self, n):
		if n==0:
			self.nextVolt = self.volt
			self.curr = 0.0
		elif n==1:
			self.nextVolt = self.volt
			self.curr = 0.0
		else:
			self.curr = 0.0

	def add(self, value):
		self.curr = self.curr + value

	def update(self, n, dt):
		if n==0:
			self.nextVolt = self.volt + self.curr*dt/self.cap
			self.y0 = self.nextVolt
		elif n==1:
			self.nextVolt = self.volt + self.curr*dt/self.cap
		else:
			self.nextVolt = self.nextVolt + self.curr*dt/self.cap
			self.y1 = self.nextVolt

	def get_error_ratio(self):
		if self.y1 == self.y0:
			return float('inf')
		else:
			return self.vtol/abs(self.y1-self.y0)

	def accept(self):
		self.volt = 2*self.y1 - self.y0

class Device(object):
	pass

class TwoTerm(Device):
	def __init__(self, p, n):
		self.p = p
		self.n = n

		self.vPrev = 0.0
		self.iPrev = 0.0

	def get_volt(self):
		return float(self.p.nextVolt-self.n.nextVolt)

	def update(self, t, dt):
		curr = self.calc_curr(t, dt)

		self.p.add(-curr)
		self.n.add(curr)

	def accept(self, t, dt):
		self.vPrev = self.get_volt()
		self.iPrev = self.calc_curr(t, dt)

class Resistor(TwoTerm):
	def __init__(self, p, n, res=1e3):
		self.res = res
		super(Resistor, self).__init__(p,n)

	def calc_curr(self, t, dt):
		return self.get_volt()/self.res

class Capacitor(TwoTerm):
	def __init__(self, p, n, cap=1e-9):
		self.cap = cap
		super(Capacitor, self).__init__(p,n)

	def calc_curr(self, t, dt):
		return self.cap*(self.get_volt()-self.vPrev)/dt

class Inductor(TwoTerm):
	def __init__(self, p, n, ind=100e-9):
		self.ind = ind
		super(Inductor, self).__init__(p,n)

	def calc_curr(self, t, dt):
		return self.iPrev + (1.0/self.ind)*self.get_volt()*dt

class VoltSource(TwoTerm):
	def __init__(self, p, n, volt=1.0, res=10e-3, tr=10e-9):
		self.volt = volt
		self.res = res
		self.tr = tr
		super(VoltSource, self).__init__(p,n)

	def calc_curr(self, t, dt):
		if t <= self.tr:
			vset = float(t)/self.tr * self.volt
		else:
			vset = self.volt

		return float(self.get_volt()-vset)/self.res

class Diode(TwoTerm):
	def __init__(self, p, n, Is=1e-15, Vlin=0.85, Vth=0.025):
		# save model parameters
		self.Is = Is
		self.Vth = Vth
		self.Vlin = Vlin

		# setup the limited exponential model
		self.setup_limexp()

		super(Diode, self).__init__(p,n)

	def raw_current(self, v):
		# Ideal diode behavior.  Should not be used directly
		# due to convergence issues at high bias
		return self.Is*(math.exp(v/self.Vth)-1.0)

	def setup_limexp(self):
		# Setup model that gives linear I-V behavior when
		# the applied voltage exceeds Vlin
		self.Ilin = self.raw_current(self.Vlin)
		self.Glin = self.Ilin/self.Vlin

	def calc_curr(self, t, dt):
		# Implement limited exponential behavior
		volt = self.get_volt()
		if volt >= self.Vlin:
			return self.Ilin + (volt-self.Vlin)*self.Glin
		else:
			return self.raw_current(volt)

if __name__ == "__main__":
	main()