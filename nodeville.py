#!/usr/bin/env python

# Steven Herbst
# sgherbst@gmail.com

# January 25, 2017

import numpy as np
import math
import matplotlib.pyplot as plt

def main():

	c = Circuit()

	n1 = c.node('n1')
	n2 = c.node('n2')
	n3 = c.node('n3')
	gnd = GndNode()

	V1 = c.add(VoltSource(n1, gnd, 10.0))

	R1 = c.add(Resistor(n1, n2, 100))
	R2 = c.add(Resistor(n2, n3, 100))
	R3 = c.add(Resistor(n3, gnd, 100))

	C1 = c.add(Capacitor(n1, gnd, 1e-9))
	C2 = c.add(Capacitor(n2, gnd, 1e-9))
	C3 = c.add(Capacitor(n3, gnd, 1e-9))

	D1 = c.add(Diode(n3, gnd))

	c.run(1e-6)
	c.plot(c.nodes)
	# c.plot_step()
	# c.plot_it()

	# show all plots
	plt.show()

class Circuit(object):

	def __init__(self):
		self.nodes = []
		self.devices = []

		# simulation parameters

		self.dt_def = 1e-12
		self.dt_min = 1e-18
		self.dt_max = 1

		self.dt_shrink_factor = 0.5
		self. dt_grow_factor = 1.01

		self.it_max = 100

	def node(self, name):
		self.nodes.append(Node(name))
		return self.nodes[-1]

	def add(self, device):
		self.devices.append(device)
		return device

	def run(self, tstop):
		# time and voltage storage
		self.t_vec = []
		self.volt_vec = {}
		for node in self.nodes:
			self.volt_vec[node.name] = []

		# keep track of loop parameters for debugging purposes
		self.it_vec = []
		self.dt_vec = [] 

		# initialize loop variables
		t = 0
		dt = self.dt_def

		while t<tstop:
			it = 0

			# start time step for nodes and devices
			for node in self.nodes:
				node.start_step()

			# iteratively solve nodes at this timestep
			while it<self.it_max:

				# start iteration for nodes and devices
				for node in self.nodes:
					node.start_iter()

				flag = False

				# device will raise a flag if currents exceed a certain amount
				# or changes more than a certain amount
				for device in self.devices:
					device.update(dt)
					flag = flag or device.flag(dt)

				# node will raise a flag if voltage change exceeds allowed range
				# or if the current error exceeds allowed range
				for node in self.nodes:
					node.update(dt)
					flag = flag or node.flag(dt)

				# if no flags were raised, break out of the loop
				if not flag:
					break
				else:
					it=it+1
			
			# if there's a flag on exit, reduce the timestep and try again
			if flag:
				dt = dt*self.dt_shrink_factor
				if dt < self.dt_min:
					raise Exception('Timestep too small.')
				else:
					continue
			
			# otherwise accept the timestep
			for node in self.nodes:
				node.end_step(dt)
			for device in self.devices:
				device.end_step(dt)

			# update timestep
			self.t_vec.append(t)
			t = t+dt
			dt = min(dt*self.dt_grow_factor, self.dt_max)

			# store dt and iteration count for debugging purposes
			self.dt_vec.append(dt)
			self.it_vec.append(it)

			# then save all the node voltages
			for node in self.nodes:
				self.volt_vec[node.name].append(node.volt)

		# convert time and voltage lists to arrays
		print len(self.t_vec)
		self.t_vec = np.asarray(self.t_vec)
		self.dt_vec = np.asarray(self.dt_vec)
		self.it_vec = np.asarray(self.it_vec)
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

		ax1.plot(self.tvec, self.svec, '-r', label='dt')

		ax1.set_xlabel('Time (s)')
		ax1.set_ylabel('Timestep (s)')
		ax1.set_title('Transient simulation')
		ax1.legend(loc='lower right')

	def plot_it(self):
		fig1 = plt.figure()
		ax1 = fig1.add_subplot(111)

		ax1.plot(self.tvec, self.ivec, '-r', label='it')

		ax1.set_xlabel('Time (s)')
		ax1.set_ylabel('Iter count')
		ax1.set_title('Transient simulation')
		ax1.legend(loc='lower right')

class GndNode(object):
	def __init__(self):
		self.volt = 0.0
		self.nextVolt = 0.0

	def add(self, value):
		pass

class Node(object):
	def __init__(self, name, maxDelta=10e-3, maxErr=1e-7):
		self.name = name
		self.maxDelta = maxDelta
		self.maxErr = maxErr
		self.volt = 0
		self.nextVolt = 0

	def start_step(self):
		self.vmin = self.volt-self.maxDelta
		self.vmax = self.volt+self.maxDelta
		self.nextVolt = self.vmin
		self.count = 0

	def start_iter(self):
		self.curr = 0

	def add(self, value):
		self.curr = self.curr + value

	def update(self, dt):
		# first check current at max and min points
		# then do a binary search in between
		if self.count==0:
			self.i_vmin = self.curr
			self.nextVolt = self.vmax
			self.count = self.count + 1
			return
		elif self.count==1:
			self.i_vmax = self.curr
		else:
			if self.i_vmin < self.i_vmax:
				if self.curr < 0:
					self.i_vmin = self.curr
					self.vmin = self.nextVolt
				else:
					self.i_vmax = self.curr
					self.vmax = self.nextVolt
			else:
				if self.curr < 0:
					self.i_vmax = self.curr
					self.vmax = self.nextVolt
				else:
					self.i_vmin = self.curr
					self.vmin = self.nextVolt
		
		# Next search point in between min and max values
		self.nextVolt = 0.5*(self.vmin+self.vmax)
		self.count = self.count + 1

	def flag(self, dt):
		return abs(self.curr) > self.maxErr

	def end_step(self, dt):
		self.volt = self.nextVolt

class TwoTerm(object):
	def __init__(self, p, n, iMax=10):
		self.p = p
		self.n = n
		self.vPrev = 0.0
		self.iPrev = 0.0
		self.iMax = iMax

	def calc_volt(self, dt):
		return float(self.p.nextVolt-self.n.nextVolt)

	def update(self, dt):
		curr = self.calc_curr(dt)

		self.p.add(-curr)
		self.n.add(curr)

	def end_step(self, dt):
		self.vPrev = self.calc_volt(dt)
		self.iPrev = self.calc_curr(dt)

	def flag(self, dt):
		return abs(self.calc_curr(dt))>self.iMax

class Resistor(TwoTerm):
	def __init__(self, p, n, res=1e3):
		self.res = res
		super(Resistor, self).__init__(p,n)

	def calc_curr(self, dt):
		return self.calc_volt(dt)/self.res

class Capacitor(TwoTerm):
	def __init__(self, p, n, cap=1e-9):
		self.cap = cap
		super(Capacitor, self).__init__(p,n)

	def calc_curr(self, dt):
		return self.cap*(self.calc_volt(dt)-self.vPrev)/dt

class Inductor(TwoTerm):
	def __init__(self, p, n, ind=100e-9):
		self.ind = ind
		super(Inductor, self).__init__(p,n)

	def calc_curr(self, dt):
		return self.iPrev + (1.0/self.ind)*(self.calc_volt(dt)-self.vPrev)*dt

class VoltSource(TwoTerm):
	def __init__(self, p, n, volt=1.0, res=1):
		self.volt = volt
		self.res = res
		super(VoltSource, self).__init__(p,n)

	def calc_curr(self, dt):
		return float(self.calc_volt(dt)-self.volt)/self.res

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

	def calc_curr(self, dt):
		# Implement limited exponential behavior
		volt = self.calc_volt(dt)
		if volt >= self.Vlin:
			return self.Ilin + (volt-self.Vlin)*self.Glin
		else:
			return self.raw_current(volt)

if __name__ == "__main__":
	main()