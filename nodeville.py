#!/usr/bin/env python

# Steven Herbst
# sgherbst@gmail.com

# January 25, 2017

import numpy as np
import matplotlib.pyplot as plt

def main():

	c = Circuit()

	n1 = c.node('n1')
	n2 = c.node('n2')
	n3 = c.node('n3')
	gnd = GndNode()

	V1 = c.add(VoltSource(n1, gnd, 1.0))

	R1 = c.add(Resistor(n1, n2, 100))
	R2 = c.add(Resistor(n2, n3, 100))
	R3 = c.add(Resistor(n3, gnd, 100))

	C1 = c.add(Capacitor(n1, gnd, 1e-9))
	C2 = c.add(Capacitor(n2, gnd, 1e-9))
	C3 = c.add(Capacitor(n3, gnd, 1e-9))
	

	c.run(1e-6)
	c.plot(c.nodes)
	c.plot_step()
	c.plot_it()

	plt.show()

class Circuit:

	def __init__(self):
		self.nodes = []
		self.devices = []

	def node(self, name):
		self.nodes.append(Node(name))
		return self.nodes[-1]

	def add(self, device):
		self.devices.append(device)
		return device

	def run(self, tstop):
		# create empty lists to hold time and node voltages
		self.tvec = []
		self.svec = []
		self.ivec = []
		self.vvec = {}
		for node in self.nodes:
			self.vvec[node.name] = []

		# initialize loop variables
		t = 0
		dt = 1e-12
		dt_min = 1e-18
		dt_shrink_factor = 0.5
		dt_grow_factor = 1.01
		itmax = 100

		while t<tstop:
			it = 0

			# start time step for nodes and devices
			for node in self.nodes:
				node.start_step()
			for device in self.devices:
				device.start_step()

			# iteratively solve nodes at this timestep
			while it<itmax:

				# start iteration for nodes and devices
				for node in self.nodes:
					node.start_iter()
				for device in self.devices:
					device.start_iter()

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
				dt = dt*dt_shrink_factor
				if dt<dt_min:
					raise Exception('Timestep too small.')
				else:
					continue
			
			# otherwise accept the timestep
			for node in self.nodes:
				node.end_step()
			for device in self.devices:
				device.end_step()

			# update timestep
			self.tvec.append(t)
			self.svec.append(dt)
			self.ivec.append(it)
			t = t+dt
			dt = dt*dt_grow_factor

			# then save all the node voltages
			for node in self.nodes:
				self.vvec[node.name].append(node.volt)

		# convert time and voltage lists to arrays
		print len(self.tvec)
		self.tvec = np.asarray(self.tvec)
		self.svec = np.asarray(self.svec)
		self.ivec = np.asarray(self.ivec)
		for node in self.nodes:
			self.vvec[node.name] = np.asarray(self.vvec[node.name])

	def plot(self, nodes):
		fig1 = plt.figure()
		ax1 = fig1.add_subplot(111)

		styles = ['r-', 'g-', 'b-', 'c-', 'm-', 'k-', 'y-']
		for node, style in zip(nodes,styles):
			ax1.plot(self.tvec, self.vvec[node.name], style, label=node.name)

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

class GndNode:
	def __init__(self):
		self.volt = 0.0
		self.nextVolt = 0.0

	def add(self, value):
		pass

class Node:
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

	def end_step(self):
		self.volt = self.nextVolt

class Resistor:
	def __init__(self, p, n, res=1e3, imax=10):
		self.p = p
		self.n = n
		self.res = res
		self.imax = imax

	def calc_volt(self):
		return float(self.p.nextVolt-self.n.nextVolt)

	def calc_curr(self, dt):
		return self.calc_volt()/self.res

	def start_step(self):
		pass

	def start_iter(self):
		pass

	def update(self, dt):
		curr = self.calc_curr(dt)

		self.p.add(-curr)
		self.n.add(curr)

	def end_step(self):
		pass

	def flag(self, dt):
		return abs(self.calc_curr(dt))>self.imax

class Capacitor:
	def __init__(self, p, n, cap=1e-9, imax=10):
		self.p = p
		self.n = n
		self.cap = cap
		self.imax = imax
		self.vPrev = 0.0

	def calc_volt(self):
		return float(self.p.nextVolt-self.n.nextVolt)

	def calc_curr(self, dt):
		return self.cap*(self.calc_volt()-self.vPrev)/dt

	def start_step(self):
		pass

	def start_iter(self):
		pass

	def update(self, dt):
		# calculate current
		curr = self.calc_curr(dt)

		# update nodes
		self.p.add(-curr)
		self.n.add(curr)

	def end_step(self):
		self.vPrev = self.calc_volt()

	def flag(self, dt):
		return abs(self.calc_curr(dt))>self.imax

class VoltSource:
	def __init__(self, p, n, volt=1.0, res=1):
		self.p = p
		self.n = n
		self.volt = volt
		self.res = res
		self.imax = 10.0*volt/res

	def calc_volt(self):
		return float(self.p.nextVolt-self.n.nextVolt)

	def calc_curr(self, dt):
		return float(self.volt-self.calc_volt())/self.res

	def start_step(self):
		pass

	def start_iter(self):
		pass

	def update(self, dt):
		curr = self.calc_curr(dt)

		self.p.add(curr)
		self.n.add(-curr)

	def end_step(self):
		pass

	def flag(self, dt):
		return abs(self.calc_curr(dt))>self.imax

if __name__ == "__main__":
	main()