from ctypes import *
import numpy
import matplotlib.pyplot as plt
from distpyinterface import *
class DISTlib:
	def __init__(self):
		self.dist = cdll.LoadLibrary("../demo/buildDemo/libhello.so")

	def initializedistribution(self,num):
		self.dist.initializedistribution(c_int(num))

	def setEmittance12(self, e1, e2):
		self.dist.setemitt12(c_double(e1),c_double(e2))

	def setEmittance3(self,e3):
		self.dist.setemitt3(c_double(e3))

	def setscan_para_diagonal(self,variable, variable_type, space_type, start_value, stop_value):
		self.dist.setscan_para_diagonal(c_int(variable), c_int(variable_type), c_int(space_type), c_double(start_value), c_double(stop_value))

	def settotalsteps(self, totalsetps):
		self.dist.settotalsteps(c_int(totalsetps))

	def setEandMass(self,energy0, mass0):
		self.dist.sete0andmass0(c_double(energy0),c_double(mass0))

	def  get6trackcoord(self):
		totlength = c_int(0)
		npart=pointer(totlength)

		self.dist.getarraylength(npart)
		print(npart.contents.value)

		double6 = c_double * npart.contents.value
		x1 =  double6(0)
		x2 =  double6(0)
		x3 =  double6(0)
		x4 =  double6(0)
		x5 =  double6(0)
		x6 =  double6(0)
		xd = []
		pxd = []
		yd = []
		yxd = []
		deltas = []
		dp = []
		self.dist.get6trackcoord(x1,x2,x3,x4,x5,x6, npart)
		for i in range(0,npart.contents.value):
			xd.append(x1[i])
			pxd.append(x2[i])
			yd.append(x3[i])
			yxd.append(x4[i])
			deltas.append(x5[i])
			dp.append(x6[i])
		return [xd, pxd, yd, yxd, deltas, dp]


	def gettasunitcoord(self):
		totlength = c_int(0)
		npart=pointer(totlength)

		self.dist.getarraylength(npart)
		print(npart.contents.value)

		double6 = c_double * npart.contents.value
		x1 =  double6(0)
		x2 =  double6(0)
		x3 =  double6(0)
		x4 =  double6(0)
		x5 =  double6(0)
		x6 =  double6(0)
		xd = []
		pxd = []
		yd = []
		yxd = []
		deltas = []
		dp = []
		self.dist.getunconvertedcoord(x1,x2,x3,x4,x5,x6, npart)
		for i in range(0,npart.contents.value):
			xd.append(x1[i])
			pxd.append(x2[i])
			yd.append(x3[i])
			yxd.append(x4[i])
			deltas.append(x5[i])
			dp.append(x6[i])
		return [xd, pxd, yd, yxd, deltas, dp]


	def settasmatrix(self,  tas):
		for i in range(0,6):
			for j in range(0,6):
				self.dist.settasmatrix_element(c_double(tas[i][j]), c_int(i), c_int(j))


	def set_action_angle_cuts(self,variable, min_v, max_v):
		self.dist.setactionanglecut(c_int(variable), c_double(min_v), c_double(max_v))

			#def createtas0coupling(dist, betx, alfx, bety, alfy, dx, dpx, dy, dpy):
#	dist.createtas0coupling_(c_double(betx),c_double(alfx),c_double(bety),c_double(alfy), c_double(0), c_double(0), c_double(0), c_double(0))