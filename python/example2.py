import numpy as np
from ctypes import *
import matplotlib.pyplot as plt
from distpyinterface import *
def readtasfromdump(filename):
	file = open(filename, "r")
	splitlines = (file.readlines())
	tasstr = (str.split(splitlines[3]))
	tasstr = np.array(tasstr[3:])
	tas = tasstr.astype(np.float)
	tasm = tas.reshape(6,6)
	return tasm

def settasmatrix(dist, tas):
	for i in range(0,6):
		for j in range(0,6):
			dist.settasmatrix_element(c_double(tas[i][j]), c_int(i), c_int(j))

def setParameters(dist, index, start, stop, numb, type):
	dist.setparameter_(byref(c_int(index)), byref(c_double(start)),byref(c_double(stop)),byref(c_int(numb)),byref(c_int(type)))

'''
SIMULATION
#/                   mass, charge,   A,    Z, momentum (7Z TeV/c)
REFERENCE 938.0      1  1    1  6500000.0 
/      turns, total particles
TRACK  100000   64
NEXT
/
DISTRIBUTION
//INPUT PARAMETERS
EMITTANCE   1           2.5e-6    //Normalized [m]
EMITTANCE   2           2.503e-6  //Normalized [m]
/          BETAX           ALFX            DX            DPX       BETAY          ALFY            DY              DPY
TWISS      156.921221293   2.11599999399   0.2973032402 -0.0030478 78.0330611129 -1.08169977617  -0.1155854323 -0.0011203
CLOSEDORBIT  0.0   0.0   0.0   0.0  0.0   0.0
DUMP FILENAME BINARY NORMALIZED
// DISTRIBUTION TYPE
SUBSAMPLE  1
JX       LINEAR    2.0   4.0 / [sigma]
PHIX     CONSTANT  0
JY       LINEAR    2.0   4.0 / [sigma]
PHIX     CONSTANT  0
DELTA    CONSTANT  0.00027
SIGMA    CONSTANT  0
SUBSAMPLE  2 /second pair
JX       LINEAR    2.0   4.0 / [sigma]
PHIX     CONSTANT  0.001
JY       LINEAR    2.0   4.0 / [sigma]
PHIX     CONSTANT  0
DELTA    CONSTANT  0.00027
SIGMA    CONSTANT  0
NEXT
DUMP
POSTPR 1000 / produce fort.90 every 1000 turns
NEXT
'''

dist = DISTlib()
dist.initializedistribution(2)

myfile = "/mnt/c/Users/tobia/codes/SixTrackTobias/test/orbit6d-element-quadrupole/IP1_DUMP.dat"
tas = readtasfromdump(myfile)
#tas[4][2]=0.002
#tas[5][2]=0.1
dist.settasmatrix(tas)




pia2 = np.pi*2
momentum = (6500000)
mass = 938.0
eps = 2.0
dim = (6)
zero = (0)
e1 = (1.0)
e2 = (1.0)

dist.setEmittance12(e1,e2)
dist.setEmittance3(1)

dist.setEandMass(momentum, mass)

dist.settotalsteps(3000)
#dist.setdisttype(0)


dist.settasmatrix(tas)
dist.setscan_para_diagonal(0,0,6,zero,1);
dist.setscan_para_diagonal(1,0,4,zero,pia2);
dist.setscan_para_diagonal(2,0,6,zero,1);
dist.setscan_para_diagonal(3,0,4,zero,pia2);
dist.setscan_para_diagonal(4,3,0,1.00,0.1);
dist.setscan_para_diagonal(5,3,0,0.0,0.02);

[x,xp,y,yp,sigma,dp]=dist.gettasunitcoord()

print(x[-1],xp[-1],y[-1],yp[-1],sigma[-1],dp[-1])
plt.plot(x,xp , '.')
plt.show()


#dist.printdistsettings_()
#dist.print2fort13()

