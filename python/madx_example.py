from distpyinterface import *
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
import scipy.linalg as la
def readtasfrommadx(filename):
	inputn = np.loadtxt(filename)
	print(inputn)
	return inputn

e1 = 1
e2 = 1.2
e3 = 0.2
zero = 0
one =1
pia2=np.pi*2
dist = DISTlib()

dist.initializedistribution(2)
dist.setEmittance12(e1,e2)
dist.setEmittance3(e3)
dist.setEandMass(6500, 0.938)
dist.settotalsteps(3000)

myfile = "/home/tobias/codes/SixTrackTobias/test/orbit6d-element-quadrupole/myeigen.dat"
tas = readtasfrommadx(myfile)


print(tas)
#1st First input (coordinate)
#2nd coord type, 0 - action angle, 1 - normalized
#3rd type 4 - uniform dist, 5 - normalized,  6-reyligh
#4 start value for uniform, mean value for normalized and reyligh
#5 sstop value for uniform, sigma for normalized and reyligh
dist.settasmatrix(tas)
dist.setscan_para_diagonal(0,0,6,zero,one); 
dist.setscan_para_diagonal(1,0,4,zero,pia2); 
dist.setscan_para_diagonal(2,0,6,zero,one);
dist.setscan_para_diagonal(3,0,4,zero,pia2);
dist.setscan_para_diagonal(4,0,6,zero,one);
dist.setscan_para_diagonal(5,0,4,zero,pia2);

[x,xp,y,yp,sigma,dp]=dist.gettasunitcoord()
#print(np.std(x), np.std(xp),np.std(y), np.std(yp), np.std(sigma), np.std(dp))
#print(invtas)
#print(inv(invtas))
#print(tas)
print(x[-1],xp[-1],y[-1],yp[-1],sigma[-1],dp[-1])
plt.plot(x,xp , '.')
plt.show()