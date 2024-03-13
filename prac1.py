# Import the required libraries

import matplotlib.pyplot as plt
import random as r
import math

# Define constants
Cm = 1.0
V_Na = 50.0
V_K = -77.0
V_L = -54.387
g_Na = 120.0
g_K = 36.0
g_L = 0.3

# Define initial conditions
V = -65.0
m = 0.053
h = 0.6
n = 0.318
I_ext = 0

# Define time step 
dt = 0.001

# Define Hodgkin-Huxley equations
def dVdt(t):
	return (1 / Cm) * (I(t) - (g_L * (V - V_L)) - g_Na * m**3 * h * (V - V_Na) - g_K * n**4 * (V - V_K))

def dmdt():
	return (alpha_m() * (1 - m)) - (beta_m() * m) 

def dhdt():
	return (alpha_h() * (1 - h)) - (beta_h() * h) 

def dndt():
	return (alpha_n() * (1 - n)) - (beta_n() * n) 

def alpha_m():
	return (0.1 * (-V - 40)) / (math.exp((-V - 40) / 10) - 1) 

def alpha_h():
	return 0.07 * math.exp((-V - 65.) / 18.)

def alpha_n():
	return (0.01 * (-V - 55)) / (math.exp((-V - 55) / 10) -1)

def beta_m():
	return 4 * math.exp((-V - 65) / 18)

def beta_h():
	return 1 / (math.exp((-V - 35) / 10) + 1)

def beta_n():
	return 0.125 * math.exp((-V - 65) / 80)

def I(t):
	if t < 0.001:
		return  0.1 / (7.854 * 10**(-3))
	else:
		return 0
	
# Execute functions
t = 0
points = {}
mval = []
hval = []
nval = []
tval = []
while t < 100:
	V += dt * dVdt(t)
	m += dt * dmdt()
	h += dt * dhdt()
	n += dt * dndt()
	t += dt
	tval.append(t)
	mval.append(m)
	nval.append(n)
	hval.append(h)
	points[t] = V
plt.plot(list(points.keys()), list(points.values()), color = "black", linewidth = 0.7)
plt.xlabel('t')
plt.ylabel('V')
#plt.plot(tval, mval, linewidth = 0.7, label = 'm')
#plt.plot(tval, nval, linewidth = 0.7, label = 'n')
#plt.plot(tval, hval, linewidth = 0.7, label = 'h')
#plt.legend()
plt.show()


	
