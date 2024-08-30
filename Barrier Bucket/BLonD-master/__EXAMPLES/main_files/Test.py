import numpy as np
from Functions_01 import fourier_coefficient, sigma_modulation, Voltage_amplitude, Barrier_Bucket

C_ring = 503.04
gamma = 5
gamma_transition = 7.087
phi_rf = 0


eta = -1
N = 25 # amount of Fourirer

harm = 2
Ts = 1701
T_1 = 50 # ns
T_2 = (Ts - 2*T_1)/(2*harm) # ns
phi_r = (2*T_1) / T_2
h_r = np.pi/phi_r
m = 3
V_peak = 5000
U = 300
N = 25

h, V, phi = Barrier_Bucket (V_peak, T_1, harm, C_ring, gamma, gamma_transition, phi_rf, N)
h_acc, V_acc, phi_acc = Barrier_Bucket (U, Ts, harm, C_ring, gamma, gamma_transition, phi_rf, N)

#print (h)
#print (V)
#print (V_acc)
#print (np.add(V, V_acc))
#print (phi)

a0 = np.linspace(0.0196684, 0.0201662, 6225+1)
a1 = np.linspace(-0.0608641, -0.0669929, 6225+1)
a2 = np.linspace(0.843772, 1.09957, 6225+1)

btot = []
b=5
for i in range (0, 5):
    btot.append(b)
    i = i + 1
b=5
for i in range (0, 2):
    btot.append(-b)
    i = i + 1

#aaaa = np.linspace(2e9, 2e9+(300*2000), 2001)
#print ("a0 = ", a0)
#print ("a1 = ", a1)
#print ("a2 = ", a2)

print ("V1+V2 = ", btot)
