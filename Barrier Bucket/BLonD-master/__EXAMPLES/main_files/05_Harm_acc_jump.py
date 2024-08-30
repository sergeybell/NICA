# Copyright 2014-2017 CERN. This software is distributed under the
# terms of the GNU General Public Licence version 3 (GPL Version 3),
# copied verbatim in the file LICENCE.md.
# In applying this licence, CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.
# Project website: http://blond.web.cern.ch/

'''
Tracking near transition and at transition
incliding high order MCF

:Authors: **Sergey Kolokolchikov**
:E-mail:sergey.bell13@gmail.com
'''

#  General Imports
from __future__ import division, print_function

import os
from builtins import range

import matplotlib as mpl
import numpy as np
import scipy.constants as const

from blond.beam.beam import Beam, Proton
from blond.beam.distributions import bigaussian, parabolic
from blond.beam.profile import CutOptions, FitOptions, Profile
from blond.input_parameters.rf_parameters import RFStation
#  BLonD Imports
from blond.input_parameters.ring import Ring
from blond.monitors.monitors import BunchMonitor
from blond.plots.plot import Plot
from blond.trackers.tracker import RingAndRFTracker

# My IMPORT
from blond.plots.plot_parameters import plot_voltage_programme
from Functions_01 import fourier_coefficient, sigma_modulation, Voltage_amplitude, Barrier_Bucket
from blond.beam.coasting_beam import generate_coasting_beam

mpl.use('Agg')

this_directory = os.path.dirname(os.path.realpath(__file__)) + '/'

os.makedirs(this_directory + '../output_files/05_Harm_acc_jump_fig', exist_ok=True)


# Simulation parameters -------------------------------------------------------
# Bunch parameters
N_b = 1.0e12                 # Design Intensity in SIS100
N_p = 50000                  # Macro-particles
#tau_0 = 1200.0e-9             # Initial bunch length, 4 sigma [s]

# Machine
m_p = 0.938e9          # mass of proton in eV
C_ring = 503.04        # Machine circumference [m]

#E_inj = 5.66704e9          # Injection Energy, [eV] before JUMP procedure
#E_inj = 2.0e9              # Injection Energy, [eV] Test
#E_inj = 1.87e9              # Injection Energy, [eV] from Nuclotron at low energy
#E_inj = 5.70832e9          # Injection Energy, [eV] before JUMP
E_inj = 4.7e9          # Injection Energy, [eV] before JUMP
gamma_inj = E_inj/m_p + 1
beta_inj = np.sqrt(1-1/gamma_inj**2)
p_inj = np.sqrt((E_inj+m_p)**2-m_p**2)

gamma_tr = 7.087   # Transition gamma
alpha_0 =  0.01991272251      # 1st order mom. comp. factor stationary
alpha_0_st =  0.01991272251
alpha_1 = -0.06371671224      # 2nd order mom. comp. factor stationary
alpha_2 =  0.9684748052       # 3d order mom. comp. factor stationary

# Tracking details
N_t = 40000           # Number of turns to track
dt_plt = 1000        # Time steps between plots

#RF Parameters
Vamp = 4*25000
U = Vamp * 1/2              # Accelerating Voltage, [eV]

E_fin = E_inj + U * N_t
gamma_fin = E_fin/m_p + 1
p_fin = np.sqrt((E_fin+m_p)**2-m_p**2)
T_s = C_ring/(beta_inj * const.c) * 1e9

print (E_inj, E_fin)
#RF Parameters

# Simulation setup ------------------------------------------------------------
print("Setting up the simulation...")
print("")

alpha_0 = []
n_tr = ((gamma_tr-1)*m_p-E_inj)/U
n_tr = int(n_tr)
n_jump = 6226
d_alpha_0=0.0201662-0.0196684

a1 = - (d_alpha_0) / (3*n_jump)
a2 = d_alpha_0 / n_jump
a3 = - (d_alpha_0) / (3*n_jump)

b1 = alpha_0_st - a1*(n_tr-(2)*n_jump)
b2 = alpha_0_st - a2*(n_tr)
b3 = alpha_0_st - a3*(n_tr+(2)*n_jump)

print(n_tr)

for i in range (0, int(n_tr-2*n_jump)):
    alpha_0.append(alpha_0_st)
    #alpha_1.append(alpha_1_st)
    i = i + 1
print(i)
print(alpha_0[i-1])
for i in range (int(n_tr-2*n_jump), int(n_tr-n_jump/2)):
    alpha_0.append(a1*i+b1)
    i = i + 1
print(i)
print(alpha_0[i-1])
for i in range (int(n_tr-n_jump/2), int(n_tr+n_jump/2)):
    alpha_0.append(a2*i+b2)
    i = i + 1
print(i)
print(alpha_0[i-1])
for i in range (int(n_tr+n_jump/2), int(n_tr+2*n_jump)):
    alpha_0.append(a3*i+b3)
    #alpha_1.append(alpha_1_st)
    i = i + 1
print(i)
print(alpha_0[i-1])
for i in range (int(n_tr+2*n_jump), N_t+1):
    alpha_0.append(alpha_0_st)
    #alpha_1.append(alpha_1_st)
    i = i + 1
print(i)
print(alpha_0[i-1])

# Define general parameters
#np.linspace(E_inj, E_fin, N_t+1)
#ring = Ring(C_ring, alpha, E_inj, Proton(), N_t,
#            synchronous_data_type='kinetic energy')
#ring = Ring(C_ring, alpha_0, np.linspace(E_inj, E_fin, N_t+1), Proton(), N_t,
#            synchronous_data_type='kinetic energy',
#            alpha_1=0, alpha_2=0)
# High order MCF included
ring = Ring(C_ring, alpha_0, np.linspace(E_inj, E_fin, N_t+1), Proton(), N_t,
            synchronous_data_type='kinetic energy',
            alpha_1=alpha_1, alpha_2=alpha_2)

#Jump procedure, transition-energy increasing
#ring = Ring(C_ring, np.linspace(alpha_0, 0.0196684, N_t+1),
#            np.linspace(E_inj, E_fin, N_t+1), Proton(), N_t,
#            synchronous_data_type='kinetic energy',
#            alpha_1=np.linspace(alpha_1, -0.0608641, N_t+1),
#            alpha_2=np.linspace(alpha_2, 0.843772, N_t+1))
#Transition jump
#ring = Ring(C_ring, np.linspace(0.0196684, 0.0201662, N_t+1),
#            np.linspace(E_inj, E_fin, N_t+1), Proton(), N_t,
#            synchronous_data_type='kinetic energy',
#            alpha_1=np.linspace(-0.0608641, -0.0669929, N_t+1),
#            alpha_2=np.linspace(0.843772, 1.09957, N_t+1))


# Define beam and distribution
beam = Beam(ring, N_p, N_b)

h=22
V=100000 #100 kV RF2
harm=22
dphi=np.pi
N_tr = n_tr

Vamp=100000 #100 kV RF2
Volt = []
harmonic = []
dphi_arr = []
for i in range (0, N_tr):
    Volt.append(Vamp)
    harmonic.append(harm)
    dphi_arr.append(dphi)
    i = i + 1
for i in range (0, N_t-N_tr+1):
    Volt.append(Vamp)
    harmonic.append(harm)
    dphi_arr.append(dphi+np.pi/2+np.pi/6)
    i = i + 1
print("RF parameters set")
#print([harmonic], [V], [dphi_arr])
# Define RF station parameters and corresponding tracker
rf = RFStation(ring, [harmonic], [Volt], [dphi_arr])
#rf = RFStation(ring, [h], [V], [dphi])
#rf = RFStation(ring, h, V, phi, N_fourier)
#rf = RFStation(ring, h_acc, V_acc, phi_acc, N_fourier)
long_tracker = RingAndRFTracker(rf, beam, solver='exact')
print("RF set")

# DISTRIBUTIONS
#reinsertion='on'
#seed=1
bigaussian(ring, rf, beam, T_s*1e-9/4/harm/1.5/4, reinsertion=True, seed=1)
#parabolic(ring, rf, beam, tau_0, seed=1)
#generate_coasting_beam(beam, 0, 1.7e-6/harm, spread=3E-3,
#                           spread_type='dp/p', energy_offset=0,
#                           distribution='gaussian', user_distribution=None,
#                           user_probability=None)


# Need slices for the Gaussian fit
profile = Profile(beam, CutOptions(n_slices=100),
                  FitOptions(fit_option='gaussian'))

# Define what to save in file
bunchmonitor = BunchMonitor(ring, rf, beam,
                            this_directory + '../output_files/05_Harm_acc_jump_output_data', Profile=profile)

format_options = {'dirname': this_directory + '../output_files/05_Harm_acc_jump_fig'}
plots = Plot(ring, rf, beam, dt_plt, N_t,
             100e-9/harm*2, 900e-9/harm*2, -100e6, 100e6, xunit='s', separatrix_plot=True,
             Profile=profile, h5file=this_directory + '../output_files/05_Harm_acc_output_data',
             format_options=format_options)

# For testing purposes
test_string = ''
test_string += '{:<17}\t{:<17}\t{:<17}\t{:<17}\n'.format(
    'mean_dE', 'std_dE', 'mean_dt', 'std_dt')
test_string += '{:+10.10e}\t{:+10.10e}\t{:+10.10e}\t{:+10.10e}\n'.format(
    np.mean(beam.dE), np.std(beam.dE), np.mean(beam.dt), np.std(beam.dt))


# Accelerator map
map_ = [long_tracker] + [profile] + [bunchmonitor] + [plots]
print("Map set")
print("")

# Tracking --------------------------------------------------------------------
for i in range(1, N_t + 1):

    # Plot has to be done before tracking (at least for cases with separatrix)
    if (i % dt_plt) == 0:
        print("Outputting at time step %d..." % i)
        print("   Beam momentum %.6e eV" % beam.momentum)
        print("   Beam gamma %3.3f" % beam.gamma)
        print("   Beam beta %3.3f" % beam.beta)
        print("   Beam energy %.6e eV" % beam.energy)
        print("   Four-times r.m.s. bunch length %.4e s" % (4. * beam.sigma_dt))
        print("   Gaussian bunch length %.4e s" % profile.bunchLength)

    # Track
    for m in map_:
        m.track()

    # Define losses according to separatrix and/or longitudinal position
    #beam.losses_separatrix(ring, rf)
    #beam.losses_longitudinal_cut(0., 2.5e-9)

# For testing purposes
test_string += '{:+10.10e}\t{:+10.10e}\t{:+10.10e}\t{:+10.10e}\n'.format(
    np.mean(beam.dE), np.std(beam.dE), np.mean(beam.dt), np.std(beam.dt))
with open(this_directory + '../output_files/05_Harm_acc_jump_text_data.txt', 'w') as f:
    f.write(test_string)

print("Done!")
