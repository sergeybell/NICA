# Copyright 2014-2017 CERN. This software is distributed under the
# terms of the GNU General Public Licence version 3 (GPL Version 3),
# copied verbatim in the file LICENCE.md.
# In applying this licence, CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.
# Project website: http://blond.web.cern.ch/

'''
Stitching 2 stages of transition crossing:
– Acceleration (Accumulation) near transition
    at constant distance toward transition energy;
– Rapid Jump wo RF Barriers, beam accelerates.

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
import pylab as plt

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

# Impedance
from blond.impedances.impedance import (InducedVoltageFreq, InducedVoltageTime,
                                        InductiveImpedance,
                                        TotalInducedVoltage,
                                        InducedVoltageResonator)
from blond.impedances.impedance_sources import Resonators
from blond.plots.plot_impedance import plot_induced_voltage_vs_bin_centers

mpl.use('Agg')
this_directory = os.path.dirname(os.path.realpath(__file__)) + '/'
os.makedirs(this_directory + '../output_files/Test_07_fig', exist_ok=True)

# Simulation parameters -------------------------------------------------------
# Bunch parameters
N_b = 2.0e13                 # Design proton Intensity in NICA proton 13 GeV
N_p = 50000                  # Macro-particles

# Machine
m_p = 0.938e9          # mass of proton in eV
C_ring = 503.04        # Machine circumference [m]
#E_inj = 5.66704e9          # Injection Energy, [eV]
E_inj = 5.70832e9          # Injection Energy, [eV] before JUMP
gamma_inj = E_inj/m_p + 1
beta_inj = np.sqrt(1-1/gamma_inj**2)
p_inj = np.sqrt((E_inj+m_p)**2-m_p**2)
U_2 = 0             # Accelerating Voltage, [eV] EQUAL to Accumulation
U_3 = 300             # Accelerating Voltage, [eV] at JUMP
#gamma_tr = 7.087     # Transition gamma
#alpha_0 =  0.01991272251      # 1st order mom. comp. factor stationary
#alpha_1 = -0.06371671224      # 2nd order mom. comp. factor stationary
#alpha_2 =  0.9684748052       # 3d order mom. comp. factor stationary

alpha_0 =  0.0196684      # 1st order mom. comp. factor stationary
alpha_1 = -0.0608641      # 2nd order mom. comp. factor stationary
alpha_2 =  0.843772       # 3d order mom. comp. factor stationary


# Tracking details
N_t_2 = 5           # Number of turns to track BEFORE JUMP
N_t_3 = 5           # Number of turns to track AT JUMP
dt_plt_2 = 1        # Time steps between plots
dt_plt_3 = 1        # Time steps between plots

E_fin = E_inj + U * N_t
gamma_fin = E_fin/m_p + 1
p_fin = np.sqrt((E_fin+m_p)**2-m_p**2)

print (E_inj, E_fin)
#RF Parameters
N_fourier = 25 # Amount of Fourirer
m = 3
harm = 2
#T_s = 1701
T_s = C_ring/(beta_inj * const.c) * 1e9
print(T_s)
# BARRIER BUCKET
V_peak = 5000
T_1 = 50 # ns
T_2 = (T_s - 2*T_1)/(2*harm) # ns
phi_r = (2*T_1) / T_2
h_r = np.pi/phi_r
phi_rf = np.pi

# ACCELERATING
T_acc = T_s
phi_acc = 0

h_2, V_2, phi_2 = Barrier_Bucket(V_peak, T_1, harm, C_ring, gamma_inj, gamma_tr, phi_rf, N_fourier)

# Simulation setup ------------------------------------------------------------
print("Setting up the simulation...")
print("")

# Define general parameters
ring_2 = Ring(C_ring, alpha_0, E_inj, Proton(), N_t,
            synchronous_data_type='kinetic energy',
            alpha_1=alpha_1, alpha_2=alpha_2)
ring_3 = Ring(C_ring, np.linspace(0.0196684, 0.0201662, N_t+1),
            np.linspace(E_inj, E_fin, N_t+1), Proton(), N_t,
            synchronous_data_type='kinetic energy',
            alpha_1=np.linspace(-0.0608641, -0.0669929, N_t+1),
            alpha_2=np.linspace(0.843772, 1.09957, N_t+1))

# Define beam and distribution
beam = Beam(ring, N_p, N_b)

# Define RF station parameters and corresponding tracker
h_3 = 1
V_3 = 10
phi_rf_3 = 0

rf_2 = RFStation(ring_2, h_2, V_2, phi_2, N_fourier)
rf_3 = RFStation(ring_3, [h_3], [V_3], [phi_3])

long_tracker = RingAndRFTracker(rf, beam, solver='exact')
print (h_2, V_2, phi_2, N_fourier)
print (h_3, V_3, phi_3)

# DISTRIBUTIONS
#bigaussian(ring, rf, beam, T_s*1e-9/4/harm/1.7, 1e7, reinsertion='on', seed=1)
generate_coasting_beam(beam, 0, 1.7e-6/harm, spread=3E-3,
                           spread_type='dp/p', energy_offset=0,
                           distribution='gaussian', user_distribution=None,
                           user_probability=None)

# Need slices for the Gaussian fit
profile = Profile(beam, CutOptions(n_slices=100),
                  FitOptions(fit_option='gaussian'))

# Define what to save in file
bunchmonitor = BunchMonitor(ring, rf, beam,
                            this_directory + '../output_files/Test_07_output_data', Profile=profile)

format_options = {'dirname': this_directory + '../output_files/Test_07_fig'}
plots = Plot(ring, rf, beam, dt_plt, N_t,
             -250e-9/harm*4, 1000e-9/harm*4, -100e6, 100e6, xunit='s', separatrix_plot=True,
             Profile=profile, h5file=this_directory + '../output_files/Test_07_output_data',
             format_options=format_options)

# For testing purposes
test_string = ''
test_string += '{:<17}\t{:<17}\t{:<17}\t{:<17}\n'.format(
    'mean_dE', 'std_dE', 'mean_dt', 'std_dt')
test_string += '{:+10.10e}\t{:+10.10e}\t{:+10.10e}\t{:+10.10e}\n'.format(
    np.mean(beam.dE), np.std(beam.dE), np.mean(beam.dt), np.std(beam.dt))

# Accelerator map
map_ = [profile] + [total_ind_volt_ZoN] + [total_ind_volt_time] +\
        [long_tracker] + [bunchmonitor] +[plots]
#map_ = [profile] + [total_ind_volt_ZoN] +\
#        [long_tracker] + [bunchmonitor] +[plots]
print("Map set")
print("")

# Tracking --------------------------------------------------------------------
for i in range(N_t):

    # Plot has to be done before tracking (at least for cases with separatrix)
    if (i % dt_plt) == 0:
        print("Outputting at time step %d..." % i)
        print("   Beam momentum %.6e eV" % beam.momentum)
        print("   Beam gamma %3.3f" % beam.gamma)
        print("   Beam beta %3.3f" % beam.beta)
        print("   Beam energy %.6e eV" % beam.energy)
        print("   Four-times r.m.s. bunch length %.4e s" % (4. * beam.sigma_dt))
        print("   Gaussian bunch length %.4e s" % profile.bunchLength)

    if (i % dt_plt) == 0:
        plot_induced_voltage_vs_bin_centers(total_ind_volt_ZoN, figure_index=i,
                                            style='.', dirname=this_directory + '../output_files/Test_07_fig/1')
        plot_induced_voltage_vs_bin_centers(total_ind_volt_time, figure_index=i,
                                            style='.', dirname=this_directory + '../output_files/Test_07_fig/2')
    # Track
    for m in map_:
        m.track()
    # Define losses according to separatrix and/or longitudinal position
    #beam.losses_separatrix(ring, rf)
    #beam.losses_longitudinal_cut(0., 2.5e-9)

# For testing purposes
test_string += '{:+10.10e}\t{:+10.10e}\t{:+10.10e}\t{:+10.10e}\n'.format(
    np.mean(beam.dE), np.std(beam.dE), np.mean(beam.dt), np.std(beam.dt))
with open(this_directory + '../output_files/Test_07_test_data.txt', 'w') as f:
    f.write(test_string)

plt.figure()
plt.plot(profile.bin_centers * 1e9, total_ind_volt_time.induced_voltage,
         lw=2, alpha=0.75, label='Resonator time domain')
plt.plot(profile.bin_centers * 1e9, total_ind_volt_ZoN.induced_voltage,
         lw=2, alpha=0.75, label=r'Z/n = -10 $\Omega$')
plt.xlabel('Time [ns]')
plt.ylabel('Induced voltage [V]')
plt.legend(loc=2, fontsize='medium')
plt.savefig(this_directory + '../output_files/Test_07_fig/fig.png')

print("Done!")
