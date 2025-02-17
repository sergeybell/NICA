import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

IBS_time = pd.read_table('IBS_data.dat', sep="\s+")
structure = pd.read_table('IBS_structure.dat', sep="\s+")
beam = pd.read_table('IBS_beam.dat', sep="\s+")

# do same but attach it to the dataframe
IBS_time['t_all'] = IBS_time.apply(lambda row: 1/((1/np.abs(row.t_h)) + (1/np.abs(row.t_v))), axis=1)
IBS_time['t_l_abs'] = IBS_time.apply(lambda row: np.abs(row.t_l), axis=1)
IBS_time['t'] = IBS_time[['t_l_abs','t_h', 't_v']].min(axis=1)
print(IBS_time)

#print(IBS_time)
#temp_X = parameters_X.merge(particles_X, on=['I'])
#data_X = temp_X.merge(tss_X, on=['I'])
#IBS_time[IBS_time['structure'] == 'Regular'].plot('kinetic_energy', y = ['t_l', 't_h', 't_v'])
#IBS_time[IBS_time['structure'] == 'Resonant'].plot('kinetic_energy', y = ['t_l', 't_h', 't_v'])
#IBS_time[IBS_time['structure'] == 'Regular_and_Resonant'].plot('kinetic_energy', y = ['t_l', 't_h', 't_v'])

regular = IBS_time[IBS_time['structure'] == 'Regular']
resonant = IBS_time[IBS_time['structure'] == 'Resonant']
regular_and_resonant = IBS_time[IBS_time['structure'] == 'Regular_and_Resonant']
print(regular)

plt.plot(regular['kinetic_energy'], regular['t_all'], 'k')
plt.plot(resonant['kinetic_energy'], resonant['t_all'], 'b')
plt.plot(regular_and_resonant['kinetic_energy'], regular_and_resonant['t_all'], 'r')
plt.legend(['Regular', 'Resonant', 'Regular and resonant'])
plt.title("Intrabeam Scattering")
plt.xlabel('Energy, GeV/u')
plt.ylabel('IBS time, sec')
plt.grid()
plt.xlim(1, 4.5)
plt.ylim(0, 3000)
plt.show()

plt.plot(regular['kinetic_energy'], regular['t'], 'k')
plt.plot(resonant['kinetic_energy'], resonant['t'], 'b')
plt.plot(regular_and_resonant['kinetic_energy'], regular_and_resonant['t'], 'r')
plt.legend(['Regular', 'Resonant', 'Regular and resonant'])
plt.title("Intrabeam Scattering")
plt.xlabel('Energy, GeV/u')
plt.ylabel('IBS time, sec')
plt.grid()
plt.xlim(1, 4.5)
plt.ylim(0, 3000)
plt.show()

plt.plot(regular['kinetic_energy'], regular['t_h'], 'k')
plt.plot(resonant['kinetic_energy'], resonant['t_h'], 'b')
plt.plot(regular_and_resonant['kinetic_energy'], regular_and_resonant['t_h'], 'r')
plt.legend(['Regular', 'Resonant', 'Regular and resonant'])
plt.title("IBS time in horisontal plane")
plt.show()

plt.plot(regular['kinetic_energy'], regular['t_v'], 'k')
plt.plot(resonant['kinetic_energy'], resonant['t_v'], 'b')
plt.plot(regular_and_resonant['kinetic_energy'], regular_and_resonant['t_v'], 'r')
plt.legend(['Regular', 'Resonant', 'Regular and resonant'])
plt.title("IBS time in vertical plane")
plt.show()

plt.plot(regular['kinetic_energy'], regular['t_l'], 'k')
plt.plot(resonant['kinetic_energy'], resonant['t_l'], 'b')
plt.plot(regular_and_resonant['kinetic_energy'], regular_and_resonant['t_l'], 'r')
plt.legend(['Regular', 'Resonant', 'Regular and resonant'])
plt.title("IBS time in longitudinal plane")
plt.show()
