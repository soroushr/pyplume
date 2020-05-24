import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append('..')
import pyplume
import matplotlib

plt.rc('font', family='serif')
#plt.rcParams['axes.grid'] = True

E_0 = 0.09 

SAVE_DIR = "./figs/"
if not os.path.exists(SAVE_DIR):
    os.mkdir(SAVE_DIR)

## Model Setup
h_i = 1000. # ice thickness [m]
h_w = 900. # ocean water depth [m]
#discharge = 1000.  # discharge [m^3/s]
#discharges = np.arange(1000, 2000+1, 500)
discharges = np.arange(500, 2500+1, 250)

# Generate lists
depth = np.arange(0, h_w + 1., 1.)
salinity = np.linspace(31.,34.,len(depth))
# linear stratification
fjord_bottom_temps = np.arange(-2, 4+0.1, 0.5)
#fjord_bottom_temps = np.arange(-2,2+0.1,1)

array_discharge = []
for discharge in discharges:

    array_temperature = []
    for fjord_bottom_temp in fjord_bottom_temps:
        temperature = np.linspace(-1, fjord_bottom_temp, len(depth))

        ambient = pyplume.Ambient(h_w, salinity, temperature, depth=depth)

        radius, velocity = pyplume.inlet(h_i, h_w, discharge)

        plume = pyplume.calc_plume(velocity, radius, h_w, ambient,
                                      t_0=-0.5, MELT=True, E_0=E_0)

        if np.any(np.isnan(plume['t_p'])):
            idx = min(np.argwhere(np.isnan(plume['t_p'])))
            plume_top_temperature = plume['t_p'][idx[0]-1]
            plume_peak = plume['z'][idx[0]-1]
        else:
            plume_top_temperature = plume['t_p'][-1]
            plume_peak = np.asarray(plume['z']).max()

        if plume_peak < h_w:
            edgecolor = "none" #'#7f6400'
        else:
            if plume_top_temperature < 0:
                edgecolor = 'black'
            else:
                edgecolor = "#d3d3d3"

        array_temperature.append([plume_top_temperature, edgecolor])
        print('plume height is: {}'.format(plume_peak) )
        print("Q: {}, T: {}".format(discharge, fjord_bottom_temp) )

    array_temperature = np.asarray(array_temperature)

    scatter_plot = plt.scatter(fjord_bottom_temps, discharge*np.ones_like(fjord_bottom_temps),
                                    c=array_temperature[:,0].astype(float), cmap='bwr',
                                    vmin=-2., vmax=2, linewidth=2, s=250, alpha=0.8,
                                    marker='o', edgecolor=array_temperature[:,1])

plt.colorbar(scatter_plot, label='temperature of plume top ($^{\circ}$C)')
plt.title("supercooling as a function of $Q$ and $T_{a0}$")
plt.xlabel("$T_{a0}$ (surface temperature fixed at -2$^{\circ}$C)")
plt.ylabel("$Q$, subglacial discharge (m$^3$/s)")
plt.savefig(SAVE_DIR + 'scatter-test.png', dpi=300)
print("temperature saved")
plt.clf()
