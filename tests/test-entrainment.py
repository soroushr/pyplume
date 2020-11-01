import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append('..')
import pyplume

plt.rc('font', family='serif')

SAVE_DIR = "./figs/"
if not os.path.exists(SAVE_DIR):
    os.mkdir(SAVE_DIR)

## Model Setup
h_i = 500. # ice thickness [m]
h_w = 450. # ocean water depth [m]
discharge = 500.  # discharge [m^3/s]

# Generate lists
depth = np.arange(0, h_w + 1., 1.)
salinity = np.linspace(30.,33.,len(depth))
#temperature = np.linspace(-2.,4.,len(depth))

## 75% 4 degrees, linearly to -2
t_top, t_mid, t_bot = -2, 0, 0.5

temperature_top = np.linspace(t_top, t_mid, int(len(depth)/4)+1)
temperature_bot = np.linspace(t_mid, t_bot, int(3*len(depth)/4))
temperature = np.hstack([temperature_top, temperature_bot])

# Make the profiles into an Ambient object
ambient = pyplume.Ambient(h_w, salinity, temperature, depth=depth)

radius, velocity = pyplume.inlet(h_i, h_w, discharge)

# a range for entrainment
#E_0s = np.arange(0.01, 0.09+0.01, 0.01)
E_0s = [0.09]

for E_0 in E_0s[::-1]:
    plume = pyplume.calc_plume(velocity, radius, h_w, ambient,
                               t_0=-0.5, MELT=True, E_0=E_0)

    f,ax = plt.subplots(1,4, sharey=True, figsize=(15,8))

    if np.any(np.isnan(plume['t_p'])):
        idx = min(np.argwhere(np.isnan(plume['t_p'])))
        plume_top = plume['t_p'][idx[0]-1]
    else:
        idx = len(plume['t_p'])
        plume_top = plume['t_p'][-1]

    if plume_top > 0:
        color='red'
    else:
        color='blue'

    # radius
    ax[0].plot(plume['b_p'], plume['z'], color)
    ax[0].set_xlim(0, 120)
    ax[0].set_ylim(0, h_w)

    ax[0].set_xlabel('Radius (m)')
    ax[0].set_ylabel('Height above source (m)')

    # velocity
    ax[1].plot(plume['w_p'], plume['z'], color)
    ax[1].set_xlim(0, 9)
    ax[1].set_xlabel('Velocity (m/s)')

    # temperature
    ax[2].plot(plume['t_p'], plume['z'], color)
    ax[2].plot(ambient.temperature, ambient.z, 'k:')
    ax[2].set_xlim(-2, 4)
    ax[2].set_xlabel('Temperature ($^{\circ}$C)')

    # salinity
    ax[3].plot(plume['s_p'], plume['z'], color, label='Plume Profile')
    ax[3].plot(ambient.salinity, ambient.z, 'k:', label= 'Ambient Profile')
    ax[3].set_xlim(0, 35)
    ax[3].set_xlabel('Salinity (g/kg)')

    # Add in the neutral buoyancy height
    nb_height = max(plume['z'])#[~np.isnan(plume['b_p'])])
    for axis in ax:
        axis.plot([-1000, 1000], [nb_height, nb_height], '--', label = 'Neutral Buoyancy')
        axis.grid(linestyle='dotted')
        
    ax[3].legend(loc=2)

    plt.gcf().suptitle("entrainment={:.2f}, T_p={:.4f}".format(E_0, plume_top))
    plt.savefig(SAVE_DIR + 'E_0_{:.2f}.png'.format(E_0), dpi=300)
    print("entrainment {:.2f} saved - {}".format(E_0, color))
    plt.clf()
