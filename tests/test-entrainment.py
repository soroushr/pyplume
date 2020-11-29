import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append('..')
import pyplume
import gsw

plt.rc('font', family='serif')

SAVE_DIR = "./figs/entrainment/"
if not os.path.exists(SAVE_DIR):
    os.mkdir(SAVE_DIR)

def get_pressure(given_depth):
    return given_depth*(1027.*9.81*1.e-4)

## Model Setup
h_i = 400. # ice thickness [m]
h_w = 350. # ocean water depth [m]
discharge = 10000.  # discharge [m^3/s]

# Generate lists
depth = np.arange(0, h_w + 1., 1.)
#salinity = np.linspace(30.,33.,len(depth))
#temperature = np.linspace(-2.,4.,len(depth))

## 75% 4 degrees, linearly to -2
t_top, t_mid, t_bot = 0, 1, 2

temperature_top = np.linspace(t_top, t_mid, int(len(depth)/4)+1)
temperature_bot = np.linspace(t_mid, t_bot, int(3*len(depth)/4))
temperature = np.hstack([temperature_top, temperature_bot])

s_top, s_mid, s_bot = 30, 32.0, 34.0

salinity_top = np.linspace(s_top, s_mid, int(len(depth)/4)+1)
salinity_bot = np.linspace(s_mid, s_bot, int(3*len(depth)/4))
salinity = np.hstack([salinity_top, salinity_bot])

# modify top layer's temperature and salinity
# to match Mary Louise Timmermanns' suggestion
# "prime" the surface for 10 meters
priming_depth = 50

intrusion_depths = [0]

for int_depth in intrusion_depths:
    salinity[int_depth:int_depth+priming_depth] = 0.
    temperature[int_depth:int_depth+priming_depth] = gsw.t_freezing(p = get_pressure(int_depth),
                                                                    SA = 10,
                                                                    saturation_fraction = 0)

# Make the profiles into an Ambient object
ambient = pyplume.Ambient(h_w, salinity, temperature, depth=depth)

# get ambient freezing profile
ambient_pmp = gsw.t_freezing(p = ambient.pressure,
                               SA = ambient.salinity,
                               saturation_fraction = 0)

radius, velocity = pyplume.inlet(h_i, h_w, discharge)

# a range for entrainment
#E_0s = np.arange(0.01, 0.09+0.01, 0.01)
E_0s = [0.01, 0.05, 0.09]

for E_0 in E_0s[::-1]:
    plume = pyplume.calc_plume(velocity, radius, h_w, ambient,
                               t_0=-0.5, MELT=True, E_0=E_0)

    pressure = [ get_pressure(ff) for ff in plume['z'][::-1] ]
    plume_pmp = gsw.t_freezing(p = pressure,
                                   SA = plume['s_p'],
                                   saturation_fraction = 0)

    f,ax = plt.subplots(1,5, sharey=True, figsize=(15,8))

    if np.any(np.isnan(plume['t_p'])):
        idx = min(np.argwhere(np.isnan(plume['t_p'])))
        plume_top = plume['t_p'][idx[0]-1]
        plume_pmp_top = plume_pmp[idx[0]-1]
    else:
        idx = len(plume['t_p'])
        plume_top = plume['t_p'][-1]


    print("temperature PMP of plume top is: {}".format(plume_pmp_top) )
    print("temperature of plume top is: {}".format(plume_top) )

    if plume_top > plume_pmp_top:
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
    ax[1].set_xlim(0, 10)
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

    # pressure melting temperature profile
    ax[4].plot(plume['t_p']-plume_pmp, plume['z'], color, label='Plume Profile')
    ax[4].plot(ambient.temperature-ambient_pmp, ambient.z, 'k:', label= 'Ambient Profile')
    ax[4].set_xlim(-1, 5)
    ax[4].set_xlabel(r'T-T$_{pmp}$ ($^{\circ}$C)')

    # Add in the neutral buoyancy height
    nb_height = max(plume['z'])#[~np.isnan(plume['b_p'])])
    for axis in ax:
        axis.plot([-1000, 1000], [nb_height, nb_height], '--', label = 'Neutral Buoyancy')
        axis.grid(linestyle='dotted')

    ax[-1].legend(loc='best')

    plt.gcf().suptitle("entrainment={:.2f}, T_p={:.4f}".format(E_0, plume_top))
    plt.savefig(SAVE_DIR + 'E_0_{:.2f}.png'.format(E_0), dpi=300)
    print("entrainment {:.2f} saved - {}".format(E_0, color))
    plt.clf()
