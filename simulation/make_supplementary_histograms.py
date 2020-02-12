import matplotlib.pyplot as plt
import numpy as np

from formulae import theory_moments, estimate_general
from params import DEFAULT_PARAMS
from settings import GLOB_PSS_BOUND
from trajectory_plotting import plot_traj_and_mean_sd, plot_means, plot_vars, plot_hist, plot_estimation
from trajectory_simulate import multitraj
from trajectory_analysis import get_moment_timeseries, theory_moments

plt.style.use('parameters.mplstyle')

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

num_traj = 6000
num_steps = 2000
init_bound = 0.0
# model specification
params = DEFAULT_PARAMS
"""
model = 'mode_1'
print(model)
# simulate trajectories
traj_array, times_array = multitraj(num_traj, bound_fraction=init_bound, num_steps=num_steps, model=model, params=params)
# compute moments from data
simdata_1, moment_times_1 = get_moment_timeseries(traj_array, times_array, params, model)
# expectations from theory
theory_curves_1 = theory_moments(moment_times_1, init_bound, method="generating", model=model, p=params)

model = 'mode_2'
# simulate trajectories
traj_array, times_array = multitraj(num_traj, bound_fraction=init_bound, num_steps=num_steps, model=model, params=params)
# compute moments from data
simdata_2, moment_times_2 = get_moment_timeseries(traj_array, times_array, params, model)
# expectations from theory
theory_curves_2 = theory_moments(moment_times_2, init_bound, method="generating", model=model, p=params)
print(model)
"""
model = 'combined'
# simulate trajectories
traj_array, times_array = multitraj(num_traj, bound_fraction=init_bound, num_steps=num_steps, model=model, params=params)
# compute moments from data
simdata_c, moment_times_c = get_moment_timeseries(traj_array, times_array, params, model)
# expectations from theory
theory_curves_c = theory_moments(moment_times_c, init_bound, method="generating", model=model, p=params)
print(model)

model = 'kpr'
# simulate trajectories
traj_array, times_array = multitraj(num_traj, bound_fraction=init_bound, num_steps=num_steps, model=model, params=params)
# compute moments from data
simdata_kpr, moment_times_kpr = get_moment_timeseries(traj_array, times_array, params, model)
# expectations from theory
theory_curves_kpr = theory_moments(moment_times_kpr, init_bound, method="generating", model=model, p=params)
# specify histogram timepoints
print(model)

del traj_array, times_array

#min_time = np.min( [moment_times_1[-1],  moment_times_2[-1],  moment_times_c[-1],  moment_times_kpr[-1]])
min_time = np.min( [moment_times_kpr[-1],  moment_times_c[-1]])

#_, idx1 = find_nearest(moment_times_1, min_time)
#_, idx2 = find_nearest(moment_times_2, min_time)
_, idxc = find_nearest(moment_times_c, min_time)
_, idxkpr = find_nearest(moment_times_kpr, min_time)

print(idxkpr, idxc)
print(min_time)

# model dependent plotting
"""
plot_hist(moment_times_1, simdata_1['distribution_n'], idx1, 'mode_1', state_label='n',
                  theory_mean=theory_curves_1['mean_n'], theory_var=theory_curves_1['var_n'], show=False)

plot_hist(moment_times_2, simdata_2['distribution_m'], idx2, 'mode_2', state_label='m',
                  theory_mean=theory_curves_2['mean_m'], theory_var=theory_curves_2['var_m'], show=False)
"""
plot_hist(moment_times_c, simdata_c['distribution_n'], idxc, 'combined', state_label='n',
                  theory_mean=theory_curves_c['mean_n'], theory_var=theory_curves_c['var_n'], show=False)
plot_hist(moment_times_c, simdata_c['distribution_m'], idxc, 'combined', state_label='m',
                  theory_mean=theory_curves_c['mean_m'], theory_var=theory_curves_c['var_m'], show=False)

plot_hist(moment_times_kpr, simdata_kpr['distribution_n'], idxkpr, 'kpr', state_label='n',
                  theory_mean=theory_curves_kpr['mean_n'], theory_var=theory_curves_kpr['var_n'], show=False)
plot_hist(moment_times_kpr, simdata_kpr['distribution_m'], idxkpr, 'kpr', state_label='m',
                  theory_mean=theory_curves_kpr['mean_m'], theory_var=theory_curves_kpr['var_m'], show=False)

print(min_time)
