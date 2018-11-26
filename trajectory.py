import numpy as np
import matplotlib.pyplot as plt

"""
State for mode 1, 2: [bool, n]
    - bool == 1 implies BOUND, 0 implies UNBOUND
    - n is number of downstream molecules produced (due to being bound)
"""

# globals
num_rxn_mode1 = 3
state_size_mode1 = 2

# model parameters
c = 4.0
k_on = 5.0
k_off = 1.0
k_p = 10.0
x = k_on*c/k_off

# initial conditions
n0 = 0.0
m0 = 0.0
p0 = 0.0
INIT_COND_MODE1 = [n0, p0]

# misc
NUM_STEPS = 100

# rxn update dictionaries for each mode
update_dict_mode1 = {0: np.array([1.0, 0.0]),
                     1: np.array([-1.0, 0.0]),
                     2: np.array([0.0, 1.0])}


def propensities_mode1(state):
    propensities = np.zeros(num_rxn_mode1)
    propensities[0] = k_on * c * (1 - state[0])  # bind (0.0 if already bound)
    propensities[1] = k_off * state[0]           # unbind (0.0 if already unbound)
    propensities[2] = k_p * state[0]             # produce one n molecule
    return propensities


def update_state(state_timeseries, rxn_idx, step):
    current_state = state_timeseries[step]
    increment = update_dict_mode1[rxn_idx]
    state_timeseries[step + 1, :] = current_state + increment
    return state_timeseries


def simulate_traj(num_steps=NUM_STEPS, init_cond=INIT_COND_MODE1):
    times = np.zeros(num_steps)
    traj = np.zeros((num_steps, state_size_mode1))
    traj[0, :] = init_cond
    times[0] = 0.0
    for step in xrange(num_steps-1):
        # generate two U[0,1] random variables
        r1, r2 = np.random.random(2)

        # generate propensities and their partitions
        alpha = propensities_mode1(traj[step])
        alpha_partitions = np.zeros(len(alpha)+1)
        alpha_sum = 0.0
        for i in xrange(len(alpha)):
            alpha_sum += alpha[i]
            alpha_partitions[i + 1] = alpha_sum

        # find time to first reaction
        tau = np.log(1 / r1) / alpha_sum

        # pick a reaction
        r2_scaled = alpha_sum * r2
        for rxn_idx in xrange(len(alpha)):
            if alpha_partitions[rxn_idx] <= r2_scaled < alpha_partitions[rxn_idx + 1]:  # i.e. rxn_idx has occurred
                break

        # update state
        traj = update_state(traj, rxn_idx, step)
        times[step+1] = times[step] + tau
    return traj, times


def multitraj(num_traj, num_steps=NUM_STEPS):
    traj_array = np.zeros((num_steps, state_size_mode1, num_traj))
    times_array = np.zeros((num_steps, num_traj))

    for k in xrange(num_traj):
        traj, times = simulate_traj(num_steps=num_steps)
        traj_array[:, :, k] = traj
        times_array[:, k] = times
    return traj_array, times_array


def get_state_at_t(traj, times, t):
    step = 0
    while times[step] < t:
        step += 1
    return traj[step, :]


def get_mean_timeseries(traj_array, times_array):
    num_traj = np.shape(traj_array)[-1]
    dt = np.mean(times_array[1,:])
    endtime = np.min(times_array[-1,:])
    mean_times = np.arange(0.0,endtime, dt)
    mean_vals = np.zeros(len(mean_times))
    for idx, t in enumerate(mean_times):
        statesum = 0.0
        for k in xrange(num_traj):
            state_at_t = get_state_at_t(traj_array[:,:,k], times_array[:,k], t)
            statesum += state_at_t[1]
        mean_vals[idx] = statesum / k
    return mean_vals, mean_times


if __name__ == '__main__':
    num_traj = 100
    traj_array, times_array = multitraj(num_traj)
    for k in xrange(num_traj):
        times_k = times_array[:, k]
        traj_k = traj_array[:, 1, k]
        plt.plot(times_k, traj_k, '--', lw=1)
    mean_vals, mean_times = get_mean_timeseries(traj_array, times_array)
    plt.plot(mean_times, mean_vals, 'k', lw=2)
    plt.title('Mode 1 <n>(t) for %d trajectories' % num_traj)
    plt.xlabel('time')
    plt.ylabel('n')
    plt.show()
