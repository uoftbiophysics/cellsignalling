import matplotlib.pyplot as plt
import numpy as np

from formulae import theory_cumulants
from trajectory_simulate import multitraj


def get_state_at_t(traj, times, t):
    step = 0
    while times[step] <= t:
        step += 1
    return traj[step - 1, :]


def get_mean_var_timeseries(traj_array, times_array):
    num_traj = np.shape(traj_array)[-1]
    dt = np.mean(times_array[1, :])
    endtime = np.min(times_array[-1, :])
    moment_times = np.arange(0.0, endtime, dt)
    mean_vals = np.zeros(len(moment_times))
    var_vals = np.zeros(len(moment_times))
    for idx, t in enumerate(moment_times):
        statesum = 0.0
        statesquaresum = 0.0

        for k in xrange(num_traj):
            state_at_t = get_state_at_t(traj_array[:, :, k], times_array[:, k], t)
            statesum += state_at_t[1]
            statesquaresum += state_at_t[1]**2

        mean_vals[idx] = statesum / num_traj
        var_vals[idx] = statesquaresum / num_traj - mean_vals[idx]**2
    return mean_vals, var_vals, moment_times


if __name__ == '__main__':
    # settings
    num_traj = 500
    num_steps = 200
    init_bound = 1.0

    # compute
    traj_array, times_array = multitraj(num_traj, bound_fraction=init_bound, num_steps=num_steps)
    mean_vals, var_vals, moment_times = get_mean_var_timeseries(traj_array, times_array)
    mean_vals_direct, var_vals_direct = theory_cumulants(moment_times, "direct", init_bound=init_bound)
    mean_vals_gen, var_vals_gen = theory_cumulants(moment_times, "generating", init_bound=init_bound)

    # plot trajectories
    for k in xrange(num_traj):
        times_k = times_array[:, k]
        traj_k = traj_array[:, 1, k]
        plt.plot(times_k, traj_k, '--', lw=0.5, alpha=0.5)
    # plot trajectory cumulants
    plt.plot(moment_times, mean_vals, 'k', lw=2, label="data")
    plt.plot(moment_times, mean_vals + np.sqrt(var_vals), '--k', lw=2)
    plt.plot(moment_times, mean_vals - np.sqrt(var_vals), '--k', lw=2)
    # plot theory cumulants (direct)
    """
    plt.plot(moment_times, mean_vals_direct, 'b', lw=2, label="direct")
    plt.plot(moment_times, mean_vals_direct + np.sqrt(var_vals_direct), '--b', lw=2)
    plt.plot(moment_times, mean_vals_direct - np.sqrt(var_vals_direct), '--b', lw=2)
    """
    # plot theory cumulants (generating function)
    plt.plot(moment_times, mean_vals_gen, 'r', lw=2, label="generating")
    plt.plot(moment_times, mean_vals_gen + np.sqrt(var_vals_gen), '--r', lw=2)
    plt.plot(moment_times, mean_vals_gen - np.sqrt(var_vals_gen), '--r', lw=2)
    # decorate
    plt.title('Mode 1 <n>(t) +- sqrt(var(t)) for %d trajectories' % num_traj)
    plt.xlabel('time')
    plt.ylabel('n')
    plt.legend()
    plt.show()

    # only means
    plt.figure()
    plt.scatter(moment_times, mean_vals, s=4.0, c='k', marker='s', label="data", alpha=0.5)
    plt.scatter(moment_times, mean_vals_direct, s=4.0, c='b', marker='s', label="direct", alpha=0.5)
    plt.scatter(moment_times, mean_vals_gen, s=4.0, c='r', marker='s', label="generating", alpha=0.5)
    plt.title('Mode 1 mean(n) for %d trajectories' % num_traj)
    plt.xlabel('time')
    plt.ylabel('<n>(t)')
    plt.legend()
    plt.show()

    # only vars
    plt.figure()
    plt.plot(moment_times, var_vals, '--k', lw=2, label="data")
    plt.plot(moment_times, var_vals_direct, '--b', lw=2, label="direct")
    plt.plot(moment_times, var_vals_gen, '--r', lw=2, label="generating")
    plt.title('Mode 1 var(n) for %d trajectories' % num_traj)
    plt.xlabel('time')
    plt.ylabel('var(n)')
    plt.legend()
    plt.show()
