import matplotlib.pyplot as plt
import numpy as np

from formulae import theory_moments
from trajectory_simulate import multitraj


def get_state_at_t(traj, times, t):
    step = 0
    while times[step] <= t:
        step += 1
    return traj[step - 1, :]


def get_moment_timeseries(traj_array, times_array):
    """
    Returns, as a dict, multiple output timeseries (aligned to moment_times) depending on the model:
        - mode_1: mean(n)(t), var(n)(t), None, None, None
        - mode_2: None, None, mean(m)(t), var(m)(t), None
        - combined: mean(n)(t), var(n)(t), mean(m)(t), var(m)(t), cov(n,m)(t)
    """
    # prepare moment times
    num_traj = np.shape(traj_array)[-1]
    dt = np.mean(times_array[1, :])
    endtime = np.min(times_array[-1, :])
    moment_times = np.arange(0.0, endtime, dt)
    # prepare value dict
    moment_curves = {'mean_n': None,
                     'var_n': None,
                     'mean_m': None,
                     'var_m': None,
                     'cov_nm': None}

    if model in ['mode_1', 'mode_2']:
        mean_1or2 = {'mode_1': 'mean_n', 'mode_2': 'mean_m'}[model]
        var_1or2 = {'mode_1': 'var_n', 'mode_2': 'var_m'}[model]
        moment_curves[mean_1or2] = np.zeros(len(moment_times))
        moment_curves[var_1or2] = np.zeros(len(moment_times))
        for idx, t in enumerate(moment_times):
            statesum = 0.0
            statesquaresum = 0.0
            for k in xrange(num_traj):
                state_at_t = get_state_at_t(traj_array[:, :, k], times_array[:, k], t)
                statesum += state_at_t[1]
                statesquaresum += state_at_t[1]**2
            moment_curves[mean_1or2][idx] = statesum / num_traj
            moment_curves[var_1or2][idx] = statesquaresum / num_traj - moment_curves[mean_1or2][idx]**2

    else:
        assert model == 'combined'
        moment_curves = {key: np.zeros(len(moment_times)) for key in moment_curves.keys()}
        for idx, t in enumerate(moment_times):
            statesum_n = 0.0
            statesquaresum_n = 0.0
            statesum_m = 0.0
            statesquaresum_m = 0.0
            stateprod_nm = 0.0
            for k in xrange(num_traj):
                state_at_t = get_state_at_t(traj_array[:, :, k], times_array[:, k], t)
                statesum_n += state_at_t[1]
                statesquaresum_n += state_at_t[1] ** 2
                statesum_m += state_at_t[2]
                statesquaresum_m += state_at_t[2] ** 2
                stateprod_nm += state_at_t[1] * state_at_t[2]

            moment_curves['mean_n'][idx] = statesum_n / num_traj
            moment_curves['mean_m'][idx] = statesum_m / num_traj
            moment_curves['var_m'][idx] = statesquaresum_n / num_traj - moment_curves['mean_n'][idx]**2
            moment_curves['var_m'][idx] = statesquaresum_m / num_traj - moment_curves['mean_m'][idx]**2
            moment_curves['cov_nm'][idx] = stateprod_nm / num_traj - \
                                           moment_curves['mean_n'][idx] * moment_curves['mean_m'][idx]
    return moment_curves, moment_times


if __name__ == '__main__':
    # settings
    model = 'mode_1'
    num_traj = 200
    num_steps = 200
    init_bound = 1.0

    # compute
    traj_array, times_array = multitraj(num_traj, bound_fraction=init_bound, num_steps=num_steps, model=model)
    data_moments, moment_times = get_moment_timeseries(traj_array, times_array)
    # theory
    theory_curves_gen = theory_moments(moment_times, init_bound, method="generating", model=model)
    mean_vals_gen = theory_curves_gen['mean_n']
    var_vals_gen = theory_curves_gen['var_n']

    # plot trajectories
    for k in xrange(num_traj):
        times_k = times_array[:, k]
        traj_k = traj_array[:, 1, k]
        plt.plot(times_k, traj_k, '--', lw=0.5, alpha=0.5)
    # plot trajectory cumulants
    plt.plot(moment_times, data_moments['mean_n'], 'k', lw=2, label="data")
    plt.plot(moment_times, data_moments['mean_n'] + np.sqrt(data_moments['var_n']), '--k', lw=2)
    plt.plot(moment_times, data_moments['mean_n'] - np.sqrt(data_moments['var_n']), '--k', lw=2)
    # plot theory cumulants (generating function)
    plt.plot(moment_times, mean_vals_gen, 'r', lw=2, label="generating")
    plt.plot(moment_times, mean_vals_gen + np.sqrt(var_vals_gen), '--r', lw=2)
    plt.plot(moment_times, mean_vals_gen - np.sqrt(var_vals_gen), '--r', lw=2)
    # decorate
    plt.title('%s <n>(t) +- sqrt(var(t)) for %d trajectories' % (model, num_traj))
    plt.xlabel('time')
    plt.ylabel('n')
    plt.legend()
    plt.show()

    # only means
    plt.figure()
    plt.scatter(moment_times, data_moments['mean_n'], s=4.0, c='k', marker='s', label="data", alpha=0.5)
    plt.scatter(moment_times, mean_vals_gen, s=4.0, c='r', marker='s', label="generating", alpha=0.5)
    plt.title('%s mean(n) for %d trajectories' % (model, num_traj))
    plt.xlabel('time')
    plt.ylabel('<n>(t)')
    plt.legend()
    plt.show()

    # only vars
    plt.figure()
    plt.plot(moment_times, data_moments['var_n'], '--k', lw=2, label="data")
    plt.plot(moment_times, var_vals_gen, '--r', lw=2, label="generating")
    plt.title('%s var(n) for %d trajectories' % (model, num_traj))
    plt.xlabel('time')
    plt.ylabel('var(n)')
    plt.legend()
    plt.show()
