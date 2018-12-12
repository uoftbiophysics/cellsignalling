import matplotlib.pyplot as plt
import numpy as np
import os

from settings import FOLDER_OUTPUT


def plot_traj_and_mean_sd(traj_array, times_array, moment_times, model, data_mean=None, data_var=None, theory_mean=None,
                          theory_var=None, state_label='n', state_idx=1, title='', show=True):
    # plot individual trajectories
    num_traj = traj_array.shape[-1]
    for k in xrange(num_traj):
        times_k = times_array[:, k]
        traj_k = traj_array[:, state_idx, k]
        plt.plot(times_k, traj_k, '--', lw=0.5, alpha=0.5)
    # plot trajectory moments
    if data_mean is not None:
        plt.plot(moment_times, data_mean, 'k', lw=2, label="data")
        if data_var is not None:
            plt.plot(moment_times, data_mean + np.sqrt(data_var), '--k', lw=2)
            plt.plot(moment_times, data_mean - np.sqrt(data_var), '--k', lw=2)
    # plot theory moments (generating function)
    if theory_mean is not None:
        plt.plot(moment_times, theory_mean, 'r', lw=2, label="generating")
        if theory_var is not None:
            plt.plot(moment_times, theory_mean + np.sqrt(theory_var), '--r', lw=2)
            plt.plot(moment_times, theory_mean - np.sqrt(theory_var), '--r', lw=2)
    # decorate
    plt.title(title)
    plt.xlabel('time')
    plt.ylabel(state_label)
    plt.legend()
    plt.savefig(FOLDER_OUTPUT + os.sep + 'traj_%s_%s.png' % (model, state_label))
    if show:
        plt.show()
    return plt.gca()


def plot_means(moment_times, data_mean, model, theory_mean=None, title='', state_label='n', show=True):
    plt.scatter(moment_times, data_mean, s=4.0, c='k', marker='s', label="data", alpha=0.5)
    plt.scatter(moment_times, theory_mean, s=4.0, c='r', marker='s', label="generating", alpha=0.5)
    plt.title(title)
    plt.xlabel('time')
    plt.ylabel('<%s>(t)' % state_label)
    plt.legend()
    plt.savefig(FOLDER_OUTPUT + os.sep + 'traj_%s_mean_%s.png' % (model, state_label))
    if show:
        plt.show()
    return plt.gca()


def plot_vars(moment_times, data_var, model, theory_var=None, title='', state_label='n', show=True):
    plt.plot(moment_times, data_var, '--k', lw=2, label="data")
    if theory_var is not None:
        plt.plot(moment_times, theory_var, '--r', lw=2, label="generating")
    plt.title(title)
    plt.xlabel('time')
    plt.ylabel('Var(%s)' % state_label)
    plt.legend()
    plt.savefig(FOLDER_OUTPUT + os.sep + 'traj_%s_var_%s.png' % (model, state_label))
    if show:
        plt.show()
    return plt.gca()


def plot_hist(moment_times, data_observations, step, model, state_label='n', theory_mean=None, theory_var=None, show=True):
    time = moment_times[step]
    hist_data = data_observations[step, :]
    hist_bins = np.arange(0, hist_data.max() + 1.5) - 0.5
    num_traj = len(hist_data)
    count, bins, _ = plt.hist(data_observations[step, :]+0.01, bins=hist_bins)  # TODO cleanup
    if theory_mean is not None:
        assert theory_var is not None
        mu = theory_mean[step]
        sigma = np.sqrt(theory_var[step])
        normal_at_bins = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2))
        plt.plot(bins, num_traj * normal_at_bins, '--k')
    plt.title('Histogram (%d traj) for %s, %s at step:%d, time:%.2f' % (num_traj, model, state_label, step, time))
    plt.savefig(FOLDER_OUTPUT + os.sep + 'hist_%s_%s_%d_%.2f.png' % (model, state_label, step, time))
    if show:
        plt.show()
    return plt.gca()


def plot_estimation(moment_times, estimate_data, model, params, theory_curves, est_label='x', show=True):
    # TODO care normalize var and implement theory
    print "CARE NORMALIZE EST VAR IN PLOT"
    p = params

    def true_model_value():
        if est_label == 'x':
            model_value = p.k_on * p.c / p.k_off
        elif est_label == 'c':
            model_value = p.c
        else:
            assert est_label == 'k_off'
            model_value = p.k_off
        return model_value

    def est_theory(t):
        if model == 'mode_1':
            assert est_label == 'x'
            mu_n = theory_curves['mean_n'][idx]
            x_star = mu_n/(p.k_p*t - mu_n)
            est_mu_t = x_star
            est_var_t = (1/mu_n) * ((1 + x_star)**2 + 2 * p.k_p / p.k_off)

        elif model == 'mode_2':
            assert est_label == 'x'
            print "WARNING, est_theory() in plot_estimation() not implemented for mode_2"
            assert 1 == 2
            est_mu_t = -1
            est_var_t = -1

        else:
            assert est_label in ['c', 'k_off']
            mu_n = theory_curves['mean_n'][idx]
            mu_m = theory_curves['mean_m'][idx]
            k_off_star = (p.k_p / p.k_off) * mu_m / mu_n
            c_star = k_off_star * mu_n / (p.k_p * t - mu_n)
            x_star = c_star * p.k_on / k_off_star
            factor = (1 + x_star) / (x_star * k_off_star * t)
            if est_label == 'c':
                est_mu_t = c_star
                term_1 = (1 + 3 * x_star ** 2) / (1 + x_star) ** 2
                term_2 = (k_off_star / p.k_p) * x_star * (1 + x_star - x_star**2 + x_star**3)
                est_var_t = c_star ** 2 * factor * (term_1 + term_2)
            else:
                est_mu_t = k_off_star
                term_1 = (3 + x_star ** 2) / (1 + x_star) ** 2
                term_2 = (k_off_star / p.k_p) * ((1 + x_star**2) - 2 * (1 - x_star))
                est_var_t = k_off_star ** 2 * factor * (term_1 + term_2)
        return est_mu_t, est_var_t

    estimate_mean_at_t = np.zeros(moment_times.shape)
    estimate_var_at_t = np.zeros(moment_times.shape)
    theory_mean_at_t = np.zeros(moment_times.shape)
    theory_var_at_t = np.zeros(moment_times.shape)
    model_value = true_model_value()
    for idx in xrange(len(moment_times)):
        est_mu_t, est_var_t = est_theory(moment_times[idx])
        theory_mean_at_t[idx] = est_mu_t
        theory_var_at_t[idx] = est_var_t
        estimate_mean_at_t[idx] = np.mean(estimate_data[idx, :])
        estimate_var_at_t[idx] = np.var(estimate_data[idx, :])

    plt.close()
    fig = plt.figure(figsize=(8, 6))
    plt.clf()
    plt.suptitle('Estimation performance over time (%d traj)' % estimate_data.shape[-1])
    # plot mean over time
    ax1 = fig.add_subplot(121)
    ax1.plot(moment_times, estimate_mean_at_t, '-b', label='sample mean of estimate')
    ax1.plot(moment_times, theory_mean_at_t, '--b', label='theory')
    ax1.plot(moment_times, [model_value for _ in moment_times], '--k', label='true value')
    ax1.set_title('mean of estimate for %s' % est_label)
    ax1.set_xlabel('t')
    ax1.set_ylabel('<x_guess>')
    ax1.legend()
    # plot var over time
    ax2 = fig.add_subplot(122)
    ax2.plot(moment_times, estimate_var_at_t, '-b', label='sample var of estimate')
    ax2.plot(moment_times, theory_var_at_t, '--b', label='theory')
    ax2.set_title('variance of estimate for %s' % est_label)
    ax2.set_xlabel('t')
    ax2.set_ylabel('Var(x_guess)')
    ax2.legend()
    # save and show
    plt.savefig(FOLDER_OUTPUT + os.sep + 'est_%s_%s.png' % (model, est_label))
    if show:
        plt.show()
    return plt.gca()
