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
