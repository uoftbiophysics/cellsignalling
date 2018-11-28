import matplotlib.pyplot as plt
import numpy as np

from params import DEFAULT_PARAMS
from settings import DEFAULT_MODEL, VALID_MODELS, NUM_RXN, UPDATE_DICTS, NUM_STEPS, INIT_CONDS, STATE_SIZE


def update_state(state_timeseries, rxn_idx, step, model):
    """
    Update the state based on the increment associated with the rxn_idx of specified model
    """
    current_state = state_timeseries[step]
    increment = UPDATE_DICTS[model][rxn_idx]
    state_timeseries[step + 1, :] = current_state + increment
    return state_timeseries


def propensities(state, model, params=DEFAULT_PARAMS):
    """
    Computes propensitiers of each rxn based on current state.
    Reaction labelling follows order defined in UPDATE_DICTS.
    """
    p=params
    assert model in VALID_MODELS
    propensities = np.zeros(NUM_RXN[model])
    if model == 'mode_1':
        propensities[0] = p.k_on * p.c * (1 - state[0])  # bind (0.0 if already bound)
        propensities[1] = p.k_off * state[0]           # unbind (0.0 if already unbound)
        propensities[2] = p.k_p * state[0]             # produce one n molecule
    elif model == 'mode_2':
        # TODO
        return 0
    elif model == 'combined':
        # TODO
        return 0
    return propensities


def simulate_traj(num_steps=NUM_STEPS, init_cond=None, model=DEFAULT_MODEL, params=DEFAULT_PARAMS):
    """
    Simulate single trajectory for num_steps
    - model: model string
    - params: params object storing k_on, k_off, c etc
    """
    # prep trajectories
    times = np.zeros(num_steps)
    traj = np.zeros((num_steps, STATE_SIZE[model]))
    # setup init conds
    if init_cond is not None:
        assert len(init_cond) == STATE_SIZE[model]
    else:
        init_cond = INIT_CONDS[model]
    traj[0, :] = init_cond
    times[0] = 0.0
    for step in xrange(num_steps-1):
        # generate two U[0,1] random variables
        r1, r2 = np.random.random(2)
        # generate propensities and their partitions
        alpha = propensities(traj[step], model, params=params)
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
        traj = update_state(traj, rxn_idx, step, model)
        times[step+1] = times[step] + tau
    return traj, times


def multitraj(num_traj, num_steps=NUM_STEPS, bound_fraction=0.0, model=DEFAULT_MODEL, params=DEFAULT_PARAMS):
    """
    Return:
    - traj_array: num_steps x STATE_SIZE[model] x num_traj
    - times_array: num_steps x num_traj
    """
    traj_array = np.zeros((num_steps, STATE_SIZE[model], num_traj))
    times_array = np.zeros((num_steps, num_traj))
    # prep init cond of varying bound states (such that average is bound_fraction
    init_cond_base = INIT_CONDS[model]
    draws = np.random.binomial(1, bound_fraction, num_traj)
    # simulate k trajectories
    for k in xrange(num_traj):
        init_cond_base[0] = draws[k]
        traj, times = simulate_traj(num_steps=num_steps, init_cond=init_cond_base, model=model, params=params)
        traj_array[:, :, k] = traj
        times_array[:, k] = times
    return traj_array, times_array


if __name__ == '__main__':
    # settings
    model = 'mode_1'
    num_traj = 200
    num_steps = 100
    init_bound = 1.0
    # compute
    traj_array, times_array = multitraj(num_traj, bound_fraction=init_bound, num_steps=num_steps)
    # plot trajectories
    for k in xrange(num_traj):
        times_k = times_array[:, k]
        traj_k = traj_array[:, 1, k]
        plt.plot(times_k, traj_k, '--', lw=0.5, alpha=0.5)
    # decorate
    plt.title('Mode 1 - %d trajectories' % num_traj)
    plt.xlabel('time')
    plt.ylabel('n')
    plt.legend()
    plt.show()
