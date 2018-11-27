import matplotlib.pyplot as plt
import numpy as np

"""
State for mode 1, 2: [bool, n]
    - bool == 1 implies BOUND, 0 implies UNBOUND
    - n is number of downstream molecules produced (due to being bound)
    
TODO 
- protect global constants
- more comments
"""

# globals
num_rxn_mode1 = 3
state_size_mode1 = 2
state_size_mode2 = 2

# model parameters
c = 4.0
k_on = 5.0
k_off = 1.0
k_p = 2.0

# model identities
x = k_on * c / k_off
pss = x / (1 + x)
r = k_on * c + k_off

# initial conditions (for single trajectory)
n0 = 0.0
m0 = 0.0
p0 = 0
assert p0 in [0, 1]
INIT_COND_MODE1 = [p0, n0]

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


def multitraj(num_traj, num_steps=NUM_STEPS, bound_fraction=0.0):
    traj_array = np.zeros((num_steps, state_size_mode1, num_traj))
    times_array = np.zeros((num_steps, num_traj))

    if bound_fraction in [0.0, 1.0]:
        init_conds = [[bound_fraction, n0] for _ in xrange(num_traj)]
    else:
        draws = np.random.binomial(1, bound_fraction, num_traj)
        #print(draws)
        init_conds = [[draws[k], n0] for k in xrange(num_traj)]

    for k in xrange(num_traj):
        traj, times = simulate_traj(num_steps=num_steps, init_cond=init_conds[k])
        traj_array[:, :, k] = traj
        times_array[:, k] = times
    return traj_array, times_array


def get_state_at_t(traj, times, t):
    step = 0
    while times[step] < t:
        step += 1
    return traj[step, :]


def get_mean_var_timeseries(traj_array, times_array):
    num_traj = np.shape(traj_array)[-1]
    dt = np.mean(times_array[1,:])
    endtime = np.min(times_array[-1,:])
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
        var_vals[idx] = statesquaresum/num_traj -  mean_vals[idx]**2
    return mean_vals, var_vals, moment_times


def theory_cumulants(moment_times, label, init_n=0.0, init_p1=0.0):
    assert label in ["direct", "generating"]
    mean_vals = np.zeros(moment_times.shape[0])
    var_vals = np.zeros(moment_times.shape[0])
    # identities
    delta1 = init_p1 - pss

    if label == "direct":
        for idx, t in enumerate(moment_times):
            mean_vals[idx] = k_p * pss * t + init_n + k_p * delta1 * (1-np.exp(-r*t)) / r
            var_vals[idx] = k_p**2 * pss * t ** 2 + 2 * k_p * init_n * t + 2 * k_p * delta1 * (r*t-1+np.exp(-r*t)) / r**2 \
                            + mean_vals[idx] - mean_vals[idx]**2
    else:
        for idx, t in enumerate(moment_times):
            #assert init_p1 == pss  # var given only for init_p1==pss)
            x1 = c*k_on*k_p*(np.exp(-r*t) - 1 + r*t) / r**2
            x2 = k_p*(k_off - np.exp(-r*t)*k_off + k_on*c*(r*t)) / r**2
            mean_vals[idx] = x1*(1 - init_p1) + x2*init_p1
            init_p0 = 1-init_p1
            var_vals[idx] = (k_p/r**4) * (r**2 * ((1 - np.exp(-r*t))*k_off*init_p1 + c**2 *k_on**2 *(init_p0 + init_p1)*t + c*k_on*(k_off*init_p1*t + init_p0*(-1 + np.exp(-r*t) + k_off*t))) - k_p*((1 - np.exp(-r*t))*k_off*init_p1 + c**2 *k_on**2 *(init_p0+init_p1)*t + c*k_on*(k_off *init_p1*t + init_p0*(-1 + np.exp(-r*t) + k_off*t)))**2 + np.exp(-r*t)*k_p*((c*k_on)**4 * np.exp(r*t)*(init_p0 + init_p1)*t**2 - 2*k_off**2*init_p1*(1 - np.exp(r*t) + k_off*t) + 2*c**3*np.exp(r*t)*k_on**3 * t*(k_off*init_p1*t + init_p0*(-1 + k_off*t)) + 2*c*k_off*k_on*(init_p0*(2 + k_off*t + np.exp(r*t)*(-2 + k_off*t)) + init_p1*(2 - k_off*t + 2*np.exp(r*t)*(-1 + k_off*t))) + (k_on*c)**2 *(np.exp(r*t)*k_off*init_p1*t*(4 + k_off*t) + init_p0*(-2 + 2*k_off*t + np.exp(r*t)*(2 + k_off**2 * t**2)))))
           # var_vals[idx] = k_p*pss*t + 2*k_p**2*t/k_off*x/(1+x)**3*(1 + (np.exp(-r*t) - 1)/r)
    return mean_vals, var_vals


if __name__ == '__main__':
    # settings
    num_traj = 500
    init_bound = pss

    # compute
    traj_array, times_array = multitraj(num_traj, bound_fraction=init_bound)
    mean_vals, var_vals, moment_times = get_mean_var_timeseries(traj_array, times_array)
    mean_vals_direct, var_vals_direct = theory_cumulants(moment_times, "direct", init_p1=init_bound)
    mean_vals_gen, var_vals_gen = theory_cumulants(moment_times, "generating", init_p1=init_bound)

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
    plt.plot(moment_times, mean_vals_direct, 'b', lw=2, label="direct")
    plt.plot(moment_times, mean_vals_direct + np.sqrt(var_vals_direct), '--b', lw=2)
    plt.plot(moment_times, mean_vals_direct - np.sqrt(var_vals_direct), '--b', lw=2)
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
