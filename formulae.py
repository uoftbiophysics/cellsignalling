import numpy as np

from params import DEFAULT_PARAMS
from settings import DEFAULT_MODEL, VALID_MODELS


def theory_cumulants(moment_times, label, init_n=0.0, init_bound=0.0, p=DEFAULT_PARAMS, model=DEFAULT_MODEL):
    assert label in ["direct", "generating"]
    assert model in VALID_MODELS

    mean_vals = np.zeros(moment_times.shape[0])
    var_vals = np.zeros(moment_times.shape[0])
    if model == 'combined':
        cov_vals = np.zeros(moment_times.shape[0])

    # unpack params
    c, k_on, k_off, k_p, x, pss, r, d_n, d_m = p.unpack()
    # identities
    delta1 = init_bound - pss

    if label == "direct":
        C1 = pss ** 3 + 1 / r * (delta1 * pss + x / (1 + x) ** 3)
        C2 = pss ** 2 - delta1 * (1 - x) / (1 + x) ** 2
        C3 = delta1 * (1 / x) * (1 + x**2) / (1 + x) ** 2
        for idx, t in enumerate(moment_times):
            mean_vals[idx] = k_p * pss * t + init_n + k_p * delta1 * (1-np.exp(-r*t)) / r
            var_vals[idx] = 2 * k_p**2 * (C1 * (t - 1/r * (1 - np.exp(-r*t))) +
                                          C2 * t**2 / 2 +
                                          C3 * (1 - np.exp(-r*t) * (r*t + 1)) / r**2) + \
                            mean_vals[idx] - mean_vals[idx]**2
    else:
        init_p0 = 1 - init_bound
        for idx, t in enumerate(moment_times):
            x1 = c*k_on*k_p*(np.exp(-r*t) - 1 + r*t) / r**2
            x2 = k_p*(k_off - np.exp(-r*t)*k_off + k_on*c*(r*t)) / r**2
            mean_vals[idx] = x1 * (1 - init_bound) + x2 * init_bound

            var_term_1 = r ** 2 * ((1 - np.exp(-r * t)) * k_off * init_bound +
                                   c ** 2 * k_on ** 2 * (init_p0 + init_bound) * t +
                                   c * k_on * (k_off * init_bound * t + init_p0 * (-1 + np.exp(-r * t) + k_off * t)))
            var_term_2 = k_p * ((1 - np.exp(-r * t)) * k_off * init_bound + c ** 2 * k_on ** 2 * (init_p0 + init_bound) * t +
                                c * k_on * (k_off * init_bound * t + init_p0 * (-1 + np.exp(-r * t) + k_off * t))) ** 2
            var_term_3 = (c * k_on) ** 4 * (init_p0 + init_bound) * t ** 2 - 2 * k_off ** 2 * init_bound * (
                        np.exp(-r * t) - 1 + np.exp(-r*t) * k_off * t)
            var_term_4 = 2 * (k_on*c) ** 3 * t * (k_off * init_bound * t + init_p0 * (-1 + k_off * t))
            var_term_5 = 2 * c * k_off * k_on * (
                    init_p0 * (2*np.exp(-r*t) + k_off * t*np.exp(-r*t) + (-2 + k_off * t)) + init_bound * (
                            2*np.exp(-r*t) - k_off * t*np.exp(-r*t) + 2 * (-1 + k_off * t)))
            var_term_6 = (k_on * c) ** 2 * (
                    k_off * init_bound * t * (4 + k_off * t) + init_p0 * (
                            -2*np.exp(-r*t) + 2 * k_off * t*np.exp(-r*t) + (2 + (k_off*t)** 2)))
            var_term_7 = var_term_3 + var_term_4 + var_term_5 + var_term_6

            var_vals[idx] = (k_p / r ** 4) * (var_term_1 - var_term_2 + k_p * var_term_7)
    return mean_vals, var_vals
