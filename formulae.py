import numpy as np

from params import DEFAULT_PARAMS
from settings import DEFAULT_MODEL, VALID_MODELS


def theory_cumulants(moment_times, bound_fraction, method='generating', init_n=0.0, init_m=0.0,
                     p=DEFAULT_PARAMS, model=DEFAULT_MODEL):
    """
    Returns, as a dict, multiple output timeseries (aligned to moment_times) depending on the model:
        - mode_1: mean(n)(t), var(n)(t), None, None, None
        - mode_2: None, None, mean(m)(t), var(m)(t), None
        - combined: mean(n)(t), var(n)(t), mean(m)(t), var(m)(t), cov(n,m)(t)
    """
    assert method in ["direct", "generating"]
    assert model in VALID_MODELS

    # unpack params
    c, k_on, k_off, k_p, x, pss, r, d_n, d_m = p.unpack()
    assert d_n == 0.0  # no theory for this case
    assert d_m == 0.0  # no theory for this case

    # declare output
    theory_curves = {'mean_n': None,
                     'var_n': None,
                     'mean_m': None,
                     'var_m': None,
                     'cov_nm': None}

    if model == 'mode_1':
        theory_curves['mean_n'] = np.zeros(moment_times.shape[0])
        theory_curves['var_n'] = np.zeros(moment_times.shape[0])

        if method == "direct":
            delta1 = bound_fraction - pss
            C1 = pss ** 3 + 1 / r * (delta1 * pss + x / (1 + x) ** 3)
            C2 = pss ** 2 - delta1 * (1 - x) / (1 + x) ** 2
            C3 = delta1 * (1 / x) * (1 + x**2) / (1 + x) ** 2
            for idx, t in enumerate(moment_times):
                mean_n_val = k_p * pss * t + init_n + k_p * delta1 * (1-np.exp(-r*t)) / r
                var_n_val = 2 * k_p**2 * (C1 * (t - 1/r * (1 - np.exp(-r*t))) +
                                           C2 * t**2 / 2 +
                                           C3 * (1 - np.exp(-r*t) * (r*t + 1)) / r**2) + \
                            mean_n_val - mean_n_val**2
                theory_curves['mean_n'][idx] = mean_n_val
                theory_curves['var_n'][idx] = var_n_val
        else:
            init_p0 = 1 - bound_fraction
            for idx, t in enumerate(moment_times):
                x1 = c*k_on*k_p*(np.exp(-r*t) - 1 + r*t) / r**2
                x2 = k_p*(k_off - np.exp(-r*t)*k_off + k_on*c*(r*t)) / r**2
                mean_n_val = x1 * (1 - bound_fraction) + x2 * bound_fraction

                var_term_1 = r ** 2 * ((1 - np.exp(-r * t)) * k_off * bound_fraction +
                                       c ** 2 * k_on ** 2 * (init_p0 + bound_fraction) * t +
                                       c * k_on * (k_off * bound_fraction * t + init_p0 * (-1 + np.exp(-r * t) + k_off * t)))
                var_term_2 = k_p * ((1 - np.exp(-r * t)) * k_off * bound_fraction + c ** 2 * k_on ** 2 * (init_p0 + bound_fraction) * t +
                                    c * k_on * (k_off * bound_fraction * t + init_p0 * (-1 + np.exp(-r * t) + k_off * t))) ** 2
                var_term_3 = (c * k_on) ** 4 * (init_p0 + bound_fraction) * t ** 2 - 2 * k_off ** 2 * bound_fraction * (
                            np.exp(-r * t) - 1 + np.exp(-r*t) * k_off * t)
                var_term_4 = 2 * (k_on*c) ** 3 * t * (k_off * bound_fraction * t + init_p0 * (-1 + k_off * t))
                var_term_5 = 2 * c * k_off * k_on * (
                        init_p0 * (2*np.exp(-r*t) + k_off * t*np.exp(-r*t) + (-2 + k_off * t)) + bound_fraction * (
                                2*np.exp(-r*t) - k_off * t*np.exp(-r*t) + 2 * (-1 + k_off * t)))
                var_term_6 = (k_on * c) ** 2 * (
                        k_off * bound_fraction * t * (4 + k_off * t) + init_p0 * (
                                -2*np.exp(-r*t) + 2 * k_off * t*np.exp(-r*t) + (2 + (k_off*t)** 2)))
                var_term_7 = var_term_3 + var_term_4 + var_term_5 + var_term_6

                var_n_val = (k_p / r ** 4) * (var_term_1 - var_term_2 + k_p * var_term_7)
                theory_curves['mean_n'][idx] = mean_n_val
                theory_curves['var_n'][idx] = var_n_val

    elif model == 'mode_2':
        assert method == 'generating'
        theory_curves['mean_m'] = np.zeros(moment_times.shape[0])
        theory_curves['var_m'] = np.zeros(moment_times.shape[0])

        for idx, t in enumerate(moment_times):
            # TODO
            theory_curves['mean_m'][idx] = 0.0
            theory_curves['var_m'][idx] = 0.0

    else:
        assert method == 'generating'
        theory_curves['mean_n'] = np.zeros(moment_times.shape[0])
        theory_curves['var_n'] = np.zeros(moment_times.shape[0])
        theory_curves['mean_m'] = np.zeros(moment_times.shape[0])
        theory_curves['var_m'] = np.zeros(moment_times.shape[0])
        theory_curves['cov_nm'] = np.zeros(moment_times.shape[0])

        for idx, t in enumerate(moment_times):
            # TODO
            theory_curves['mean_n'][idx] = 0.0
            theory_curves['var_n'][idx] = 0.0
            theory_curves['mean_m'][idx] = 0.0
            theory_curves['var_m'][idx] = 0.0
            theory_curves['cov_nm'][idx] = 0.0

    return theory_curves
