import numpy as np

"""
State for mode 1, 2: [bool, n]
    - bool == 1 implies BOUND, 0 implies UNBOUND
    - n is number of downstream molecules produced (due to being bound)

TODO
- more comments
- fix or remove theory curves for 'direct'
- convert model structures into class
"""

# model parameters
GLOB_C = 4.0
GLOB_K_ON = 5.0
GLOB_K_OFF = 40.0
GLOB_K_P = 2.0
GLOB_DEG_N = 0.0
GLOB_DEG_M = 0.0

# initial conditions (for single trajectory)
GLOB_N0 = 0.0
GLOB_M0 = 0.0
GLOB_BOUND_BOOL = 0
assert GLOB_BOUND_BOOL in [0, 1]

# misc
NUM_STEPS = 100

# models
VALID_MODELS = ['mode_1', 'mode_2', 'combined']
DEFAULT_MODEL = 'mode_1'

# model structures
NUM_RXN = {'mode_1': 4,
           'mode_2': 3,
           'combined': 5}
STATE_SIZE = {'mode_1': 2,
              'mode_2': 2,
              'combined': 3}
# init cond for each model
INIT_CONDS = {'mode_1': [GLOB_BOUND_BOOL, GLOB_N0],
              'mode_2': [GLOB_BOUND_BOOL, GLOB_M0],
              'combined': [GLOB_BOUND_BOOL, GLOB_N0, GLOB_M0]}
# reaction event update dictionary for each model
UPDATE_DICTS = {
    'mode_1': {0: np.array([1.0, 0.0]),    # binding
               1: np.array([-1.0, 0.0]),   # unbinding
               2: np.array([0.0, 1.0]),    # production
               3: np.array([0.0, -1.0])},  # degradation n
    'mode_2': {0: np.array([1.0, 1.0]),      # bind + GPCR event
               1: np.array([-1.0, 0.0]),     # unbind
               2: np.array([0.0, -1.0])},    # degradation m
    'combined': {0: np.array([1.0, 0.0, 1.0]),     # binding + GPCR event
                 1: np.array([-1.0, 0.0, 0.0]),    # unbinding
                 2: np.array([0.0, 1.0, 0.0]),     # production of n
                 3: np.array([0.0, -1.0, 0.0]),    # degradation n
                 4: np.array([0.0, 0.0, -1.0])}    # degradation m
}
