import numpy as np
import os

"""
State for mode 1, 2: [bool, n]
    - bool == 1 implies BOUND, 0 implies UNBOUND
    - n is number of downstream molecules produced (due to being bound)

TODO
- more comments
- fix or remove theory curves for 'direct'
- convert model structures into class
"""

# project level
FOLDER_OUTPUT = "output"
if not os.path.exists(FOLDER_OUTPUT):
    os.makedirs(FOLDER_OUTPUT)

# model parameters
GLOB_C = 1.0
GLOB_K_ON = 1.0
GLOB_K_OFF = 50.0
GLOB_K_P = 80.0
GLOB_DEG_N = 0.0
GLOB_DEG_M = 0.0
GLOB_K_F = 10.0

# defined
GLOB_X = GLOB_C * GLOB_K_ON / GLOB_K_OFF
GLOB_PSS_BOUND = GLOB_X / (1 + GLOB_X)

# initial conditions (for single trajectory)
GLOB_N0 = 0.0
GLOB_M0 = 0.0
GLOB_BOUND_BOOL = 0
assert GLOB_BOUND_BOOL in [0, 1]

# misc
NUM_STEPS = 100

# models
VALID_MODELS = ['mode_1', 'mode_2', 'combined', 'kpr', 'two_ligand_kpr']
DEFAULT_MODEL = 'mode_1'

# model structures
NUM_RXN = {'mode_1': 4,
           'mode_2': 3,
           'combined': 5,
           'kpr': 7,
           'two_ligand_kpr': 12}
STATE_SIZE = {'mode_1': 2,
              'mode_2': 2,
              'combined': 3,
              'kpr': 4,
              'two_ligand_kpr': 5}
# init cond for each model
INIT_CONDS = {'mode_1': [GLOB_BOUND_BOOL, GLOB_N0],
              'mode_2': [GLOB_BOUND_BOOL, GLOB_M0],
              'combined': [GLOB_BOUND_BOOL, GLOB_N0, GLOB_M0],
              'kpr': [GLOB_BOUND_BOOL, 0, GLOB_N0, GLOB_M0],  # TODO handle init cond for kpr
              'two_ligand_kpr': [GLOB_BOUND_BOOL, GLOB_N0, GLOB_M0, GLOB_N0, GLOB_M0]}
# reaction event update dictionary for each model
UPDATE_DICTS = {
    'mode_1': {0: np.array([1.0, 0.0]),  # binding
               1: np.array([-1.0, 0.0]),  # unbinding
               2: np.array([0.0, 1.0]),  # production
               3: np.array([0.0, -1.0])},  # degradation n
    'mode_2': {0: np.array([1.0, 1.0]),  # bind + GPCR event
               1: np.array([-1.0, 0.0]),  # unbind
               2: np.array([0.0, -1.0])},  # degradation m
    'combined': {0: np.array([1.0, 0.0, 1.0]),  # binding + GPCR event
                 1: np.array([-1.0, 0.0, 0.0]),  # unbinding
                 2: np.array([0.0, 1.0, 0.0]),  # production of n
                 3: np.array([0.0, -1.0, 0.0]),  # degradation n
                 4: np.array([0.0, 0.0, -1.0])},  # degradation m
    'kpr': {0: np.array([1.0, 0.0, 0.0, 0.0]),  # binding
            1: np.array([-1.0, 0.0, 0.0, 0.0]),  # unbinding
            2: np.array([-1.0, 1.0, 0.0, 1.0]),  # kpr forward step + GPCR event
            3: np.array([0.0, -1.0, 0.0, 0.0]),  # fall off
            4: np.array([0.0, 0.0, 1.0, 0.0]),  # produce n
            5: np.array([0.0, 0.0, -1.0, 0.0]),  # degradation n
            6: np.array([0.0, 0.0, 0.0, -1.0])},  # degradation m
    'two_ligand_kpr':
    #			     [state, n1, m1, n2,  m2 ]
        {0: np.array([1.0, 0.0, 1.0, 0.0, 0.0]),  # ligand #1 binding
         1: np.array([3.0, 0.0, 1.0, 0.0, 0.0]),  # ligand #2 binding
         2: np.array([-1.0, 0.0, 0.0, 0.0, 0.0]),  # unbinding of ligand #1 from state 1
         3: np.array([0.0, 1.0, 0.0, 0.0, 0.0]),  # produce n1
         4: np.array([1.0, 0.0, 0.0, 0.0, 1.0]),  # kpr forward step + GPCR event
         5: np.array([-2.0, 0.0, 0.0, 0.0, 0.0]),  # unbinding of ligand #1 from P2
         6: np.array([0.0, 0.0, 0.0, 1.0, 0.0]),  # produce n2
         7: np.array([-3.0, 0.0, 0.0, 0.0, 0.0]),  # unbinding of ligand #2 from P1
         8: np.array([0.0, 1.0, 0.0, 0.0, 0.0]),  # produce n1
         9: np.array([1.0, 0.0, 0.0, 0.0, 1.0]),  # kpr forward step + GPCR event
         10: np.array([-4.0, 0.0, 0.0, 0.0, 0.0]),  # unbinding of ligand #2 from P2
         11: np.array([0.0, 0.0, 0.0, 1.0, 0.0])}  # produce n2
}
