import numpy as np
from sswm import SM_C_H

full_param_name_list = ['T_GS', 'rf_alpha', 'rf_lambda', 'Eo', 'Td',
                        'LAI', 'RAI', 'hc', 'Zm',
                        'Ps0', 'b', 'Ks', 'n', 's_h', 's_fc',
                        'k_xl_max', 'Px50', 'Pg50']

soil_tex_list = [
    [11.4,  0.405, 1.28E-06, 0.482, 'CLAY'],
    [8.52,  0.63,  2.45E-06, 0.476, 'CLAY LOAM'],
    [7.12,  0.299, 6.30E-06, 0.42,  'SANDY CLAY LOAM'],
    [5.39,  0.478, 6.95E-06, 0.451, 'LOAM'],           # index 3 = default
    [4.9,   0.218, 3.47E-05, 0.435, 'SANDY LOAM'],
    [4.38,  0.09,  1.56E-04, 0.41,  'LOAMY SAND'],
    [4.05,  0.121, 1.76E-04, 0.395, 'SAND'],
]

def get_baseline_params(soil_i=3):
    b, Ps0, Ksat, n, soiltex = soil_tex_list[soil_i]
    Ps0  = Ps0 * -9.8067e-3          # m  → MPa
    Ksat = Ksat * 12 * 60 * 60       # m/s → m/half-day
    s_h  = (-10   / Ps0) ** (-1 / b)
    s_fc = (-0.03 / Ps0) ** (-1 / b)
    return {
        'Zm':       0.5,
        'Td':       12 * 60 * 60,
        'Eo':       0.004,
        'rf_alpha': 0.01,
        'rf_lambda':0.35,
        'b':        b,
        'Ps0':      Ps0,
        'Ks':       Ksat,
        'n':        n,
        's_h':      s_h,
        's_fc':     s_fc,
        'RAI':      10,
        'dr':       0.0005,
        'hc':       20,
        'k_xl_max': 0.0008,
        'Px50':    -2.5,
        'Pg50':    -1.5,
        'LAI':      2,
        'T_GS':     180,
        'a': 3,
        'g_plus':   0.0006,
        'a_logistic': 0.0003
    }