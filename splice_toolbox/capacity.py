'''
Aditya Jhunjhunwala
04/13/2022

Functions for the capacity part of the Monte Carlo
'''

import numpy as np
from scipy import stats


# ----------------------------------------------------------------------------------------------------------------------
# %% CVN to K1d functions (distribution)
# ----------------------------------------------------------------------------------------------------------------------

def CVNtoK1dKimberly(CVN, nK=1):
    '''
    Function that returns K_Id value based on modified Barsom and Rolfe correlation (Stillmaker et al., 2016)

    :param CVN: Charpy impact energy in ft-lb
    :param nK: no of sample points if carrying out Monte Carlo
    :return: median dynamic fracture toughness capacity K_Id in ksi-sqrt(in)
    '''

    E = 29000 * 1000.0
    K1d_med = np.sqrt(7.576 * E * CVN) / 1000
    bias = 1.0
    cov = 0.11

    if nK == 1:
        K1d = np.array([K1d_med])
    else:
        K1d = stats.norm.rvs(loc=bias * K1d_med, scale=cov * K1d_med, size=nK)

    return K1d


# ----------------------------------------------------------------------------------------------------------------------
# %% K1d to K0 using temperature shift and the master curve median fit
# ----------------------------------------------------------------------------------------------------------------------

def K1dtoK0(K1d, Td, Tc, fy, strain_rate=0):
    '''
    Use Barsom and Rolfe temperature shift along with Master Curve to determine K0 at given strain rate and temperature

    :param K1d: Median dynamic fracture toughness at Td (ksi-sqrt(in); can be a ndarray of values
    :param Td: Temperature at which K1d is measured (in F); single value
    :param Tc: Temperature at which static fracture is to be evaluated (in F); single value
    :param fy: Yield stress of material for temperature shift evaluation (in ksi); single value
    :param strain_rate: Strain rate for temeprature shift evaluation; single value
    :return: K0 at Tc (in ksi-sqrt(in))
    '''

    fac = 1.0988435

    if strain_rate <= 10.0 ** -5:
        Tshift = -(215.0 - 1.5 * fy)  # in F
    elif strain_rate <= 10.0 ** -3 and strain_rate > 10.0 ** -5:
        T_inter_from_static = (150.0 - fy) * (10.0 ** -3) ** 0.17
        T_from_static = (np.log10(strain_rate) + 5.0) * T_inter_from_static / 2.0
        Tshift = -(215.0 - 1.5 * fy - T_from_static)  # in F
    elif strain_rate <= 10.0 and strain_rate > 10.0 ** -3:
        T_inter_from_static = min((150.0 - fy) * (strain_rate) ** 0.17, 215.0 - 1.5 * fy)
        Tshift = -(215.0 - 1.5 * fy - T_inter_from_static)  # in F
    else:
        Tshift = 0.0

    Tshifted = ((Td + Tshift) - 32.0) * 5 / 9  # in C

    K1d_SI = K1d * fac

    if isinstance(K1d, int) or isinstance(K1d, float):
        K1d_SI = max(K1d_SI, 30.5)
    else:
        K1d_SI[np.where(K1d_SI <= 30.0)] = 30.5

    Tmed = Tshifted - np.log((K1d_SI - 30.0) / 70.0) / 0.019
    K0 = 31.0 + 77.0 * np.exp(((Tc - 32.0) * 5 / 9 - Tmed) * 0.019)  # in MPa-m^0.5

    return K0 / fac


# ----------------------------------------------------------------------------------------------------------------------
# %% BS code correlations
# ----------------------------------------------------------------------------------------------------------------------

def CVNtoK0_BS1(CVN, T_CVN, Tc, flag=False, npts=1):
    """
    Determining K0 at Tc from CVN at T_CVN using BS 7910 J.2.1 method

    :param CVN: Charpy impact energy (in ft-lb)
    :param T_CVN: Temperature at which CVN is measures (in F)
    :param Tc: Temperature at which static fracture toughness, K0 is to be determined (in F)
    :param flag: True - for unbiased estimate, False - for original Lower Bound
    :param npts: Number of data points if performing Monte Carlo. Use only when flag is True
    :return: K0 at Tc (in ksi-sqrt(in))
    """

    fac = 1.0988435

    Tadj = 68.7  # bias
    stdDev = 31.0  # std dev around adjusted value

    CVN_SI = 1.35582 * CVN
    T_CVN_SI = (T_CVN - 32.0) * 5.0 / 9.0
    Tc_SI = (Tc - 32.0) * 5.0 / 9.0
    K_at_TCVN = 12 * np.sqrt(CVN_SI)  # in MPa-m^0.5
    if flag:
        if npts == 1:
            T0 = np.array([T_CVN_SI - np.log((K_at_TCVN - 30.0) / 70.0) / 0.019 - Tadj])
        else:
            T0 = T_CVN_SI - np.log((K_at_TCVN - 30.0) / 70.0) / 0.019 - Tadj + stats.norm.rvs(loc=0, scale=stdDev,
                                                                                              size=npts)
    else:
        T0 = np.array([T_CVN_SI - np.log((K_at_TCVN - 30.0) / 70.0) / 0.019])

    K0 = 31.0 + 77.0 * np.exp((Tc_SI - T0) * 0.019)  # in MPa-m^0.5
    return K0 / fac


def CVNtoK0_BS2(T27, Tc, flag=False, npts=1):
    '''
    Determining K0 at Tc using BS 7910 J.2.2 method

    :param T27: Temperature at CVN = 27J (in F)
    :param Tc: Temperature at which static fracture toughness, K0 is to be determined (in F)
    :param flag: True - for unbiased estimate, False - for original Lower Bound
    :param npts: Number of data points if performing Monte Carlo. Use only when flag is True
    :return: K0 at Tc (in ksi-sqrt(in))
    '''
    Tadj = 49.1
    stdDev = 32.8
    fac = 1.0988435
    T27_SI = (T27 - 32.0) * 5.0 / 9.0
    Tc_SI = (Tc - 32.0) * 5.0 / 9.0

    if flag:
        if npts == 1:
            T0 = np.array([T27_SI - 18 + 25 - Tadj])
        else:
            T0 = T27_SI - 18 + 25 - Tadj + stats.norm.rvs(loc=0, scale=stdDev, size=npts)
    else:
        T0 = np.array([T27_SI - 18 + 25])
    K0 = 31.0 + 77.0 * np.exp((Tc_SI - T0) * 0.019)  # in MPa-m^0.5
    return K0 / fac
