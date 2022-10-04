'''
Aditya Jhunjhunwala
04/13/2022

Functions for the demand part of the Monte Carlo
'''

import numpy as np
from scipy.interpolate import interp1d
from scipy import stats
import os

def read_inp_data(filename):
    # reads input file and gives as a list of list
    f = open(filename)
    a = []
    for line in f:
        stripped_line = line.strip()
        line_list = stripped_line.split()
        temp_list = np.array([float(val) for val in line_list])
        a.append(temp_list)
    return a


def processLoadInputRS(tu, bf, twu, bw, rtran, story_ht, splice_ht, Mx1, Mz1, Py, Mx2, Mz2):
    depth_top_sec = 2 * tu + 2 * rtran + bw
    Iz = 2 * ((bf * tu ** 3) / 12 + (bf * tu) * (bw / 2 + rtran + tu / 2) ** 2) + (twu * (bw + rtran * 2) ** 3) / 12
    Sz = Iz / (depth_top_sec / 2)
    Ix = 2 * ((tu * bf ** 3) / 12) + ((bw + rtran * 2) * twu ** 3) / 12
    Sx = Ix / (bf / 2)
    Area = 2 * (bf * tu) + (bw + 2 * rtran) * twu + (4 * rtran ** 2 - np.pi * rtran ** 2)

    splice_ht_top = story_ht - splice_ht

    # load calculation
    Vx1 = (Mz1 - Mz2) / story_ht
    Vz1 = (Mx2 - Mx1) / story_ht

    Mz_splice = abs(Mz1 - Vx1 * splice_ht_top)
    Mx_splice = abs(Mx1 + Vz1 * splice_ht_top)

    # splice flange in complete tension check
    x1 = (depth_top_sec / 2)
    x2 = x1 - tu

    sigma_max = Py / Area + Mz_splice * x1 / Iz
    sigma_min = Py / Area + Mz_splice * x2 / Iz
    sigma_avg = (sigma_min + sigma_max) / 2.0
    sigma_secondary = Mx_splice / Sx
    Gx = (sigma_max - sigma_avg) / sigma_avg
    grady = 100.0 * sigma_secondary / sigma_avg

    sigma_lowest = sigma_min - sigma_secondary
    if sigma_lowest <= 0.0:
        flag_tension = False
    else:
        flag_tension = True

    return sigma_avg, Gx, grady, flag_tension


def toughness_demand_kimberly(tu, tl, a_crack, K):
    if tu > tl:
        tu = tl
    eta = a_crack / tl
    xi = tu / tl

    A1 = 4.31
    A2 = 0.247
    A3 = 1.12
    B1 = 4.25
    B2 = 0.39
    C1 = 0.00000178
    C2 = -0.000499
    C3 = 0.0177
    C4 = -0.000000917
    C5 = 0.000272
    C6 = -0.00659
    D1 = 0.0003
    D2 = -0.0766
    D3 = 3.46

    f1 = A1 * eta ** 2 + A2 * eta + A3
    f2 = B1 * (xi - 1) ** 2 + B2 * (xi - 1) + 1.0
    g1 = (C1 * K ** 3 + C2 * K ** 2 + C3 * K) * eta + (C4 * K ** 3 + C5 * K ** 2 + C6 * K) + 1.0
    g2 = (D1 * K ** 2 + D2 * K + D3) * (xi - 1) ** 2 + 1.0

    factor = np.sqrt(np.pi * a_crack) * xi * f1 * f2 * g1 * g2
    sigma = K / factor

    return sigma


# ----------------------------------------------------------------------------------------------------------------------
# %% functions for direct estimation & processing of data
# ----------------------------------------------------------------------------------------------------------------------

def process_j(outputFileNodeSet, outputFileJ, outputFileStress=None):
    '''
    Processes the abaqus output files to numpy arrays

    :param outputFileNodeSet: file containing node data along the crack front
    :param outputFileJ: file containing J integral for different load factor along the crack front
    :param outputFileStress: file containing the S22 for different load factor along the crack front
    :return: (z_over_l, s, kData_processed, bf)
    '''

    nodeData = np.loadtxt(outputFileNodeSet)
    zCoord_raw = nodeData[:, 3]
    bf = max(zCoord_raw)
    nCoord = zCoord_raw.size
    idx = np.arange(0, nCoord, 2)
    z_over_l = zCoord_raw[idx] / max(zCoord_raw)

    jData = np.loadtxt(outputFileJ)
    s = jData[:, 0]
    kData_raw = np.sqrt(np.abs(jData[:, 1:]) * 29000 / (1 - 0.3 * 0.3))

    if outputFileStress is not None:
        if os.path.exists(outputFileStress):
            stressData = np.loadtxt(outputFileStress)
            stressAlongCrack = stressData[:, 1:]
            id_tension = np.where(stressAlongCrack < 0)
            kData_raw[id_tension] = 0.0
        else:
            print('[Warning] Output file for stress along the crack does not exist.')
            print('K field for compression region will also be considered in Pf calculations.')
            pass

    kData_processed = kData_raw[:, idx]
    return z_over_l, s, kData_processed, bf


def fracture_probability(z_over_L, K, K0, Kmin, bf):
    '''
    returns mean fracture probability for single instance of K vs. x/L for an array of K0

    :param z_over_L: array of z/L along the crack front
    :param K: array of K-field value (in ksi-sqrt(in)) for the z/L values
    :param K0: array of K0 for 1T specimen (in ksi-sqrt(in))
    :param Kmin: min threshold value of fracture toughness
    :param bf: width of flange (in in.)
    :return: single value of Probability of fracture
    '''

    zCoord = z_over_L * bf
    if isinstance(K0, np.ndarray) and K0.size > 1:
        K0 = np.reshape(K0, (-1, 1))
    Knorm = np.maximum(np.zeros(K.shape), K - Kmin) / (K0 - Kmin)
    Pf = 1 - np.exp(-np.trapz(Knorm ** 4, zCoord))
    Pf = np.average(Pf)
    return Pf


def fracture_fragility_individual(z_over_L, s, kData, bf, K0, Kmin, nPts=15):
    '''
    returns stress vs. Pf for a single geometry

    :param z_over_L: array of z/L along the crack front
    :param s: array of stress value or load factor
    :param kData: 2D array with K-field data at different s value and z/L (in ksi-sqrt(in))
    :param bf: width of flange (in in.)
    :param K0: array of K0 for 1T specimen (in ksi-sqrt(in))
    :param Kmin: min threshold value of fracture toughness
    :param nPts: number of equally spaced points between s[-1] and s[0] to be used for evaluation
    :return: (stress, Pf) array of probability of fracture values at different s values
    '''

    stress = np.linspace(s[0], s[-1], nPts)
    lin_fit = interp1d(s, kData, axis=0)
    K = lin_fit(stress)

    Pf = np.zeros(stress.shape)

    for j in range(len(stress)):
        Pf[j] = fracture_probability(z_over_L, K[j], K0, Kmin, bf)

    return stress, Pf


# ----------------------------------------------------------------------------------------------------------------------
# %% functions for generating random values for geometry for Monte Carlo
# ----------------------------------------------------------------------------------------------------------------------

def flange_thk(tf, n=1):
    # normal distribution
    # mean = 1.01, std Dev = 1% of mean
    bias = 1.01
    cov = 0.01
    if n == 1:
        tf_arr = np.array([tf])
    else:
        tf_arr = stats.norm.rvs(loc=bias * tf, scale=cov * tf, size=n)

    return tf_arr


def crack_width(a, n=1):
    # normal distribution
    # mean = 1.01, std Dev = 1% of mean
    bias = 1.01
    cov = 0.03
    if n == 1:
        a_arr = np.array([a])
    else:
        a_arr = stats.norm.rvs(loc=bias * a, scale=cov * a, size=n)

    return a_arr


# ----------------------------------------------------------------------------------------------------------------------
# %% Response surface based toughness demands - epistemic uncertainty - not to be used in the current Toolbox!
# ----------------------------------------------------------------------------------------------------------------------

def f_el_e(eta, xi, tu, bf, n=1):
    '''
    Elastic SIF for edge location
    input parameters can be array or individual
    if array and n!= 1 then error added to entire array mean
    :param eta:  a/tu
    :param xi: tu/tl
    :param tu: upper flange thickness
    :param bf: flange width
    :param n: number of samples if the input variables are individual
    :return: elastic SIF ratio - size: min(n,len(eta)) in case input is array
    '''
    xi_1 = 1 - xi
    bftf = bf / tu
    fel_e_mean = 1.0 - 0.593631 * xi_1 + (1.152114 - 1.412829 * xi_1 + 0.041723 * bftf) * eta + (
            0.072369 * bftf - 0.236603) * eta ** 2
    fel_e_stdDev = 0.045
    bias = 1.0
    if n == 1:
        fel_e = fel_e_mean * bias
    else:
        if isinstance(eta, int) or isinstance(eta, float):
            numSample = n
        else:
            numSample = len(eta)
        fel_e = fel_e_mean * bias + stats.norm.rvs(loc=0.0, scale=fel_e_stdDev, size=numSample)
    return fel_e


def f_el_c(eta, xi, tu, bf, n=1):
    # elastic SIF for centre
    xi_1 = 1 - xi
    etac = (eta * tu + 1.25) / (tu + 1.25)
    fel_c_mean = 1.0 - 0.8383 * xi_1 + (6.4483 - 0.9376 * tu) * etac + (-15.212 + 2.007 * tu) * etac ** 2 + (
            10.5027 - 1.1877 * tu) * etac ** 3
    fel_c_stdDev = 0.06
    bias = 1.012
    if n == 1:
        fel_c = fel_c_mean * bias
    else:
        if isinstance(eta, int) or isinstance(eta, float):
            numSample = n
        else:
            numSample = len(eta)
        fel_c = bias * fel_c_mean + stats.norm.rvs(loc=0.0, scale=fel_c_stdDev, size=numSample)
    return fel_c


def f_pl_e(eta, xi, tu, bf, s, n=1):
    # factor for plastic K values at edge
    s_app = s / (1.0 - 1.1 * eta) / 55.0
    fpl_e_mean = xi * (6.01779 * eta - 16.78847 * eta ** 2 + 17.12514 * eta ** 3 - 6.59487 * eta ** 4) * s_app ** 2
    fpl_e_cov = 0.15
    bias = 1.00
    if n == 1:
        fpl_e = fpl_e_mean * bias
    else:
        if isinstance(eta, int) or isinstance(eta, float):
            numSample = n
        else:
            numSample = len(eta)
        fpl_e = bias * fpl_e_mean * (1.0 + stats.norm.rvs(loc=0.0, scale=fpl_e_cov, size=numSample))
    return fpl_e


def f_pl_c(eta, xi, tu, bf, s, n=1):
    # factor for plastic K values at centre
    s_app = s / (1.0 - 1.1 * eta) / 55.0
    fpl_c_mean = xi * (17.12129 * eta - 121.16106 * eta ** 2 + 390.41342 * eta ** 3 - 650.45682 * eta ** 4 +
                       541.69554 * eta ** 5 - 178.88997 * eta ** 6 - 0.03592 * tu * eta) * s_app ** 2
    fpl_c_cov = 0.15
    bias = 1.00
    if n == 1:
        fpl_c = fpl_c_mean * bias
    else:
        if isinstance(eta, int) or isinstance(eta, float):
            numSample = n
        else:
            numSample = len(eta)
        fpl_c = bias * fpl_c_mean * (1.0 + stats.norm.rvs(loc=0.0, scale=fpl_c_cov, size=numSample))
    return fpl_c


def phi_Gx_e(eta, s_avg, Gx, n=1):
    # factor to modify the equivalent stress
    phi_mean = 1.0 - 1.0196992 * Gx + 0.826854 * Gx * eta - 0.0054162 * Gx * s_avg
    phi_stdDev = 0.01
    bias = 1.00
    if n == 1:
        phi = phi_mean * bias
    else:
        if isinstance(eta, int) or isinstance(eta, float):
            numSample = n
        else:
            numSample = len(eta)
        phi = bias * phi_mean + stats.norm.rvs(loc=0.0, scale=phi_stdDev, size=numSample)
    return phi


def phi_Gx_c(eta, s_avg, Gx, n=1):
    # factor to modify the equivalent stress
    phi_mean = 1.0 - 1.4118048 * Gx + 1.3187604 * Gx * eta - 0.0057795 * Gx * s_avg
    phi_stdDev = 0.01
    bias = 1.00
    if n == 1:
        phi = phi_mean * bias
    else:
        if isinstance(eta, int) or isinstance(eta, float):
            numSample = n
        else:
            numSample = len(eta)
        phi = bias * phi_mean + stats.norm.rvs(loc=0.0, scale=phi_stdDev, size=numSample)
    return phi


def phi_Gy_e(K, Gy):
    val250 = 1 - Gy + 0.8 * Gy ** 2
    if isinstance(K, int) or isinstance(K, float):
        if K <= 100.0:
            phi = 1.0
        else:
            phi = 1.0 - (K - 100.0) * (1.0 - val250) / 150.0
    else:
        if isinstance(Gy, int) or isinstance(Gy, float):
            val250 = np.ones(len(K)) * val250
        idx = np.where(K >= 100.0)[0]
        phi = np.ones(K.shape)
        phi[idx] = 1.0 - (K[idx] - 100.0) * (1.0 - val250[idx]) / 150.0
    return phi


def toughness_demand_rs(tu, tl, a_crack, bf, sigma, Gx, Gy, pl=1.0, ind_loc='e', n=1):
    # type - edge or e | centre or c
    # pl - plastic factor
    # sigma : sigma is the average stress in the flange
    # Gx and Gy are along thickness and along width gradient of stress in tne flange
    xi = tu / tl
    eta = a_crack / tu
    if ind_loc == 'e' or ind_loc == 'edge':
        phi_gx = phi_Gx_e(eta, sigma, Gx, n)
        ss = sigma * phi_gx  # stress star
        Kel = f_el_e(eta, xi, tu, bf, n) * ss * np.sqrt(np.pi * a_crack)
        fac_pl = f_pl_e(eta, xi, tu, bf, ss, n)
        K = Kel * np.sqrt(1 + pl * fac_pl ** 2)
        phi_gy = phi_Gy_e(K, Gy)  # correction for out of plane bending
        K = K * phi_gy
    else:
        phi_gx = phi_Gx_c(eta, sigma, Gx, n)
        ss = sigma * phi_gx  # stress star
        Kel = f_el_c(eta, xi, tu, n) * ss * np.sqrt(np.pi * a_crack)
        fac_pl = f_pl_c(eta, xi, tu, bf, ss, n)
        K = Kel * np.sqrt(1 + pl * fac_pl ** 2)
        # no correction for out of plane bending in centre !

    return K


def process_K1_mean_rs(s_arr, Gx, grad, tu, tl, a_crack, bf, pl_ind=1.0):
    # enter negative grad if you want low to high
    z_over_L = np.arange(0.025, 1.0, 0.025)
    if isinstance(grad, int) or isinstance(grad, float):
        Gy = np.ones(s_arr.shape) * grad / 100.0
    else:
        Gy = grad / 100.0

    stress_start = s_arr * (1 + Gy)
    stress_end = s_arr * (1 - Gy)
    stress_mid = s_arr

    Kstart = toughness_demand_rs(tu, tl, a_crack, bf, stress_start, Gx, Gy, pl=pl_ind, ind_loc='e')
    Kend = toughness_demand_rs(tu, tl, a_crack, bf, stress_end, Gx, Gy, pl=pl_ind, ind_loc='e')
    Kmid = toughness_demand_rs(tu, tl, a_crack, bf, stress_mid, Gx, Gy, pl=pl_ind, ind_loc='c')
    Kstart = np.reshape(Kstart, (-1, 1))
    Kend = np.reshape(Kend, (-1, 1))
    Kmid = np.reshape(Kmid, (-1, 1))
    lin_fit = interp1d(np.array([z_over_L[0], 0.5, z_over_L[-1]]), np.hstack((Kstart, Kmid, Kend)), axis=1)
    kData_processed = lin_fit(z_over_L)

    return z_over_L, s_arr, kData_processed, bf


def fracture_probability_monte_carlo_rs(stress, Gx, grad, tu, tl, a_crack, bf, K0, Kmin, n, pl_ind=1.0):
    '''
    :param stress: int or float
    :param Gx: int or float
    :param grad: int or float
    :param tu: int or float or ndarray
    :param tl: int or float or ndarray
    :param a_crack: int or float or ndarray
    :param bf: int or float
    :param K0: int or float or ndarray
    :param Kmin: constant
    :param n: n=1 fitting errors is not considered so the input has to be ndarray to consider monte carlo on other params
            n>1 fitting errors is considered. if input params are ndarray - n = min(n, len(ndarray))
    :param pl_ind: 1.0 to consider plastic J in calculation
    :return: Pf[n,] for the given loading condition
    '''
    # for a given loading i.e. stress, Gx and grad are int or float
    # if tu, tl and a_crack are fixed - enter n=nSim > 1 to consider the fitting errors else n = no of samples you want
    # if tu, tl and a_crack are arrays -
    # enter n = 1 to give out the mean values else n=2 for fitter errors to be included as well
    # K0 is either the shape of tu or a constant
    # bf is always constant
    # returns Pf for each realization
    z_over_L = np.arange(0.025, 1.0, 0.025)

    if isinstance(grad, int) or isinstance(grad, float):
        Gy = grad / 100.0
    else:
        print('insert gradient to be a single value (not array)')
        return 0

    stress_start = stress * (1 + Gy)
    stress_end = stress * (1 - Gy)
    stress_mid = stress

    Kstart = toughness_demand_rs(tu, tl, a_crack, bf, stress_start, Gx, Gy, pl=pl_ind, ind_loc='e', n=n)
    Kend = toughness_demand_rs(tu, tl, a_crack, bf, stress_end, Gx, Gy, pl=pl_ind, ind_loc='e', n=n)
    Kmid = toughness_demand_rs(tu, tl, a_crack, bf, stress_mid, Gx, Gy, pl=pl_ind, ind_loc='c', n=n)
    Kstart = np.reshape(Kstart, (-1, 1))
    Kend = np.reshape(Kend, (-1, 1))
    Kmid = np.reshape(Kmid, (-1, 1))
    lin_fit = interp1d(np.array([z_over_L[0], 0.5, z_over_L[-1]]), np.hstack((Kstart, Kmid, Kend)), axis=1)
    kData = lin_fit(z_over_L)

    zCoord = z_over_L * bf

    if isinstance(K0, list) or isinstance(K0, np.ndarray):
        K0 = np.reshape(K0, (-1, 1))

    Knorm = np.maximum(np.zeros(kData.shape), kData - Kmin) / (K0 - Kmin)
    rv = np.trapz(Knorm ** 4, zCoord, axis=1)
    Pf = 1 - np.exp(-rv)
    Pf = np.average(Pf)
    return Pf
