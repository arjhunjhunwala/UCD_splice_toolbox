'''
Aditya Jhunjhunwala
04/13/2022

additional function required for the toolbox
'''

import json
import numpy as np


def createJsonFile(keys, values, fileName):
    """
    Creates a json file with single dictionary give the keys and values as list or ndarray
    """
    data_dict = dict(zip(keys, values))
    with open(fileName, 'w') as f:
        json.dump(data_dict, f)


def generateTrueStressStrainRO(fy, E, n, fileName):
    '''
    function to generate true stress strain of the ramberg-osgood material
    given the fy, E, nu and hardening parameter for the RO material
    the plastic stress and strain for abaqus input is stores in fileName
    the generated values go upto true strain of ~200%
    '''
    ey = fy / E
    incr = np.cumsum(
        np.hstack((np.zeros(1), np.ones(5) * 0.2, np.ones(20) * 0.5, np.ones(20), np.ones(10) * 2.0, np.ones(5) * 5.0)))
    sigma_engg = fy + incr
    strain_engg = ey * (sigma_engg / fy) ** n
    strain_true = np.log(1 + strain_engg)
    strain_true_plastic = np.round(strain_true - sigma_engg / E, 6)
    strain_true_plastic[0] = 0.0
    sigma_true = np.round((1 + strain_engg) * sigma_engg, 3)
    output = np.vstack((sigma_true, strain_true_plastic)).transpose()
    fmt = ('%.3f', '%.6f')
    np.savetxt(fileName, output, fmt=fmt)

# import pandas as pd
# filename = 'aisc-shapes-database-v15.0.xlsx'
# data = pd.read_excel(filename, sheet_name='Database v15.0')
# dataDict = data.reset_index().to_dict(orient='list')
# dataFile = 'aisc_shapes.json'
# with open(dataFile, 'w') as f:
#     json.dump(dataDict, f)
