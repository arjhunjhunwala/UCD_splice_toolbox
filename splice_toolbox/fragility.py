'''
Aditya Jhunjhunwala
R0 - 04/15/2022
R1 - 07/11/2022 - Figure font and sizes changed

function to generate fracture fragility of the splice given inputs
'''

# %% inbuilt and custom imports
import os
import sys
import matplotlib.pyplot as plt

from splice_toolbox.additional_functions import *
from splice_toolbox.capacity import *
from splice_toolbox.demand import *

font_size = 12
plt.rcParams.update({'font.family': 'Times New Roman', 'mathtext.fontset': 'dejavuserif'})
lineStyle = ['-', (0, (5, 2)), (0, (7, 3, 3, 5)), (0, (5, 5)), (0, (2, 2))]
color = ['k', 'r', 'b']
marker_type = ['o', '^', 's', '<']


def fragility(modelName, upper_section, lower_section, flange_pen, story_ht, splice_ht, sigma_ys, sigma_yw, CVN, T_CVN,
              T_LAST, Kmin, Mz1, Mx1, Mz2, Mx2, Py, nCore=4, flag_abaqus=False, flag_geomUnc=False, flag_useExist=True, flag_generateImg=False):
    # ----------------------------------------------------------------------------------------------------------------------
    # %% Creating text file for input to abaqus file generation
    # ----------------------------------------------------------------------------------------------------------------------
    with open('splice_toolbox\\aisc_shapes.json') as f:
        aisc_shapes = json.load(f)
    shape_array = np.array(aisc_shapes['AISC_Manual_Label'])

    idx_us = np.where(shape_array == upper_section)[0][0]
    idx_ls = np.where(shape_array == lower_section)[0][0]

    depth_top_sec = aisc_shapes['d'][idx_us]
    tu = aisc_shapes['tf'][idx_us]
    tl = aisc_shapes['tf'][idx_ls]
    twu = aisc_shapes['tw'][idx_us]
    twl = aisc_shapes['tw'][idx_ls]
    bf = aisc_shapes['bf'][idx_us]
    bw = aisc_shapes['T'][idx_us]
    rtran = (depth_top_sec - 2.0 * tu - bw) / 2.0
    a_crack = (1.0 - flange_pen / 100) * tu
    aw_crack = 0.5 * twu  # the web penetration is fixed in the toolbox

    if tu > tl:
        print('[ERROR] The upper section is larger than the lower section. \nRerun the program with correct inputs.')
        sys.exit()

    if story_ht <= splice_ht:
        print('[ERROR] Splice location is greater than/ equal to story ht. \nRerun the program with correct inputs.')
        sys.exit()

    # geometry json file
    geom_keys = ['modelName', 'tu', 'tl', 'bf', 'a_crack', 'twu', 'twl', 'bw', 'aw_crack', 'rtran', 'story_ht',
                 'splice_ht', 'sigma_ys', 'sigma_yw', 'flag_generateImg']
    geom_values = [modelName, tu, tl, bf, a_crack, twu, twl, bw, aw_crack, rtran, story_ht * 12.0, splice_ht * 12.0, sigma_ys, sigma_yw, flag_generateImg]
    geomDataFile = 'geomData.json'
    createJsonFile(geom_keys, geom_values, geomDataFile)

    # loading json file
    load_keys = ['Mx1', 'Mz1', 'Py', 'Mx2', 'Mz2']
    load_values = [Mx1, Mz1, Py, Mx2, Mz2]
    loadDataFile = 'loadData.json'
    createJsonFile(load_keys, load_values, loadDataFile)

    # saving the weld material and the steel material true stress strain properties as text file
    generateTrueStressStrainRO(sigma_ys, 29000.0, 10, 'steelPlasticTrue.txt')
    generateTrueStressStrainRO(sigma_yw, 29000.0, 10, 'weldPlasticTrue.txt')

    # ----------------------------------------------------------------------------------------------------------------------
    # %% Check for entire flange in tension- removed in R1 of toolbox
    # ----------------------------------------------------------------------------------------------------------------------

    # Processing loads
    # sigma_avg, Gx, grady, flag_ok = processLoadInputRS(tu, bf, twu, bw, rtran, story_ht, splice_ht, Mx1, Mz1, Py, Mx2,
    #                                                    Mz2)
    # if flag_ok == False:
    #     print(
    #         'The toolbox is valid if one of the flange is entirely in tension. \nPartial tension and compression cases '
    #         'are not dealt with in this version -  to be updated soon!.')
    #     sys.exit()

    # ----------------------------------------------------------------------------------------------------------------------
    # %% Create the variables for monte carlo
    # ----------------------------------------------------------------------------------------------------------------------

    nMC = 1000
    npts = 30

    # ----------------------------------------------------------------------------------------------------------------------
    # %% Generating the fragility with and without the uncertainty
    # Step 1 : Generate the full model (abaqus)
    # Step 2: Run the full model (abaqus)
    # Step 3: Post process the odb file to obtain J integral
    # Step 4: generate the fragility

    # Note: before Step 1 to 3 - check if existing files are there for the same model name
    # ----------------------------------------------------------------------------------------------------------------------

    if flag_abaqus == True:
        outputODB = modelName + '.odb'
        outputFileNodeSet = modelName + '-NodeCoord.txt'
        outputFileJ = modelName + '-J-2.txt'
        outputFileStress = modelName + '-Stress.txt'

        if flag_useExist and (os.path.exists(outputODB) or (os.path.exists(outputFileNodeSet) and os.path.exists(outputFileJ))):
            # pop up to check if user wants to use existing results or run a new simulation
            print('Using existing simulation run.')
            if os.path.exists(outputFileNodeSet) and os.path.exists(outputFileJ):
                pass
            else:
                # Post-process the odb to obtain the node set and J integral results
                fileName = 'splice_toolbox\\postProcess.py'
                command = 'start /wait cmd /c abaqus cae noGUI="' + fileName + '"'
                os.system(command)
        else:
            # Generate the model .inp file
            fileName = 'splice_toolbox\\fullModel.py'
            print('starting generation of inp file')
            command = 'start /wait cmd /c abaqus cae noGUI="' + fileName + '"'
            os.system(command)
            # Analyse the model to obtain the .odb file
            command = 'start /wait cmd /c call abaqus job=' + modelName + ' cpus=' + str(int(nCore)) + ' interactive'
            os.system(command)
            # Post-process the odb to obtain the node set and J integral results
            fileName = 'splice_toolbox\\postProcess.py'
            command = 'start /wait cmd /c abaqus cae noGUI="' + fileName + '"'
            os.system(command)

        z_over_L, s, kData, bf = process_j(outputFileNodeSet, outputFileJ, outputFileStress)

        Pf_noMC = np.zeros((npts, len(CVN)))
        Pf_MC = np.zeros((npts, len(CVN)))

        for i in range(len(CVN)):
            Kmed = CVNtoK1dKimberly(CVN[i], 1)
            K0 = K1dtoK0(Kmed, T_CVN[i], T_LAST[i], sigma_ys)
            Kmed_MC = CVNtoK1dKimberly(CVN[i], nMC)
            K0_MC = K1dtoK0(Kmed_MC, T_CVN[i], T_LAST[i], sigma_ys)
            stress_noMC, Pf_noMC[:, i] = fracture_fragility_individual(z_over_L, s, kData, bf, K0, Kmin, nPts=npts)
            stress_MC, Pf_MC[:, i] = fracture_fragility_individual(z_over_L, s, kData, bf, K0_MC, Kmin, nPts=npts)

    # ----------------------------------------------------------------------------------------------------------------------
    # %% Generating the fragility with and without the uncertainty - using the response surface functions - removed in R1
    # ----------------------------------------------------------------------------------------------------------------------

    # sigma_avg, Gx, grady, flag_ok = processLoadInputRS(tu, bf, twu, bw, rtran, story_ht, splice_ht, Mx1, Mz1, Py, Mx2,
    #                                                    Mz2)

    if flag_abaqus == False:
        print('No abaqus feature coming soon...')
        sys.exit()
        # LM = np.linspace(0.0, 1.0, npts + 1)
        # Pf_MC = np.zeros((len(LM), len(CVN)))
        # Pf_noMC = np.zeros((len(LM), len(CVN)))
        # if flag_geomUnc == 1:
        #     tu_MC = flange_thk(tu, nMC)
        #     tl_MC = flange_thk(tl, nMC)
        #     a_MC = flange_thk(a_crack, nMC)
        # for j in range(len(LM)):
        #     stress = sigma_avg * LM[j]
        #     for i in range(len(CVN)):
        #         Kmed = CVNtoK1dKimberly(CVN[i], 1)
        #         K0 = K1dtoK0(Kmed, T_CVN[i], T_LAST[i], sigma_ys)
        #         Kmed_MC = CVNtoK1dKimberly(CVN[i], nMC)
        #         K0_MC = K1dtoK0(Kmed_MC, T_CVN[i], T_LAST[i], sigma_ys)
        #         if flag_geomUnc == 1:
        #             Pf_MC[j, i] = fracture_probability_monte_carlo_rs(stress, Gx, grady, tu_MC, tl_MC, a_MC, bf, K0_MC, Kmin,
        #                                                            nMC)
        #         else:
        #             Pf_MC[j, i] = fracture_probability_monte_carlo_rs(stress, Gx, grady, tu, tl, a_crack, bf, K0_MC, Kmin, nMC)
        #         Pf_noMC[j, i] = fracture_probability_monte_carlo_rs(stress, Gx, grady, tu, tl, a_crack, bf, K0, Kmin, 1)

    # ----------------------------------------------------------------------------------------------------------------------
    # %% Plot the fragility
    # # fragility plots without the Monte carlo are not displayed in R1
    # ----------------------------------------------------------------------------------------------------------------------

    # fig1 = plt.figure()
    # axes1 = fig1.subplots()

    fig2 = plt.figure()
    axes2 = fig2.subplots()

    if flag_abaqus == True:
        for i in range(len(CVN)):
            # axes1.plot(stress_noMC, Pf_noMC[:,i], lw=2.0, label='CVN={:.1f} @ {:.1f}F, LAST={:.1f}F'.format(CVN[i], T_CVN[i], T_LAST[i]))
            plot_label = 'CVN={:.1f} @ {:.1f}F, LAST={:.1f}F'.format(CVN[i], T_CVN[i], T_LAST[i])
            axes2.plot(stress_MC, Pf_MC[:,i], lw=2.0, label=plot_label)
            # axes1.set_xlim(left=0.0, right=s[-1])
            axes2.set_xlim(left=0.0, right=s[-1])
        data_out_noMC = np.hstack((stress_noMC.reshape(-1,1), Pf_noMC))
        data_out_MC = np.hstack((stress_noMC.reshape(-1,1), Pf_MC))
    # else:
    #     for i in range(len(CVN)):
    #         axes1.plot(LM, Pf_noMC[:, i], 'k', lw=2.0, label='CVN={:.1f} @ {:.1f}F, LAST={:.1f}F'.format(CVN[i], T_CVN[i], T_LAST[i]))
    #         axes2.plot(LM, Pf_MC[:, i], 'r', lw=2.0, label='CVN={:.1f} @ {:.1f}F, LAST={:.1f}F'.format(CVN[i], T_CVN[i], T_LAST[i]))
    #         axes1.set_xlim(left=0.0, right=LM[-1])
    #         axes2.set_xlim(left=0.0, right=LM[-1])
    #     data_out_noMC = np.hstack((LM.reshape(-1, 1), Pf_noMC))
    #     data_out_MC = np.hstack((LM.reshape(-1, 1), Pf_MC))

    # axes1.set_xlabel('Load factor', fontsize=font_size)
    # axes1.set_ylabel('Probability of fracture', fontsize=font_size)
    # axes1.set_title('Top section: ' + upper_section + '; Bottom section: ' + lower_section + '\nFlange PJP: ' + str(
    #     flange_pen) + ' %' + '\nNo Monte Carlo')
    # axes1.legend(fontsize=font_size)
    # axes1.grid(which='major')
    # axes1.grid(which='minor', color='#EEEEEE', linewidth=0.5)
    # axes1.minorticks_on()
    # axes1.set_ylim([0, 1])

    axes2.set_xlabel('Load factor', fontsize=font_size)
    axes2.set_ylabel('Probability of fracture', fontsize=font_size)
    axes2.set_title('Top section: ' + upper_section + '; Bottom section: ' + lower_section + '\nFlange PJP: ' + str(
        flange_pen) + ' %', fontsize=font_size)
    axes2.legend(loc=2, fontsize=font_size)
    axes2.grid(which='major')
    axes2.grid(which='minor', color='#EEEEEE', linewidth=0.5)
    axes2.minorticks_on()
    axes2.set_ylim([0, 1])

    return fig2, data_out_noMC, data_out_MC


def fragility_kimberly(upper_section, lower_section, flange_pen, story_ht, splice_ht, sigma_ys, sigma_yw, CVN, T_CVN,
                       T_LAST, Kmin, Mz1, Mx1, Mz2, Mx2, Py):
    with open('splice_toolbox\\aisc_shapes.json') as f:
        aisc_shapes = json.load(f)
    shape_array = np.array(aisc_shapes['AISC_Manual_Label'])

    idx_us = np.where(shape_array == upper_section)[0][0]
    idx_ls = np.where(shape_array == lower_section)[0][0]

    depth_top_sec = aisc_shapes['d'][idx_us]
    tu = aisc_shapes['tf'][idx_us]
    tl = aisc_shapes['tf'][idx_ls]
    twu = aisc_shapes['tw'][idx_us]
    twl = aisc_shapes['tw'][idx_ls]
    bf = aisc_shapes['bf'][idx_us]
    bw = aisc_shapes['T'][idx_us]
    rtran = (depth_top_sec - 2.0 * tu - bw) / 2.0
    a_crack = (1.0 - flange_pen / 100) * tu
    aw_crack = 0.5 * twu  # the web penetration is fixed in the toolbox

    # ----------------------------------------------------------------------------------------------------------------------
    # %% Check for entire flange in tension
    # ----------------------------------------------------------------------------------------------------------------------

    # Processing loads
    sigma_avg, Gx, grady, flag_ok = processLoadInputRS(tu, bf, twu, bw, rtran, story_ht, splice_ht, Mx1, Mz1, Py, Mx2,
                                                       Mz2)
    if flag_ok == False:
        print(
            'The toolbox is valid if the entire tension flange is in tension. \nPartial tension and compression cases '
            'are not dealt with in this version.')
        sys.exit('')

    sigma_max = sigma_avg*(1+grady/100.0)

    # ----------------------------------------------------------------------------------------------------------------------
    # %% Create the variables for monte carlo
    # ----------------------------------------------------------------------------------------------------------------------

    nMC = 1000
    npts = 30

    # Generate capacity instances for Monte Carlo
    Kmed = CVNtoK1dKimberly(CVN, 1)
    K0 = K1dtoK0(Kmed, T_CVN, T_LAST, sigma_ys)
    Kmed_05 = CVNtoK1dKimberly(CVN,1)*(1.0 + stats.norm.ppf(0.05)*0.11)
    K0_05 = K1dtoK0(Kmed_05, T_CVN, T_LAST, sigma_ys)
    Kmed_MC = CVNtoK1dKimberly(CVN, nMC)
    K0_MC = K1dtoK0(Kmed_MC, T_CVN, T_LAST, sigma_ys)

    # ----------------------------------------------------------------------------------------------------------------------
    # %% Generating the fragility with and without the uncertainty
    # Step 1: determine demand from the equations developed by kimberly
    # Step 2: use the 95% K value of determine if failure takes place or not
    # ----------------------------------------------------------------------------------------------------------------------

    K_arr = np.arange(0,250)
    sigma_arr = toughness_demand_kimberly(tu, tl, a_crack, K_arr)

    LM = np.linspace(0.0, 1.0, npts + 1)
    sigma_avg_arr = sigma_avg*LM
    sigma_max_arr = sigma_max*LM
    Pf_avg = np.zeros(len(LM))
    Pf_max = np.zeros(len(LM))

    data_out = np.vstack((LM, Pf_avg, Pf_max)).transpose()

    return data_out