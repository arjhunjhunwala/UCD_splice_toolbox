"""
Splice fragility tool box

Aditya Jhunjhunwala, University of California Davis
arjhunjhunwala@ucdavis.edu

Amit Kanvinde, University of California Davis
kanvinde@ucdavis.edu

Revision history
R0 - 10/03/2022 - first version of toolbox - does not support the no abaqus feature

"""
# ----------------------------------------------------------------------------------------------------------------------
# %% Installing required modules
# ----------------------------------------------------------------------------------------------------------------------

import sys
import subprocess
import pkg_resources

python_path = sys.executable
subprocess.check_call([python_path, '-m', 'pip', 'install', '--upgrade', 'pip'], stdout=subprocess.DEVNULL)
required = {'numpy', 'scipy', 'matplotlib'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed
if missing:
    subprocess.check_call([python_path, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)

# ----------------------------------------------------------------------------------------------------------------------
# %% importing required modules
# ----------------------------------------------------------------------------------------------------------------------

import tkinter as tk
from tkinter import simpledialog, messagebox
from PIL import ImageTk, Image
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# import custom packages
from splice_toolbox.fragility import *

with open('splice_toolbox\\aisc_shapes.json') as f:
    aisc_shapes = json.load(f)
shape_array = np.array(aisc_shapes['AISC_Manual_Label'])

# ----------------------------------------------------------------------------------------------------------------------
# %% User Input GUI using tkinter
# ----------------------------------------------------------------------------------------------------------------------

# %% Root -> canvas -> frame -> widgets
root = tk.Tk()
root.option_add('*font', 'calibri 12')
root.geometry('1000x700')
root.title('Splice fracture fragility toolbox')
root['background'] = '#b8ccf4'

canvas = tk.Canvas(root, borderwidth=0, background="#b8ccf4")

window = tk.Frame(canvas, bg='#b8ccf4')
window.grid(row=0, column=0)


# %% Adding scrollbar to the canvas
def onFrameConfigure(canvas):
    """Reset the scroll region to encompass the inner frame"""
    canvas.configure(scrollregion=canvas.bbox("all"))


def _on_mousewheel(event):
    canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")


vsb = tk.Scrollbar(root, orient="vertical", command=canvas.yview)
canvas.configure(yscrollcommand=vsb.set)
vsb.pack(side="right", fill="y")
hsb = tk.Scrollbar(root, orient="horizontal", command=canvas.xview)
canvas.configure(xscrollcommand=hsb.set)
hsb.pack(side="bottom", fill="x")
canvas.pack(side="left", fill="both", expand=True)
canvas.create_window((4, 4), window=window, anchor="nw")
window.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))
window.bind_all("<MouseWheel>", _on_mousewheel)

# %% Values for Labels, fields, default values and units

pad_x = 5
pad_y = 5
label_geom1 = ['Model name', 'Story height ', 'Splice height from bottom']
label_geom2 = ['Upper section', 'Lower section']
label_geom3 = ['Flange PJP penetration']
label_geom4 = ['Web PJP penetration']
field_geom1 = ['modelName', 'story_ht', 'splice_ht']
field_geom2 = ['us', 'ls']
field_geom3 = ['flange_pen']
unit_geom1 = ['', 'ft', 'ft']
unit_geom2 = []
unit_geom3 = ['%']
unit_geom4 = ['%']
default_geom4 = [50.0]

label_material = ['Yield strength of steel', 'Yield strength of weld']
field_material = ['sigma_ys', 'sigma_yw']
default_material = [55.0, 65.0]
unit_material = ['ksi', 'ksi']

label_CVN = ['CVN of weld metal assembly', 'Temperature at CVN', 'LAST']
field_CVN = ['CVN', 'T_CVN', 'LAST']
default_CVN = ['10.0, 20.0', '70, 70', '0, 50']
unit_CVN = ['ft-lb', 'deg F', 'deg F']

label_loading = ['Mz1', 'Mx1', 'Py', 'Mz2', 'Mx2']
field_loading = ['Mz1', 'Mx1', 'Py', 'Mz2', 'Mx2']
unit_loading = ['kip in', 'kip in', 'kips', 'kip in', 'kip in']

label_abaqus = ['Use abaqus?']
field_abaqus = ['flag_abaqus']


def generateLabels(window, label_list, field_list, start_row, col_num, unit_list=[], default_val=[]):
    entries = {}
    if len(label_list) != len(field_list):
        print('Field list length not same as label list - using label list as dictionary keys')
        field_list = label_list

    for i in range(len(label_list)):
        field = field_list[i]
        # create labels
        la1 = tk.Label(window, text=label_list[i], bg='#b8ccf4')
        la1.grid(row=start_row + i, column=col_num, sticky='w', padx=pad_x, pady=pad_y)

        ent = tk.Entry(window)
        ent.grid(row=start_row + i, column=col_num + 1, sticky='w', padx=pad_x, pady=pad_y)
        if len(default_val) == len(label_list):
            ent.insert(0, default_val[i])
        entries[field] = ent

        if len(unit_list) != 0:
            la2 = tk.Label(window, text=unit_list[i], bg='#b8ccf4')
            la2.grid(row=start_row + i, column=col_num + 2, sticky='w', padx=pad_x, pady=pad_y)

    return entries


# %% Insert UC davis and peer symbol

start_row = 0
fileName = 'splice_toolbox\\ucdavis.png'
img = Image.open(fileName)
img = img.resize((175, 50), Image.ANTIALIAS)
img = ImageTk.PhotoImage(img)
la = tk.Label(window, image=img, bg='#b8ccf4')
la.image = img
la.grid(row=start_row, column=0, padx=pad_x)

fileName = 'splice_toolbox\\peer.png'
img = Image.open(fileName)
img = img.resize((100, 100), Image.ANTIALIAS)
img = ImageTk.PhotoImage(img)
la = tk.Label(window, image=img, bg='#b8ccf4')
la.image = img
la.grid(row=start_row, column=4, rowspan=2, padx=pad_x, sticky='e')

# %% Generating the labels, entries and buttons

# geometry entries 1 - model name and height
start_row += 1
entry_geom1 = generateLabels(window, label_geom1, field_geom1, start_row, 0, unit_geom1)
start_row += len(label_geom1)

# geometry entries 2 - section
label_list = label_geom2
field_list = field_geom2
col_num = 0
entry_section = {}
for i in range(len(label_list)):
    field = field_list[i]

    # create labels
    la1 = tk.Label(window, text=label_list[i], bg='#b8ccf4')
    la1.grid(row=start_row + i, column=col_num, sticky='w', padx=pad_x, pady=pad_y)

    ent = tk.StringVar()
    ent.set(shape_array[0])
    # create optionMenu for dropdown
    ent_menu = tk.OptionMenu(window, ent, *shape_array)
    ent_menu.config(bg='WHITE')
    ent_menu.grid(row=start_row + i, column=col_num + 1, sticky='w', padx=pad_x, pady=pad_y)
    entry_section[field] = ent

start_row += len(label_geom2)

# geometry entries 3 - flange penetration
entry_flangePen = generateLabels(window, label_geom3, field_geom3, start_row, 0, unit_geom3)
start_row += len(label_geom3)

# geometry entries 4 - web penetration - no entry !
la1 = tk.Label(window, text=label_geom4[0], bg='#b8ccf4')
la1.grid(row=start_row, column=col_num, sticky='w', padx=pad_x, pady=pad_y)
la2 = tk.Label(window, text=default_geom4[0], bg='#b8ccf4')
la2.grid(row=start_row, column=col_num + 1, sticky='w', padx=pad_x, pady=pad_y)
la3 = tk.Label(window, text=unit_geom4[0], bg='#b8ccf4')
la3.grid(row=start_row, column=col_num + 2, sticky='w', padx=pad_x, pady=pad_y)
start_row += len(label_geom4)

# material entries
entry_material = generateLabels(window, label_material, field_material, start_row, 0, unit_material, default_material)
start_row += len(label_material)
la1 = tk.Label(window, text='(Enter multiple values CVN properties and LAST as comma separated)', bg='#b8ccf4')
la1.grid(row=start_row, column=col_num, columnspan=3, sticky='w', padx=pad_x, pady=pad_y)
start_row += 1
entry_CVN = generateLabels(window, label_CVN, field_CVN, start_row, 0, unit_CVN, default_CVN)
start_row += len(label_CVN)

# loading entries
entry_loading = generateLabels(window, label_loading, field_loading, start_row, 0, unit_loading)
start_row += len(label_loading)

# abaqus entries
label_list = label_abaqus
field_list = field_abaqus
col_num = 0
entry_abaqus = {}
for i in range(len(label_list)):
    field = field_list[i]

    # create labels
    la1 = tk.Label(window, text=label_list[i], bg='#b8ccf4')
    la1.grid(row=start_row + i, column=col_num, sticky='w', padx=pad_x, pady=pad_y)

    ent = tk.StringVar()
    ent.set('No')
    # create optionMenu for dropdown
    ent_menu = tk.OptionMenu(window, ent, 'Yes', 'No')
    ent_menu.config(bg='WHITE')
    ent_menu.grid(row=start_row + i, column=col_num + 1, sticky='w', padx=pad_x, pady=pad_y)
    entry_abaqus[field] = ent

start_row += len(label_abaqus)

# %% sign convention image

fileName = 'splice_toolbox\\sign_convention.png'
img = Image.open(fileName)
img = img.resize((430, 500), Image.ANTIALIAS)
img = ImageTk.PhotoImage(img)
la = tk.Label(window, image=img, bg='#b8ccf4')
la.image = img
la.grid(row=2, rowspan=start_row, column=4, padx=20)


# %% Function to run analysis when next button is pressed

def runAnalysis():
    # 1) create a new popup and ask number of cores or to consider geometric uncertainty
    # default values
    flag_geomUnc = False
    flag_useExist = True
    flag_generateImg = False
    nCore = 4
    if entry_abaqus['flag_abaqus'].get() == 'Yes':
        flag_abaqus = True
        nCore = simpledialog.askinteger(title='Abaqus analysis parameters',
                                        prompt='Number of CPU cores to use for abaqus analysis?', parent=window)
        flag_useExist = messagebox.askyesno('Abaqus analysis parameters', 'Use existing analysis run if available?')
        flag_generateImg = messagebox.askyesno('Abaqus post processing parameters',
                                               'Generate stress contour images for all analysis steps?')
    else:
        flag_abaqus = False
        flag_geomUnc = messagebox.askyesno('Fitted equation consideration', 'Consider geometric uncertainties?')

    # 2) reading all input parameters for running the fragility
    # geometry inputs
    modelName = entry_geom1['modelName'].get()
    upper_section = entry_section['us'].get()
    lower_section = entry_section['ls'].get()
    flange_pen = float(entry_flangePen['flange_pen'].get())
    story_ht = float(entry_geom1['story_ht'].get())
    splice_ht = float(entry_geom1['splice_ht'].get())
    # material inputs
    sigma_ys = float(entry_material['sigma_ys'].get())
    sigma_yw = float(entry_material['sigma_yw'].get())
    CVN_str = entry_CVN['CVN'].get()
    T_CVN_str = entry_CVN['T_CVN'].get()
    T_LAST_str = entry_CVN['LAST'].get()
    CVN = [float(val.strip()) for val in CVN_str.split(',')]
    T_CVN = [float(val.strip()) for val in T_CVN_str.split(',')]
    T_LAST = [float(val.strip()) for val in T_LAST_str.split(',')]
    if len(T_LAST) != len(CVN) or len(T_CVN) != len(CVN):
        messagebox.showinfo(title=None, message='Enter same number of entries for CVN properties.')
        return
    print(CVN)
    Kmin = 20 / 1.0988435
    # loading input
    Mz1 = float(entry_loading['Mz1'].get())
    Mx1 = float(entry_loading['Mx1'].get())
    Mz2 = float(entry_loading['Mz2'].get())
    Mx2 = float(entry_loading['Mx2'].get())
    Py = float(entry_loading['Py'].get())

    # 3) run the fragility function after closing the function
    flag_run = messagebox.askyesno('Confirmation', 'Run Analysis?')
    if flag_run == False:
        return
    root.destroy()
    print('Running toolbox...')
    fig_MC, data_noMC, data_MC = fragility(modelName, upper_section, lower_section, flange_pen,
                                                     story_ht, splice_ht, sigma_ys, sigma_yw, CVN, T_CVN, T_LAST, Kmin,
                                                     Mz1, Mx1, Mz2, Mx2, Py, nCore, flag_abaqus,
                                                     flag_geomUnc, flag_useExist, flag_generateImg)
    outputFolder = os.path.join(os.getcwd(), modelName + '_Output')
    try:
        os.mkdir(outputFolder)
    except OSError as error:
        print('Note: Output folder for the model already exists. Same folder will be used.')

    figName_MC = os.path.join(outputFolder, modelName + '_fragility.png')
    fig_MC.savefig(figName_MC, dpi=600, bbox_inches='tight')
    # fig.show()
    plt.close()
    dataFileName_noMC = os.path.join(outputFolder, modelName + '_fragilityData_noMC.txt')
    dataFileName_MC = os.path.join(outputFolder, modelName + '_fragilityData_MC.txt')
    fmt = ['%.4f']
    head_txt1 = 'LM\t'
    head_txt2 = '  \t'
    head_txt3 = '  \t'
    for i in range(len(CVN)):
        head_txt1 += 'CVN={:.1f}\t'.format(CVN[i])
        head_txt2 += 'T={:.1f}\t'.format(T_CVN[i])
        head_txt3 += 'LAST={:.1f}\t'.format(T_LAST[i])
        fmt.append('%.8f')
    head_txt = head_txt1 + '\n' + head_txt2 + '\n' + head_txt3
    fmt = tuple(fmt)
    np.savetxt(dataFileName_noMC, data_noMC, fmt=fmt, header=head_txt, comments='')
    np.savetxt(dataFileName_MC, data_MC, fmt=fmt, header=head_txt, comments='')

    # display the table and figure in a new window
    root_new = tk.Tk()
    root_new.option_add('*font', 'calibri 12')
    root_new.title('Splice fracture fragility')
    chart_type = FigureCanvasTkAgg(fig_MC, root_new)
    chart_type.get_tk_widget().pack()


# %% next button and run analysis
btn_next = tk.Button(window, text='Next', command=runAnalysis)
btn_next.grid(row=start_row, column=col_num, sticky='w', padx=5, pady=5)
start_row += 1

# %%
root.mainloop()
