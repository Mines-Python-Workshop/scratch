"""
Electrochemical impedence spectroscopy plotting functions
"""

import eclabfiles as ecf
import numpy as np
import glob

from colour import Color
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker as ticker

from matplotlib.pyplot import rc
rc("text", usetex=True)
rc("font", **{"family":"sans-serif","sans-serif":["Helvetica"]},size="14")
rc("text.latex",preamble=r"\usepackage{sfmath}")

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' GENERAL PLOTTING FUNCTIONS '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def import_file(path, fn, header_rows):
    '''
    Imports a data file

    Parameters
    ----------
    path : str
        File directory path
    fn : str
        File name
    header_rows : int
        Number of rows of metadata at beginning of file

    Returns
    -------
    list
        Each column in file returned as a separate list

    '''
    return np.loadtxt(path+fn, unpack=True, dtype=float, skiprows=header_rows)

#------------------------------------------------------------------------------
def import_dir(path, filetype=None):
    '''
    Imports a file directory

    Parameters
    ----------
    path : str
        File directory path
    filetype : str, optional
        Add file extension if you only want to import one type. The default is None.

    Returns
    -------
    list_files : list (str)
        Files names

    '''
    
    '''
    path (str) -- directory path
    filetype (str, optional) -- specify file extenstion
    
    returns list of files
    '''
    if filetype is not None:
        list_files = glob.glob(path + "/" + filetype)
    else:
        list_files = glob.glob(path + "/*")
    return list_files

#------------------------------------------------------------------------------
def gradient_gen(start_hex, end_hex, num):
    '''
    Generates color gradient

    Parameters
    ----------
    start_hex : str
        Hex code for first gradient color, format "#000000"
    end_hex : str
        Hex code for final gradient color, format "#000000"
    num : int
        Number of colors to generate

    Returns
    -------
    colors_list : list (Color)
        Hex codes, will need to use .hex() to retrieve as string

    '''
    start_color = Color(start_hex)
    end_color = Color(end_hex)
    
    colors_list = list(start_color.range_to(end_color, num))
    
    return colors_list

#------------------------------------------------------------------------------
def reformat_ticks(tick_val, pos):
    '''
    Create function to reformat large axis values

    returns function to set as formatter
    i.e. ax.xaxis.set_major_formatter(FuncFormatter(new_tick_format))
    '''
    if tick_val >= 1000:
        new_tick_format = round(tick_val/1000, 1)
    elif tick_val > -1000:
        new_tick_format = round(tick_val, 1)
    elif tick_val <= -1000:
        new_tick_format = round(tick_val/1000, 1)
    else:
        new_tick_format = tick_val

    new_tick_format = str(new_tick_format)
    
    index_of_decimal = new_tick_format.find(".")
    
    if index_of_decimal != -1:
        value_after_decimal = new_tick_format[index_of_decimal+1]
        if value_after_decimal == "0":
            # remove the 0 after the decimal point since it's not needed
            new_tick_format = \
                new_tick_format[0:index_of_decimal] + \
                new_tick_format[index_of_decimal+2:]

    return new_tick_format

#------------------------------------------------------------------------------
def labelprefix(limit):
    '''
    Set axis label unit prefix

    Parameters
    ----------
    limit : float
        Axis maximum

    Returns
    -------
    prefix : str
        Unit prefix

    '''
    prefix = ''
    if limit >= 1e9:
        prefix = 'G'
    elif limit >= 1e6:
        prefix = 'M'
    elif limit >= 1e3:
        prefix = 'k'
    return prefix

#------------------------------------------------------------------------------
def save_fig(plot, path, fn):
    '''
    Save figure as a .png file

    Parameters
    ----------
    plot : figure
        Name of plot
    path : str
        File directory path
    fn : str
        File name to save to

    Returns
    -------
    None

    '''
    save = path + fn + ".png"
    plot.savefig(save, bbox_inches="tight", pad_inches=0.2, dpi=1000)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' EIS DATA MANAGEMENT 
Functions in this section:
    - import_biologic_EIS
    - sort_EIS
    - sep_EIS_cycles
    - import_cycles
    - import_param_vals '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def import_biologic_EIS(path):
    '''
    Imports EIS data from biologic .mpr file

    Parameters
    ----------
    path : str
        File directory path

    Returns
    -------
    EIS_data : dataframe
        Full EIS dataset

    '''
    file = input("Data file: ")
    filepath = path + file + ".mpr"
        
    EIS_data = ecf.to_df(filepath)
        
    return EIS_data

#------------------------------------------------------------------------------
def sort_EIS(df):
    '''
    Sorts EIS data by time

    Parameters
    ----------
    df : dataframe
        Full EIS dataset

    Returns
    -------
    None

    '''
    df.sort_values(by=["time"])

#------------------------------------------------------------------------------
def sep_EIS_cycles(df, cycle_pts):
    '''
    Separates EIS data by cycle

    Parameters
    ----------
    df : dataframe
        Full EIS dataset
    cycle_pts : int
        Number of data points per EIS cycle

    Returns
    -------
    EIS_cycles : dict
        Dictionary of EIS data where each cycle is a separate entry 
    num_cycles : int
        Total number of EIS cycles

    '''
    num_rows = len(df)
    num_cycles = int(num_rows / cycle_pts)
    
    EIS_cycles = {}
    
    for i in range(0, num_cycles):
        EIS_cycles[i] = df.truncate(before=0 + (cycle_pts*i), after=(cycle_pts-1) +  (cycle_pts*i) )
    
    print("Number of cycles: " + str(num_cycles))
    
    return EIS_cycles, num_cycles

#------------------------------------------------------------------------------
def import_cycles(path, fn, num_cycles):
    '''
    Import EIS cycle data and fit from text file

    Parameters
    ----------
    path : str
        File directory path
    fn : str
        File name, see readme for formatting information
    num_cycles : int
        Total number of EIS cycles

    Returns
    -------
    expt_real : list (float)
        Real values from observed dataset
    expt_imag : list (float)
        Imaginary values from observed dataset
    fit_real : list (float)
        Real values from calculated dataset
    fit_imag : list (float)
        Imaginary values from calculated dataset

    '''
    file_list = []
    for i in range(num_cycles):
        cycle_fn = fn + "_cycle" + str(i+1) + ".txt"
        fit_fn = fn + "_cycle" + str(i+1) + "_fit.txt"
        file_list.append((cycle_fn, fit_fn))
        
    expt_real = []
    expt_imag = []
    fit_real = []
    fit_imag = []
    
    for i in range(num_cycles):
        expt_re, expt_im = np.loadtxt(path + file_list[i][0], unpack=True, dtype=float, delimiter=" ", skiprows=1)
        fit_re, fit_im = np.loadtxt(path + file_list[i][1], unpack=True, dtype=float, delimiter=" ", skiprows=1)
        expt_real.append(expt_re)
        expt_imag.append(expt_im)
        fit_real.append(fit_re)
        fit_imag.append(fit_im)
        
    return expt_real, expt_imag, fit_real, fit_imag

#------------------------------------------------------------------------------
def import_param_vals(path, fn, num_cycles):
    '''
    Import fit parameters from text file

    Parameters
    ----------
    path : str
        File directory path
    fn : str
        File name, see readme for formatting information
    num_cycles : int
        Total number of EIS cycles

    Returns
    -------
    param_vals : list (float)
        Fit parameter values

    '''
    file_list = []
    for i in range(num_cycles):
        file = fn + "_cycle" + str(i+1) + "_params.txt"
        file_list.append(file)
        
    param_vals = []
    for i in range(num_cycles):
        param_vals.append(np.loadtxt(path + file_list[i], unpack=True, dtype=float, delimiter=" "))
        
    return param_vals


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' EIS PLOTTING 
Functions in this section:
    - plot_singlecycle
    - plot_multicycle
    - plot_R
    - plot_fit_params
    - plot_sigC '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def plot_singlecycle(expt_re, expt_im, fit_re, fit_im, fit_color, x_lim):
    '''
    Generate plot with one EIS cycle

    Parameters
    ----------
    expt_re : list (float)
        Real values from observed dataset
    expt_im : list (float)
        Imaginary values from observed dataset
    fit_re : list (float)
        Real values from calculated dataset
    fit_im : list (float)
        Imaginary values from calculated dataset
    fit_color : str
        Hex code for calculated data color, format "#000000"
    x_lim : list (float)
        Tuple with x-axis minimum and maximum

    Returns
    -------
    None

    '''
    fig, (ax) = plt.subplots(1, figsize=(6,6))
    
    marker_style = dict(marker="o", markersize=8, markerfacecolor="white", markeredgecolor="black")
    
    # plot EIS data
    ax.plot(expt_re, expt_im, color="white", label="Experimental",  **marker_style)
    ax.plot(fit_re, fit_im, color=fit_color, linewidth=3, label="Calculated")
    
    # set axis limits
    ax.set_xlim(x_lim)
    ax.set_ylim(x_lim)
    
    # format tick values
    ax.xaxis.set_major_formatter(FuncFormatter(reformat_ticks))
    ax.yaxis.set_major_formatter(FuncFormatter(reformat_ticks))
    ax.tick_params(axis="both", labelsize="14")
    
    # set axis labels
    x_label = "Z$_{\mathrm{real}}$ (" + labelprefix(x_lim[1]) + "$\Omega$)"
    y_label = "$-$Z$_{\mathrm{imag}}$ (" + labelprefix(x_lim[1]) + "$\Omega$)"
    ax.set_xlabel(x_label, fontsize=16)
    ax.set_ylabel(y_label, fontsize=16)
    
    # add legend
    ax.legend(handlelength=1, fontsize="14")
    
    return(ax)

#------------------------------------------------------------------------------
def plot_multicycle(num_cycles, expt_re, expt_im, fit_re, fit_im, x_lim, start_hex=False, end_hex=False):
    '''
    Generate plot with multiple EIS cycles

    Parameters
    ----------
    num_cycles : int
        Total number of EIS cycles
    expt_re : list (float)
        Real values from observed dataset
    expt_im : list (float)
        Imaginary values from observed dataset
    fit_re : list (float)
        Real values from calculated dataset
    fit_im : list (float)
        Imaginary values from calculated dataset
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    start_hex : str, optional
        Hex code for initial gradient color, format "#000000". The default is False.
    end_hex : str, optional
        Hex code for final gradient color, format "#000000". The default is False.

    Returns
    -------
    None

    '''
    fig, (ax)=plt.subplots(1,figsize=(7,7))
        
    # generate color gradient
    if start_hex == False:
        start_hex = "#00C6BF"
    if end_hex == False:
        end_hex = "#B430C2"
    g = gradient_gen(start_hex, end_hex, num_cycles)
    
    marker_style = dict(marker="o", markersize=5, markerfacecolor="white", markeredgecolor="black")
    
    # plot EIS data
    for i in range(num_cycles):
        ax.plot(expt_re[i], expt_im[i], color="white", label="_Experimental", linewidth=1, **marker_style)
        ax.plot(fit_re[i], fit_im[i], color=g[i].hex, linewidth=1.5, label="_Calculated")
        
    # set axis limits
    ax.set_xlim(x_lim)
    ax.set_ylim(x_lim)
    
    # format tick values
    ax.xaxis.set_major_formatter(FuncFormatter(reformat_ticks))
    ax.yaxis.set_major_formatter(FuncFormatter(reformat_ticks))
    ax.tick_params(axis="both", labelsize="14")
    
    # set axis labels
    x_label = "Z$_{\mathrm{real}}$ (" + labelprefix(x_lim[1]) + "$\Omega$)"
    y_label = "$-$Z$_{\mathrm{imag}}$ (" + labelprefix(x_lim[1]) + "$\Omega$)"
    ax.set_xlabel(x_label, fontsize=16)
    ax.set_ylabel(y_label, fontsize=16)
    
    return(ax)

#------------------------------------------------------------------------------
''' plot R value vs cycle number '''
def plot_R(num_cycles, R_vals, x_lim, y_lim, color=False, marker=False):
    '''
    Plot R value vs cycle number

    Parameters
    ----------
    num_cycles : int
        Total number of EIS cycles
    R_vals : list (float)
        Resistance values from fit parameters
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    color : str, optional
        Hex code for calculated data color, format "#000000". The default is False.
    marker : str, optional
        Plot marker shape. The default is False.

    Returns
    -------
    None

    '''
    fig, (ax)=plt.subplots(1, 1, figsize=(8,4))
    
    if color == False:
        color = "#00C6BF"
    if marker == False:
        marker == "o"
    
    cycle_count = []
    for i in range(num_cycles):
        cycle_count.append(i)
    
    # plot data
    ax.plot(cycle_count, R_vals, color=color, marker=marker, markersize="8")
    
    # set axis limits
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    
    # format x axis
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    
    # set axis labels 
    y_label = "Resistance (" + labelprefix(y_lim[1]) + "$\Omega$)"
    ax.set_xlabel("Cycle Number", fontsize=16)
    ax.set_ylabel(y_label, fontsize=16)
    
    return(ax)

#------------------------------------------------------------------------------
def plot_fit_params(plot, vals, pos):
    '''
    Plot box with fit parameter info

    Parameters
    ----------
    plot : figure
        Name of plot
    vals : list (str)
        Formatted text of fit parameters
    pos : list (float)
        (x, y) position

    Returns
    -------
    None

    '''
    props = dict(facecolor="white", edgecolor="black", pad=8, linewidth=2)
    plot.text(pos[0], pos[1], vals, bbox=props, ha="left", va="top", fontsize="14", linespacing=2)
    
#------------------------------------------------------------------------------
def plot_sigC(plot, sig_C_vals, pos, box_color):
    '''
    Plot box with sigma and C values

    Parameters
    ----------
    plot : figure
        Name of plot
    sig_C_vals : list (str)
        Formatted text of sigma and C values
    pos : list (float)
        (x, y) position
    box_color : str
        Hex code for box border color, format "#000000"

    Returns
    -------
    None

    '''
    props = dict(facecolor="white", edgecolor=box_color, pad=8, linewidth=2)
    plot.text(pos[0], pos[1], sig_C_vals, ha="left", va="top", fontsize="14", linespacing=2, bbox=props)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' EIS PLOT FORMATTING
Functions in this section:
    - expt_info
    - singlecycle_info
    - labels_list
    - params_list
    - sig_cap_text '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def expt_info(curr_dens, num_cycles):
    '''
    Generates string with EIS experiment info

    Parameters
    ----------
    curr_dens : int
        Value of current density in uA/cm^2
    num_cycles : int
        Total number of EIS cycles

    Returns
    -------
    expt_info : str
        Description of current density and total number of cycles

    '''
    expt_info = "\n".join((r"$j=$" + str(curr_dens) + r" $\mu$A$\cdot$cm$^{-2}$", str(num_cycles) + " Cycles"))
    return expt_info

#------------------------------------------------------------------------------
def singlecycle_info(curr_dens, cycle, num_cycles):
    '''
    Generates string with EIS experiment info of a single cycle

    Parameters
    ----------
    curr_dens : int
        Value of current density in uA/cm^2
    cycle : int
        Current EIS cycle number
    num_cycles : int
        Total number of EIS cycles

    Returns
    -------
    expt_info : str
        Description of current density and current cycle

    '''
    expt_info = "\n".join((r"$j=$" + str(curr_dens) + r" $\mu$A$\cdot$cm$^{-2}$", "Cycle " + str(cycle) + " of " + str(num_cycles)))
    return expt_info

#------------------------------------------------------------------------------
def labels_list(num_list):
    '''
    Create list of parameter values in scientific notation

    Parameters
    ----------
    num_list : list (float)
        Values to be formatted

    Returns
    -------
    labels : list (str)
        Reformatted text values

    '''
    nums = []
    exp_list = []
    for i in range(len(num_list)):
        if num_list[i] < 1 and num_list[i] >= 0.1:
            nums.append(round(num_list[i],2))
            exp_list.append(0)
        else:
            ret_string = "{0:.{1:d}e}".format(num_list[i], 2)
            a, b = ret_string.split("e")
            b = int(b)
            nums.append(a)
            exp_list.append(b)
    
    formatted_exp = []
    for i in range(len(exp_list)):
        x = exp_list[i]
        if x == 2:
            formatted_exp.append(r"$^{2}$")
        elif x == 3: 
            formatted_exp.append(r"$^{3}$")
        elif x == 4: 
            formatted_exp.append(r"$^{4}$")
        elif x == 5: 
            formatted_exp.append(r"$^{5}$")
        elif x == -2: 
            formatted_exp.append(r"$^{-2}$")
        elif x == -3: 
            formatted_exp.append(r"$^{-3}$")
        elif x == -4: 
            formatted_exp.append(r"$^{-4}$")
        elif x == -5: 
            formatted_exp.append(r"$^{-5}$")
        elif x == -6: 
            formatted_exp.append(r"$^{-6}$")
        elif x == -7: 
            formatted_exp.append(r"$^{-7}$")
        elif x == -8: 
            formatted_exp.append(r"$^{-8}$")
        elif x == -9: 
            formatted_exp.append(r"$^{-9}$")
        elif x == -10: 
            formatted_exp.append(r"$^{-10}$")
        elif x == -11: 
            formatted_exp.append(r"$^{-11}$")
        elif x == -12: 
            formatted_exp.append(r"$^{-12}$")
        elif x == 0:
            formatted_exp.append("0")
    
    labels = []
    for i in range(len(nums)):
        if formatted_exp[i] != "0":
            label = nums[i] + " x 10" + formatted_exp[i]
            labels.append(str(label))
        else: 
            labels.append(str(nums[i]))

    return labels

#------------------------------------------------------------------------------
''' format text for displaying fit parameters '''
def params_list(vals):
    '''
    Create list of formatted parameter values to display as text

    Parameters
    ----------
    vals : list (str)
        Values to be formatted

    Returns
    -------
    params_txt : list (str)
        Formatted text values

    '''
    params_list = []
    
    for i in range(len(vals)):
        if vals[i][1] == "R":
            if vals[i][2] == 1:
                txt = r"R$_1$ = " + vals[i][0] + r" $\Omega$"
            elif vals[i][2] == 2:
                txt = r"R$_2$ = " + vals[i][0] + r" $\Omega$"
            #params_list.append(txt)
                
        if vals[i][1] == "n":
            #n = float(vals[i][0])
            #round_n = round(n, -len(str(int(n))) + 3)
            round_n = vals[i][0]
            if vals[i][2] == 1:
                txt = r"n$_1$ = {}".format(round_n)
            elif vals[i][2] == 2:
                txt = r"n$_1$ = {}".format(round_n)
            #params_list.append(txt)
                
        if vals[i][1] == "Q":
            if vals[i][2] == 1:
                txt = r"Q$_1$ = " + vals[i][0] + r" S$^n / \Omega$"
            elif vals[i][2] == 2:
                txt = r"Q$_2$ = " + vals[i][0] + r" S$^n / \Omega$"
        
        params_list.append(txt)
    
    params_txt = "\n".join((params_list))
            
    return params_txt

#------------------------------------------------------------------------------
''' format text for displaying sigma and C values '''
def sig_cap_text(sig_val, C_val):
    '''
    Create list of formatted sigma and C values to display as text

    Parameters
    ----------
    sig_val : float
        Sigma (conductivity) value
    C_val : float
        Capacitance value

    Returns
    -------
    sig_cap_txt : str
        Text with sigma and C text

    '''
    sig_txt = r"$\sigma = $ " + sig_val + " S/cm"
    C_txt = r"$C = $ " + C_val + " F"
    
    sig_cap_txt = "\n".join(sig_txt, C_txt)
    
    return sig_cap_txt