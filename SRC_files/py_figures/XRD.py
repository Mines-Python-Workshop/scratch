"""
Diffraction plotting functions
"""

import numpy as np
from colour import Color
import glob, os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

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
''' converts from 2theta to Q'''
def tt_to_q(twotheta, wavelength):
    '''
    Converts 2theta to Q

    Parameters
    ----------
    twotheta : list (float)
        2theta values
    wavelength : float
        Instrument wavelength

    Returns
    -------
    Q : list (float)
        Corresponding Q values

    '''
    Q = 4 * np.pi * np.sin((twotheta * np.pi)/360) / wavelength
    return Q

#------------------------------------------------------------------------------
def q_to_tt(q, wavelength):
    '''
    Converts from Q to 2theta

    Parameters
    ----------
    q : list (float)
        Q values
    wavelength : float
        Instrument wavelength

    Returns
    -------
    twotheta : list (float)
        Corresponding 2theta values

    '''
    twotheta = 360 * np.pi * np.arcsin((q * wavelength) / (4 * np.pi))
    return twotheta

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
def gradient_gen_2D(top_left, top_right, bottom_left, bottom_right, rows, cols):
    '''
    Generates a 2D color gradient

    Parameters
    ----------
    top_left : str
        Hex code for top left gradient color, format "#000000"
    top_right : str
        Hex code for top right gradient color, format "#000000"
    bottom_left : str
        Hex code for bottom left gradient color, format "#000000"
    bottom_right : str
        Hex code for bottom right gradient color, format "#000000"
    rows : int
        Number of rows
    cols : int
        Number of columns

    Returns
    -------
    color_list : list (Color)
        Hex codes in 2D array format, will need to use .hex() to retrieve as string

    '''
    tl = Color(top_left)
    tr = Color(top_right)
    bl = Color(bottom_left)
    br = Color(bottom_right)
    
    start_col_colors = list(tl.range_to(bl, rows))
    end_col_colors = list(tr.range_to(br, rows))
    
    color_list = []
    for i in range(0, rows):
        c1 = start_col_colors[i]
        c2 = end_col_colors[i]
        curr_row_colors = list(c1.range_to(c2, cols))
        color_list.append(curr_row_colors)
    
    return color_list

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
def labelexp(limit):
    '''
    Sets axis label unit exponent for scaled axis

    Parameters
    ----------
    limit : float
        Axis maximum

    Returns
    -------
    exp : str
        Unit exponent

    '''
    exp = ''
    if limit >= 1e9:
        exp = "x $10^{9}$"
    elif limit >= 1e6:
        exp = "x $10^{6}$"
    elif limit >= 1e3:
        exp = "x $10^{3}$"
    return exp

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
''' XRD DATA MANAGEMENT 
Functions in this section:
    - var_dicts 
    - diff_curve
    - export_diff 
    - norm '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
''' create XRD variable dictionaries '''
def var_dicts(current_dir, header_rows, x_vals, y_vals, isQ=True):
    '''
    Create variable dictionaries for Q/2theta and intensity

    Parameters
    ----------
    current_dir : list (str)
        File names in working directory
    header_rows : int
        Number of rows of metadata at beginning of file
    x_vals : list (float)
        Q or 2theta values
    y_vals : list (float)
        Intensity values
    isQ : bool, optional
        Set to False if using units of 2theta. The default is True.

    Returns
    -------
    x_dict : dict (float)
        X-axis values sorted by corresponding file index
    y_dict : dict (float)
        Y-axis values sorted by corresponding file index

    '''

    x_dict = {}
    y_dict = {}
    
    # create keys for dictionaries with index number
    for i in current_dir:
        if isQ == True:
            globals()[f"Q_{i}"] = f"Q_{i}"
            x_dict[i] = "Q"+str(i)
        elif isQ == False:
            globals()[f"2theta_{i}"] = f"2theta_{i}"
            x_dict[i] = "2theta"+str(i)
        
        globals()[f"ints_{i}"] = f"ints_{i}"
        y_dict[i] = "ints"+str(i)
            
    for i in range(0, len(x_dict)):
        x_dict[i], y_dict[i] = np.loadtxt(current_dir[i], unpack=True, dtype=float, skiprows=header_rows)
                
    return x_dict, y_dict

#------------------------------------------------------------------------------
def diff_curve(x1, x2, y1, y2):
    '''
    Calculates the difference between two datasets

    Parameters
    ----------
    x1 : list (float)
        x-axis values of dataset 1
    x2 : list (float)
        x-axis values of dataset 2
    y1 : list (float)
        y-axis values of dataset 1
    y2 : list (float)
        y-axis values of dataset 2

    Returns
    -------
    diff_x : list (float)
        x-axis values of difference curve
    diff_y : list (float)
        y-axis values of difference curve

    '''
    diff_x_list = []
    diff_y_list = []
    for i in range(len(x1)):
        x1_val = float('%.3f'%(x1[i]))
        for j in range(len(x2)):
            x2_val = float('%.3f'%(x2[j]))
            if x1_val == x2_val:
                diff_x_list.append(x1_val)
                diff_y_list.append(y1[i]-y2[j])
    diff_x = np.array(diff_x_list)
    diff_y = np.array(diff_y_list)
    
    return diff_x, diff_y

#------------------------------------------------------------------------------
def export_diff(x, y, path, fn):
    '''
    Export difference curve to text file

    Parameters
    ----------
    x : list (float)
        x-axis values of difference curve
    y : list (float)
        y-axis values of difference curve
    path : str
        File directory path
    fn : str
        File name

    Returns
    -------
    None

    '''
    save = os.path.join(path, fn + ".txt")
    with open(save, "w") as f:
        for (x, y) in zip(x, y):
            f.write("{0} {1}\n".format(x, y))
    f.close()

#------------------------------------------------------------------------------
def norm(ints):
    '''
    Normalize intensity values

    Parameters
    ----------
    ints : list (float)
        Original intensity values

    Returns
    -------
    norm_ints : list (float)
        Normalized intensity values

    '''
    norm_ints = (ints / np.max(ints))
    return norm_ints

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' XRD PLOTTING 
Functions in this section:
    - single_fit
    - add_hkl
    - hkl_diff_subplots 
    - stacked_single_plot
    - stacked_subplots '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def single_fit(x, obs, calc, diff, x_lim, y_lim, isQ=True, Q_wl=None, 
               isNorm=False, fit_color=False):
    '''
    Generates plot with a single diffraction refinement

    Parameters
    ----------
    x : list (float)
        X-axis data (either Q or 2theta)
    obs : list (float)
        Observed intensities
    calc : list (float)
        Calculated intensities
    diff : list (float)
        Difference of observed and calculated intensities
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    isQ : bool, optional
        Set to False if using units of 2theta. The default is True.
    Q_wl : float, optional
        Set instrument wavelength. The default is None.
    isNorm : bool, optional
        Set to True if intensity data is normalized. The default is False.
    fit_color : str, optional
        Hex code for calculated data color, format "#000000". The default is False.

    Returns
    -------
    None

    '''
    
    fig, (ax)=plt.subplots(1,figsize=(7,7))
    
    
    if fit_color == False:
        fit_color = "#00B8FF"
    
    # plot XRD data
    ax.scatter(x, obs, color="black", label="Observed", marker=".", s=8)
    ax.plot(x, calc, color=fit_color, label="Calculated", linewidth=2)
    ax.plot(x, diff-np.max(diff), color="#BEBEBE", label="Difference", linewidth=1)

    # set axis limits
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)

    # format axes
    ax.yaxis.set_major_formatter(FuncFormatter(reformat_ticks))
    ax.tick_params(axis="both", labelsize="14")
    
    # set axis labels
    if isQ == False:
        x_label = r"2$\theta$ / $^{circ}$ (Cu K$\alpha$)"
    elif isQ == True:
        x_label = "Q (\AA" r"$^{-1}$, $\lambda=$" + str(Q_wl) + " \AA)"
    
    y_exp = labelexp(y_lim[1])
    if isNorm == False:
        y_label = "Intensity (counts " + y_exp + ")"
    elif isNorm== True:
        y_label = "Intensity (counts, normalized)"

    ax.set_xlabel(x_label, fontsize=16)
    ax.set_ylabel(y_label, fontsize=16)
    
    # add legend
    ax.legend(handlelength=1, fontsize="14")
    
    return(ax)

#------------------------------------------------------------------------------
def add_hkl(plot, xvals, tick_height, tick_width, tick_pos, color):
    '''
    Adds (hkl) reflection tick marks to plot

    Parameters
    ----------
    plot : figure
        Name of plot
    xvals : list (float)
        X-axis positions of (hkl) reflections
    tick_height : float
        Height in axis units
    tick_width : float
        Width in axis units
    tick_pos : float
        Position of bottom of tick marks in axis units
    color : str
        Hex code for calculated data color, format "#000000"

    Returns
    -------
    None

    '''
    plot.bar(xvals, height=tick_height, width=tick_width, bottom=tick_pos, align="center", 
             data=None, color=color)

#------------------------------------------------------------------------------
def hkl_diff_subplots(x, obs, calc, diff, hkl_vals, x_lim, data_y_lim, diff_y_lim, 
                      isQ=True, Q_wl=None, isNorm=False, fit_color=False, hkl_color=False):
    '''
    Generates XRD plot with subplots for (hkl) ticks and difference curb

    Parameters
    ----------
    x : list (float)
        X-axis data (either Q or 2theta)
    obs : list (float)
        Observed intensities
    calc : list (float)
        Calculated intensities
    diff : list (float)
        Difference of observed and calculated intensities
    hkl_vals : list (float)
        X-axis positions of (hkl) reflections
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    data_y_lim : list (float)
        Tuple with y-axis minimum and maximum for XRD data plot
    diff_y_lim : list (float)
        Tuple with y-axis minimum and maximum for difference plot
    isQ : bool, optional
        Set to False if using units of 2theta. The default is True.
    Q_wl : float, optional
        Set instrument wavelength. The default is None.
    isNorm : bool, optional
        Set to True if intensity data is normalized. The default is False.
    fit_color : str, optional
        Hex code for calculated data color, format "#000000". The default is False.
    hkl_color : str, optional
        Hex code for (hkl) color, format "#000000". The default is False.

    Returns
    -------
    None

    '''
    
    fig, (ax) = plt.subplots(3, 1, figsize=(9,6), gridspec_kw={'height_ratios': [4, 0.5, 1.25]})
    
    if fit_color == False:
        fit_color = "#00B8FF"
    if hkl_color == False:
        hkl_color = "#97DB4F"
    
    # plot XRD data
    obs = ax[0].scatter(x, obs, color="black", label="Observed", marker=".", s=8)
    calc = ax[0].plot(x, calc, color=fit_color, label="Calculated", linewidth=2)
    
    # plot hkl
    for i in range(len(hkl_vals)):
        x_range = np.array([hkl_vals[i], hkl_vals[i]])
        y_range = np.array([0.2, 0.8])
        ax[1].plot(x_range, y_range, color=hkl_color)
    
    # plot difference
    diff = ax[2].plot(x, diff, color="#BEBEBE", label="Difference", linewidth=1)

    # set axis limits
    for i in range(3):
        ax[i].set_xlim(x_lim)
        ax[i].yaxis.set_major_formatter(FuncFormatter(reformat_ticks))
        ax[i].tick_params(axis="both", labelsize="14")
        
    ax[0].set_ylim(data_y_lim)
    ax[1].set_ylim(0, 1)
    ax[2].set_ylim(diff_y_lim)
    
    ax[0].get_xaxis().set_visible(False)
    ax[1].get_xaxis().set_visible(False)
    ax[1].get_yaxis().set_visible(False)
    
    # set axis labels
    if isQ == False:
        x_label = r"2$\theta$ / $^{circ}$ (Cu K$\alpha$)"
    elif isQ == True:
        x_label = "Q (\AA" r"$^{-1}$, $\lambda=$" + str(Q_wl) + " \AA)"
    
    y_exp = labelexp(data_y_lim[1])
    if isNorm == False:
        y_label = "Intensity (counts " + y_exp + ")"
    elif isNorm== True:
        y_label = "Intensity (counts, normalized)"

    ax[2].set_xlabel(x_label, fontsize=16)
    ax[0].set_ylabel(y_label, fontsize=16)
    
    # add legend
    ax[0].legend(handlelength=1, fontsize="14")
    ax[2].legend(handlelength=1, fontsize="14")
    
    plt.subplots_adjust(hspace=0.05)
    
    return(ax)

#------------------------------------------------------------------------------
def stacked_single_plot(x_lim, y_lim, num, x_vals, y_vals, spacing, ycalc_vals=None, 
                        labels=None, label_offsets=None, start_hex=False, end_hex=False, 
                        isQ=True, Q_wl=None, isNorm=False):
    '''
    Generates plot with multiple XRD datasets on a single plot

    Parameters
    ----------
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    num : int
        Total number of datasets
    x_vals : list (float)
        X-axis data (either Q or 2theta)
    y_vals : list (float)
        Observed intensities
    spacing : float
        Vertical spacing between datasets on plot
    ycalc_vals : list (float), optional
        Calculated intensities. The default is None.
    labels : list (str), optional
        Labels for each dataset. The default is None.
    label_offsets : list (float), optional
        Tuple with offsets from x-axis maximum and vertical spacing from data
        for text labels. The default is None.
    start_hex : str, optional
        Hex code for initial gradient color, format "#000000". The default is False.
    end_hex : str, optional
        Hex code for final gradient color, format "#000000". The default is False.
    isQ : bool, optional
        Set to False if using units of 2theta. The default is True.
    Q_wl : float, optional
        Set instrument wavelength. The default is None.
    isNorm : bool, optional
        Set to True if intensity data is normalized. The default is False.

    Returns
    -------
    None

    '''
    
    fig, (ax) = plt.subplots(1,figsize=(7,7))
    
    # generate color gradient
    if start_hex == False:
        start_hex = "#00C6BF"
    if end_hex == False:
        end_hex = "#B430C2"
    g = gradient_gen(start_hex, end_hex, num)
    
    # plot data
    if ycalc_vals is None:
        for i in range(num):
            ax.plot(x_vals[i], y_vals[i] + (i * spacing), color=g[i].hex, linewidth="2")
    elif ycalc_vals is not None:
        for i in range(num):
            ax.scatter(x_vals[i], y_vals[i] + (i * spacing), color="black", label="Observed", marker=".", s=8)
            ax.plot(x_vals[i], ycalc_vals[i] + (i * spacing), color=g[i].hex, linewidth="2")
            
    # set axis limits
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    
    # format axes
    ax.yaxis.set_major_formatter(FuncFormatter(reformat_ticks))
    ax.tick_params(axis="both", labelsize="14")
    
    # set axis labels
    if isQ == False:
        x_label = r"2$\theta$ / $^{circ}$ (Cu K$\alpha$)"
    elif isQ == True:
        x_label = "Q (\AA" r"$^{-1}$, $\lambda=$" + str(Q_wl) + " \AA)"
    
    y_exp = labelexp(y_lim[1])
    if isNorm == False:
        y_label = "Intensity (counts " + y_exp + ")"
    elif isNorm== True:
        y_label = "Intensity (counts, normalized)"

    ax.set_xlabel(x_label, fontsize=16)
    ax.set_ylabel(y_label, fontsize=16)
    
    # add stack labels
    if labels is not None:
        for i in range(num):
            ax.text(x_lim[1] - label_offsets[0], label_offsets[1] + (i * spacing),
                    labels[i], color=g[i].hex, fontsize="16", ha="right", va="top")
    
    return(ax) 

#------------------------------------------------------------------------------
def stacked_subplots(x_lim, y_lim, num, x_vals, y_vals, ycalc_vals=None, diff=None,
                        labels=None, label_offsets=None, start_hex=False, 
                        end_hex=False, isQ=True, Q_wl=None, isNorm=False):
    '''
    Generates plot with a subplot for each XRD dataset

    Parameters
    ----------
    x_lim : list (float)
        Tuple with x-axis minimum and maximum
    y_lim : list (float)
        Tuple with y-axis minimum and maximum
    num : int
        Total number of datasets
    x_vals : list (float)
        X-axis data (either Q or 2theta)
    y_vals : list (float)
        Observed intensities
    ycalc_vals : list (float), optional
        Calculated intensities. The default is None.
    diff : list (float), optional
        Difference curve intensities. The default is None.
    labels : list (str), optional
        Labels for each dataset. The default is None.
    label_offsets : list (float), optional
        Tuple with offsets from x-axis and y-axis maximas for text labels. 
        The default is None.
    start_hex : str, optional
        Hex code for initial gradient color, format "#000000". The default is False.
    end_hex : str, optional
        Hex code for final gradient color, format "#000000". The default is False.
    isQ : bool, optional
        Set to False if using units of 2theta. The default is True.
    Q_wl : float, optional
        Set instrument wavelength. The default is None.
    isNorm : bool, optional
        Set to True if intensity data is normalized. The default is False.

    Returns
    -------
    None

    '''
    
    if diff is None:
        fig, (ax) = plt.subplots(num, 1, figsize=(8,8))
    elif diff is not None:
        fig, (ax) = plt.subplots(num, 1, figsize=(num*4,8))
    
    # generate color gradient
    if start_hex == False:
        start_hex = "#00C6BF"
    if end_hex == False:
        end_hex = "#B430C2"
    g = gradient_gen(start_hex, end_hex, num)
    
    # plot data
    if ycalc_vals is None:
        for i in range(num):
            ax[i].plot(x_vals[i], y_vals[i], color=g[i].hex, linewidth="2")
    elif ycalc_vals is not None:
        for i in range(num):
            ax[i].scatter(x_vals[i], y_vals[i], color="black", label="Observed", marker=".", s=8)
            ax[i].plot(x_vals[i], ycalc_vals[i], color=g[i].hex, linewidth="2")
            if diff is not None:
                ax[i].plot(x_vals[i], diff[i]-np.max(diff), color="#BEBEBE", linewidth="1")

    # set axis limits
    for i in range(num):
        ax[i].set_xlim(x_lim)
        ax[i].set_ylim(y_lim)
        ax[i].yaxis.set_major_formatter(FuncFormatter(reformat_ticks))
        ax[i].tick_params(axis="both", labelsize="14")
        if i < (num-1):
            ax[i].get_xaxis().set_visible(False)

    # set axis labels
    if isQ == False:
        x_label = r"2$\theta$ / $^{circ}$ (Cu K$\alpha$)"
    elif isQ == True:
        x_label = "Q (\AA" r"$^{-1}$, $\lambda=$" + str(Q_wl) + " \AA)"
    
    y_exp = labelexp(y_lim[1])
    if isNorm == False:
        y_label = "Intensity (counts " + y_exp + ")"
    elif isNorm== True:
        y_label = "Intensity (counts, normalized)"

    ax[-1].set_xlabel(x_label, fontsize=16)
    fig.supylabel(y_label, fontsize=16)
    
    # add stack labels
    if labels is not None:
        for i in range(num):
            ax[i].text(x_lim[1] - label_offsets[0], label_offsets[1],
                    labels[i], color=g[i].hex, fontsize="16", ha="right", va="top")
            
    plt.subplots_adjust(hspace=0.05)
    
    return(ax)
        
        
        
        
        
        
        
