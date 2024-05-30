"""
PDF plotting functions
"""

import numpy as np
import glob
from colour import Color
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

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
''' PDF PLOTTING 
Functions in this section:
    - import_PDF
    - plot_PDF '''
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


def import_PDF(path, fn, header_rows):
    '''
    Imports PDF data file

    Parameters
    ----------
    path : str
        File directory path
    fn : str
        File name, including extension
    header_rows : int
        Number of rows of metadata at beginning of file

    Returns
    -------
    r : list (float)
        r (A) values
    G : list (float)
        G(r) values
    Gdiff : list (float)
        G_diff(r) values
    Gcalc : list (float)
        G_calc(r) values

    '''
    r, G, Gdiff, Gcalc = np.loadtxt(path+fn, unpack=True, dtype=float, skiprows=header_rows)
    return r, G, Gdiff, Gcalc


def plot_PDF(r, G, Gcalc, Gdiff=False, fit_color=False, x_lim=False, y_lim=False):
    '''
    Generates plot with PDF data

    Parameters
    ----------
    r : list (float)
        r (A) values
    G : list (float)
        G(r) values
    Gcalc : list (float)
        G_calc(r) values
    Gdiff : list (float), optional
        G_diff(r) values
    fit_color : str, optional
        Hex code for calculated data color, format "#000000". The default is False.
    x_lim : list (float), optional
        Tuple with x-axis minimum and maximum. The default is False.
    y_lim : list (float), optional
        Tuple with x-axis minimum and maximum.. The default is False.

    Returns
    -------
    None

    '''
    plt.figure(figsize=(7,7))
    
    if fit_color == False:
        fit_color = "#00C6BF"
    diff_pos = np.min(G) - 0.1
    
    # plot observed G(r), calculated G(r), and difference curve
    for i in range(len(r)):
        plt.plot(r, G, color="black", label="_Observed", marker=".")
        plt.plot(r, Gcalc, color=fit_color, label="_Calculated", linewidth="2")
        if Gdiff != False:
            plt.plot(r, Gdiff+diff_pos, color="#BEBEBE", label="_Difference")
    
    # set axis limits
    if x_lim != False:
        plt.xlim(x_lim)
    if y_lim != False:
        plt.ylim(y_lim)
    
    # set axis labels
    plt.xlabel("r / " r"$\AA$")
    plt.ylabel("G(r) / " r"$\AA^{-2}$")
    
    # create legend
    G_handle = mlines.Line2D([], [], color="white", label="Observed", 
                               marker=".", mfc="black", ms=15)
    Gcalc_handle = mlines.Line2D([], [], color=fit_color, label="Calculated",
                                 linewidth="2")
    Gdiff_handle = mlines.Line2D([], [], color="#BEBEBE", label="Difference",
                                 linewidth="2")
    if Gdiff != False:
        plt.legend(handles=[G_handle, Gcalc_handle, Gdiff_handle],
                   handlelength=1, fontsize="14")
    elif Gdiff == False:
        plt.legend(handles=[G_handle, Gcalc_handle],
                   handlelength=1, fontsize="14")
    
    return(plt)