import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly 
import seaborn as sns
from scipy.optimize import curve_fit
import pdb

### Import files, and normalizes data from columns. 
def import_file(fileName, normalize=True):
    data = pd.read_csv(fileName, index_col=1)
    if 'Unnamed: 0' in data.columns.to_list():
        data.drop(columns='Unnamed: 0', inplace=True)  
    if normalize:
        data = (data-np.min(data)) / (np.max(data) - np.min(data)) #normalize 
    return data 

def separate_standard_curve(data, ref_cols):
    #standard curve data points 
    refs = data[ref_cols]
    data.drop(columns=ref_cols, inplace=True)
    return refs, data 

### Input: Beginning and ending column numbers from 96 well qPCR plate, row letters desired
### Output: labels from beg values to end w/ letters appended. 
### ie. ['A10', 'A11', 'A12', 'B10', 'B11', 'B12']
def make_labels(beg=1, end=3, rows=['A','B','C','D','E','F','G','H'], verbose=False):
    nums = range(beg, end+1)
    labels = []
    for row in rows:
        for num in nums:
            labels.append(row+str(num))
    if verbose: print("Labels Made: ", labels)
    return labels  

def get_Cqs(data, threshold=0.5, verbose=True):
    Cq = []
    #columns = list(data)
    for col in data.columns:
        #pdb.set_trace()
        res = list(map(lambda i:i>threshold, data[col])).index(True)
        if res:
            y1 = data[col].loc[res]
            y2 = data[col].loc[res+1]
        else:
            if verbose: print("ERROR in column: ", col, " Setting Cq as np.nan")
            res = np.nan
        Cq.append( res + (threshold - y1) / (y2-y1)) 
   
        #res + (.5  - data_norm['A1'].loc[res])/(data_norm['A1'].loc[res+1] - data_norm['A1'].loc[res] )
    return Cq  

def analyze_Cqs(data, threshold=0.25, skip=4, rotate = False, verbose=True):
    #TODO:  Add triplicates value? 
    Cqs = get_Cqs(data, threshold, verbose)
    if rotate:
        Cq_means = [np.mean(Cqs[val::skip]) for val in np.arange(0,skip)]
        Cq_std = [np.std(Cqs[val::skip]) for val in np.arange(0,skip)] 
    else:
        Cq_means = [ np.mean( Cqs[val: val+3]) for val in np.arange(0,len(Cqs), 3)]
        Cq_std = [ np.std( Cqs[val: val+3]) for val in np.arange(0,len(Cqs), 3)] 
    return Cqs, Cq_means, Cq_std

def plot_bar(means, stds, labels, title="Cq means", ylim=0):
    
    x_pos = np.arange(len(means))
    fig, ax = plt.subplots()
    ax.bar(x_pos, means, yerr=stds, align='center', alpha=.5)
    ax.set_ylabel('Cq')
    ax.set_xticks(x_pos)#, labels, rotation='vertical')
    ax.set_xticklabels(labels, rotation=45, ha="right")
    if ylim:
        plt.ylim(ylim,32)    
   # plt.tight_layout()
    for i, v in enumerate(means):
        plt.text(i-.15, v-2, str(round(v,2)), rotation=90, color='blue', fontweight='bold')
    plt.title(title)
    
def fitplot(data):
    data_no0 = data.drop(index=[21,22,23])
    coefs = poly.polyfit(np.log10(data_no0.Conc),  data_no0.Cq,1) 
    ffit = poly.polyval(np.log10(data_no0.Conc), coefs)
    Cq_0_mean = np.mean(data.Cq[21:24])  
    Cq_0_std = np.std(data.Cq[21:24])
    LOD = 10**((Cq_0_mean - 3*Cq_0_std - coefs[0])/coefs[1])    

    return ffit   

def cqToDNA(Cq, coefs):
    m = coefs[1]
    b = coefs[0] 
    dna = 10**((Cq-b)/m)
    return dna 


def plot_conc_highlights(plt_horiz=True, 
                         targs = ["CRP", "GDF-15", "IL-1ra", "IL-6"],
                        large=False):
    
#     LUMINEX/presented on paper
    full_targ_range = {"CRP": [2000, 100000, "b"], # CRP
                  "GDF-15": [5, 120, "orange"], #GDF-15
                  "IL-1ra": [6, 600, "g"], #IL-1ra
                  "IL-6": [.002, .300, "purple"], #IL-6
                 }
#From Dan 8/6/2020 on slack 
    full_targ_range_large = {"CRP": [861/20, 626000/20, "b"], # CRP
              "GDF-15": [4/10.0, 3500/10.0, "orange"], #GDF-15
              "IL-1ra": [5/4.0, 1300/4.0, "g"], #IL-1ra
              "IL-6": [.064/2.0, 260/2.0, "purple"], #IL-6
             }
    if large:
        full_targ_range = full_targ_range_large
    targ_range = {key: full_targ_range[key] for key in targs}
    for target in targ_range:
        minT, maxT, col = targ_range[target]
        plt.axvspan(minT, maxT , alpha=0.05, color=col) 
    if plt_horiz: plt.axhspan(10**-15.5, 10**-14.5, alpha=0.05, color="gray")

def plot_mid_targ_lines(mid_target_conc = {"CRP":811, "GDF15":4.741, 
                                           "IL1ra":2.915,
                                           "IL6":.015}, #[811, 4.741, 2.915, .064, 0.015], 
                        targs = ["CRP", "GDF15", "IL1ra", "IL6"],
                        colors = ["b", "orange", "g", "purple"],
                        ymin = 10**-19,
                        ymax = 10**-14,
                        plt_hline = True, hline_val = 10**-15,
                        plt_point_only = False,
                        marker = "x", 
                        line_style = "dashed",
                        ):
    

    if isinstance(targs, str): targs = [targs] #Ensure targets are a list even if only 1 str
    targs = [targ.replace("-", "") for targ in targs]    
        
    for i, target in enumerate(targs):
        #plot target vertical lines for every target
        if not plt_point_only: 
            plt.vlines(mid_target_conc[target], ymin, ymax, colors=colors[i], 
                       alpha=.5, linestyles=line_style)
        else:
            if len(targs) !=1: 
                plt.plot(mid_target_conc[target], hline_val, color='r', marker=marker, label="Desired Intersection" )
            else: 
                plt.plot(mid_target_conc[target], hline_val, color='r', marker=marker , label="Desired Intersection ")
           
    if plt_hline:         
        #plot horizontal line -- typically at the output DNA concentration 
        conc_range = list( map(mid_target_conc.get, targs) )
        plt.hlines(hline_val, min(conc_range)/10, max(conc_range)*10,
                   linestyles= line_style, alpha=.5)

def plot_amplification(data):
    plt.plot(data)
    plt.xlabel("Cycles")
    plt.ylabel("Intensity")
    plt.title("amplification")
    
def run_data_to_cqs_on_plate(topdir ="", data_file="", plot_data=True, cq_thresh=0.25, verbose=True):
    """
    Input: 
    Topdir: top directory like "experiments/2020_07_22_calibration/"
    data_file: qPCR Quant Amp file like "admin_2020-07-22 19-12-18_CT022252-PLA-calibration -  Quantification Amplification Results_SYBR.csv"
    plate_file: "admin_2020-07-22 19-12-18_CT022252-PLA-calibration -  End Point Results_SYBR.csv"
    Output: 
        data: original amplification results
        plate_df_cal: Summary of each well w/ Cqs 
    
    """
    #import data
    data = import_file(topdir + data_file 
                            + " -  Quantification Amplification Results_SYBR.csv")
    
    if plot_data:
        plot_amplification(data)

    #import plate data
    plate_file = topdir + data_file + " -  End Point Results_SYBR.csv"
    plate_df_cal = pd.read_csv(plate_file) #, index_col=1)

    # Get Cqs and means 
    Cqs, Cq_means, Cq_std = analyze_Cqs(data, threshold=cq_thresh, verbose=verbose)
    plate_df_cal.loc[:, 'Cqs'] = Cqs
    
    return data, plate_df_cal 


def calc_Cq_to_DNA(cal_file="1B_ctrl_output.csv"):
    # Fit Cq vs DNA output from standard curve
    refs = pd.read_csv(cal_file)
    refs_no0 = refs[refs.Target != 0]
    ref_concs_fM = refs_no0.Target.to_numpy()*1000
    log_concs = np.log10(ref_concs_fM)
    coefs = poly.polyfit(log_concs, refs_no0.Cqs.tolist(),1) 
    return coefs

def convert_df_cq_to_DNA(plate_df, coefs):
    plate_df.loc[:, 'dna_conc'] = cqToDNA(plate_df.Cqs, coefs)
    plate_df.loc[:, "dna_conc_M"] = plate_df.dna_conc * 10**-15
    return plate_df 

def plot_dna_conc(plate_df, coefs, targetname="Target", colorsplit="Content", 
                  yaxis="dna_conc_M", ylogscale=True, xlogscale=True, 
                  style="Content", figsize = (8,4)):
#     plate_df.loc[:, 'dna_conc'] = cqToDNA(plate_df.Cqs, coefs)
#     plate_df.loc[:, "dna_conc_M"] = plate_df.dna_conc * 10**-15
    plate_df = convert_df_cq_to_DNA(plate_df, coefs)
    plate_df.loc[:, "Targetp1"] = plate_df.loc[:, targetname] + 0.001 # So that we can plot 0s!! 
    plt.figure(figsize=figsize)
#     cmap = sns.cubehelix_palette(dark=.3, light=.8, as_cmap=True)
    sns.scatterplot(x='Targetp1', y=yaxis, 
                    data=plate_df, 
                    style =style,
                    hue=colorsplit, 
                   
                   )
    if xlogscale: 
        plt.xscale("log")
        plt.xlabel("log(Target) (pM)")
    if ylogscale: plt.yscale("log")

#     plt.ylabel("[DNA] M")
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    return plate_df
    
    
def plot_scatter(df, x='Target', y="dna_conc_M",hue="Probe3'", title="", 
                 xlog=True, plt_mid=True, targs = [], figsize= (4,4), 
                 style ="Content", alpha=.5, linthreshx=.01):
    """
    DEPRECATED SCATTER PLOTTING. For updated plotting w. line plots see plot_line_w_no_targ
    Input: Pandas dataframe with columns of Target concentration, "Content" for target type, and calculated dna_conc_M (in molarity). 
    Output: Plots target concentration vs DNA_conc_M in a scatter plot.
    """
#     df.loc[df.Target == 0, "Target"] = .001
    plt.figure(figsize= figsize)
    palette = ["r", "orange", "g", "c", "b", "purple", "k"] 
    palette = palette[0:len(df[hue].unique())]
    
#     sns.lineplot(x=x, y=y, data=df, alpha=alpha,
#                     hue=hue, palette = palette, style =style)
    sns.scatterplot(x=x, y=y, data=df, alpha=alpha,
                hue=hue, palette = palette, style =style)
    if xlog: 
#         linthreshx = df.loc[df.Target != 0, "Target"].min()
        plt.xscale('symlog', linthresh=linthreshx)#"symlog",  nonposx='clip')
    plt.yscale("log")
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    if plt_mid:
        plot_mid_targ_lines(targs=targs, ymin= min(df.dna_conc_M), 
                            ymax = max(df.dna_conc_M), colors=palette)
    plt.title(title)

def plot_line_w_no_targ(plate_df, x='Target', y="dna_conc_M",
                        hue="Content", title="", save_fig_name="", 
                        ylim=[], targets=['CRP' 'GDF15' 'IL1ra' 'IL1B' 'IL6'], 
                        plt_mid=True, style =None, legend=False
                       ):
#                      xlog=True, plt_mid=True, targs = [], figsize= (4,4), 
#                  style ="Content", alpha=.5, linthreshx=.01):
    """
    Input: Pandas dataframe with columns of Target concentration, "Content" for target type, and calculated dna_conc_M (in molarity). 
    Output: Plots target concentration vs DNA_conc_M in a line plot w/ std's. Also includes scatter plot of any no target controls (ie Target == 0). 
    """
    order = plate_df.Content.unique()
    palette_orig = ["b", "orange", "g", "c", "r", "purple", "k"] 
    palette = palette_orig[0:len(plate_df[plate_df.Target > 0][hue].unique())]

    sns.lineplot(data=plate_df[plate_df.Target > 0], x=x, y=y,
                 hue=hue, err_style="bars", palette=palette, legend=legend, style=style)
    sns.scatterplot(data=plate_df[plate_df.Target == 0], x=x, y=y,
                 hue=hue, palette=palette_orig, alpha=0.5, hue_order=order)

    plt.xscale('symlog', linthreshx=.01)
    plt.yscale("log")
    plt.ylabel("[DNA] M")
    plt.xlabel("[Target] pM")
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    if plt_mid:
        plot_mid_targ_lines(targs=targets,plt_point_only = False, 
                            plt_hline = True, colors= palette)
    plt.title(title)

    if ylim!=[]: plt.ylim(ylim)
    plt.savefig(save_fig_name, dpi=300, bbox_inches='tight')
    
    
### PROBE FITTING  Specific Functions  ### 
def get_DNA(log_targ_conc, a, b, c, *probe_conc): 
    #log_targ_conc: in pM
#     probe_conc
#     pdb.set_trace()
    log_dna_conc = a*log_targ_conc + b*np.log10(probe_conc) + c
    return log_dna_conc

def get_probe(log_dna_conc, log_targ_conc, a,b,c):
    """
    INPUT
    log_dna_conc: in M 
    log_targ_conc: in pM
    OUTPUT:
    probe_conc: probe concentration in pM
    """
    log_probe_conc = (log_dna_conc - a*log_targ_conc  - c)/ b
    probe_conc = 10**log_probe_conc
    return probe_conc

def convert_target_DNA_to_log(df):
    df.loc[:, 'log_target'] = np.log10(df.Targetp1)
    df.loc[:, 'log_dna'] = np.log10(df.dna_conc_M) 
    return df

def load_run_files(topdir, datafile):
    """
    :param topdir: Top directory that holds experimenttal file:
        ex: "../experiments/2022_02_16_PLA_GenVortex/"
    :param datafile: File prefix from qPCR data:
        ex: "admin_2022-02-16 14-25-49_CT022252_PLAgen_vortex"
    :return:
    """
    data, plate_df = run_data_to_cqs_on_plate(topdir, data_file, plot_data=True)
    coefs_1b = calc_Cq_to_DNA(cal_file=mydrive + "1B_ctrl_output.csv")
    plate_df = convert_df_cq_to_DNA(plate_df, coefs_1b)
    return data, plate_df


### End PROBE FITTING  Specific Functions  ### 
