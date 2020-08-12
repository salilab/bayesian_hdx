"""Various plotting functions for HDX output
"""

from __future__ import print_function
import os
import system
import numpy
import math
import analysis
#from pylab import *
import matplotlib
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Polygon
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import scipy
#from scipy.stats import gaussian_kde
#from numpy.random import normal
from scipy.integrate import simps
import tools

cdict1 = {'red': ((0.0, 0.0, 0.0),
                 (0.45, 0.0, 0.0),
                 (1.0, 1.0, 1.0)),
        'green': ((0.0, 0.0, 0.0),
                 (0.25, 0.0, 0.0),
                 (0.50, 0.0, 0.0),
                 (0.75, 0.0, 0.0),
                 (1.0, 0.0, 0.0)),
        'blue': ((0.0, 1.0, 1.0),
                 (0.55, 0.0, 0.0),
                 (1.0, 0.0, 0.0))}

def roundup(i, scale):
    return int(math.ceil(i/float(scale)))*scale

def find_minmax(lists):
    flatlist=[item for sublist in lists for item in sublist]
    minmax=(math.floor(min(numpy.array(flatlist))), math.ceil(max(numpy.array(flatlist))) )
    return minmax

def calculate_histogram(list, bins):
    # Given a list of values, calculate a histogram
    return numpy.histogram(numpy.array(list), bins=bins)

def create_sequence_figure_axes(ax, sequence, resis_per_line=40, vertical_buffer=20):
    '''
    Given an pyplot.axis object, format the axis 
    '''

def plot_sequence_overlap(datasets, sequence):
    '''
    For a set of datasets, plot the peptides as bars over a sequence.
    '''
    # Ensure datasets have identical sequences
    seqs = set([d.sequence for d in datasets])
    if len(seqs) > 1:
        raise Exception("plot_sequence_overlap: Sequences in passed datasets are not identical")

    seq = seqs[0]
    n_res = len(seq)

    fig = plt.figure()

def set_logpf_ytick_params(ax):
    # Simple function to add a standard set of ytick parameters
    # whe the axis is log Pf
    ax.set_ylabel('Log Protection Factor')
    ax.tick_params(axis='y', which='both',width=0.5, length=3)
    ax.set_yticks([0,4,8,12])
    ax.set_ylim([-2,10])

def plot_dhdx(ax, diff, z, offset=0, show_plot=True, save_plot=False, outfile="dhdx.png", outdir = ""):
    if outdir != "":
        os.makedirs(outdir, exist_ok=True)

    # Create color bar
    my_cmap=matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict1,256)
    # Normalize colormap to delta log(k) = 5 as maxes
    cNorm = matplotlib.colors.Normalize(vmin=-5., vmax=5.)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=my_cmap)

    # Allow these possible xtics.  Want about 4-6 per plot
    p_xtics = [5,10,20,50,100,200]
    p = list((numpy.array(p_xtics)-len(diff)/4.)**2)
    xtics = p_xtics[p.index(min(p))]

    xin = []
    color=[]
    for n in range(len(diff)):
        
        # Calculate residue color
        rgba=scalarMap.to_rgba(diff[n]) 
        lst=list(rgba)
        
        # calculate saturation (|Z|>3.0) is full saturation, linear scale
        if abs(z[n])>3.0:
            lst[3]=1.0
        else:
            lst[3]=abs(z[n])/3.0
        print(n, lst)
        # Append the xin and color tuples
        xin.append((n+offset+0.5,1))
        rgba=tuple(lst)
        color.append(lst)
    
    yin = (0,1)

    # Plot bar and close outfile
    ax.broken_barh(xin,yin, color=color, lw=0)
    ax.grid(True)

    ax.get_xaxis().set_ticks(range(1,len(diff)+1,xtics))
    ax.get_yaxis().set_ticks([])
    ax.set_xlim([0,len(diff)+1])
    ax.set_ylim([0,1])
    ax.set_xlabel('Residue Number', fontsize=10)

    # Add axes to bottom [left, bottom, width, height] for colormap
    #cax = fig.add_axes([0.3,0.2,0.4,0.02])
    #cb = matplotlib.colorbar.ColorbarBase(cax, cmap=my_cmap, norm=cNorm, spacing='proportional', orientation='horizontal')
    #cb.set_label('\deltaHDX log(k)', fontsize=10)
    
    if show_plot==True:
        fig.show()
    if save_plot==True:
        fig.savefig(outdir + outfile, bbox_inches=0, transparent=True)

    return my_cmap

def plot_apo_lig_dhdx(model, show_plot=True, save_plot=False, outfile="dhdx.png", outdir=None, noclobber=True):
    """ Takes a sampled model and plots the delta HDX (ligand - Apo) 
    as horizontal line plots for each liganded state 
    Also, outputs raw data for each state to individual dat files
    """     

    # Set up output directory for the plots
    if outdir is None:
        outdir = "./output/dhdx_plots/"
    else:
        outdir = outdir + "dhdx_plots/"

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    elif noclobber:
        raise Exception("Output directory ", outdir, " already exists")

    fig=plt.figure(figsize=(8,2*(len(model.states)-1)), dpi=100)
    ax=fig.add_axes([0.1,0.5,0.8,0.4])

    # calculate bounds
    
    # Create color bar
    my_cmap=matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict1,256)
    # Normalize colormap to delta log(k) = 5 as maxes
    cNorm = matplotlib.colors.Normalize(vmin=-5., vmax=5.)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=my_cmap)

    # Apo (reference) state is always state 0
    apo_state=model.states[0]

    # Get plot sizes
    nres = len(apo_state.seq)

    # Allow these possible xtics.  Want about 4-6 per plot
    p_xtics = [5,10,20,50,100,200]
    p = list((numpy.array(p_xtics)-nres/4.)**2)
    xtics = p_xtics[p.index(min(p))]

    maxx = roundup(nres, xtics)

    # Calculate avg and SD log(HDX rate) for apo state; returns vectors of length nres
    avg_apo, sd_apo=analysis.get_average_sector_values(apo_state.exp_model.exp_models, apo_state)

    nlig=-1   # This is the ligand number being processed.  It is the y-value of that state's bar in the plot.

    for s in model.states[1:]:

        # Find the high and low grid values and calculate the
        # tolerance for sectors that are too fast/slow to observe
        grid_val_hi=s.exp_model.exp_grid[-1]
        grid_val_lo=s.exp_model.exp_grid[0]
        tol = 0.05*(grid_val_hi-grid_val_lo)

        # Initialize dat file
        sname = s.state_name
        outdat = outdir + sname + "_dhdx.dat"
        f = open(outdat,'w')
        f.write("Res# Res state_name dhdx dhdx_z avg_lig avg_apo sd_lig sd_apo flag\n")

        nlig=nlig+1     # Increment ligand number
        xin=[]          # Stores residue numbers (x-values of chart)
        yin=(nlig,1.0)  # y value is nlig, bar width is 1.0
        color=[]        # stores color value - associated with xin

        # Calculate avg and SD log(HDX rate) for liganded  state
        avg_lig, sd_lig=analysis.get_average_sector_values(s.exp_model.exp_models, s)

        # Calculate Z-score
        z=analysis.calculate_zscore(apo_state.exp_model.exp_models, s.exp_model.exp_models, apo_state, s)

        # Add state name to plot
        plt.text(nres+1, nlig+0.45, s.state_name, fontsize=14)

        for n in range(len(model.seq)):
            # Write values to file
            dhdx = avg_lig[n] - avg_apo[n]
            flag=""
            if avg_apo[n] < grid_val_lo + tol and avg_apo[n] != 0:
                flag += "*LOW_APO_VAL"
            if avg_lig[n] < grid_val_lo + tol and avg_lig[n] != 0:
                flag += "*LOW_LIG_VAL"
            if avg_apo[n] > grid_val_hi - tol and avg_apo[n] != 0:
                flag += "*HI_APO_VAL"
            if avg_lig[n] > grid_val_hi - tol and avg_lig[n] != 0:
                flag += "*HI_LIG_VAL"

            if avg_lig[n]==0.0 and avg_apo[n]==0.0:
                flag = ""
    
            f.write("%i %s %s %f %f %f %f %f %f %s\n" % 
                        (n,model.seq[n], s.state_name, dhdx, z[n], avg_lig[n], avg_apo[n], sd_lig[n], sd_apo[n], flag))
 
            # Calculate residue color
            rgba=scalarMap.to_rgba(avg_lig[n]-avg_apo[n])  #my_cmap((avg_lig[n]-avg_apo[n])/3.0)
            lst=list(rgba)

            # calculate saturation (|Z|>3.0) is full saturation, linear scale
            if abs(z[n])>3.0:
                lst[3]=1.0
            else:
                lst[3]=abs(z[n])/3.0

            # Append the xin and color tuples
            xin.append((n+model.offset,1))
            rgba=tuple(lst)
            color.append(rgba)

        # Plot bar and close outfile
        ax.broken_barh(xin,yin, color=color, lw=0)
        f.close()
    ax.grid(True)

    ax.get_xaxis().set_ticks(range(0,maxx,xtics))
    ax.get_yaxis().set_ticks([])
    ax.set_xlim([0,maxx])
    ax.set_ylim([0,nlig+1])
    ax.set_xlabel('Residue Number', fontsize=10)
    ax.set_ylabel('Ligand States', fontsize=10)
    ax.set_title('Delta HDX | Target: ' + model.target_name)

    # Add axes to bottom [left, bottom, width, height] for colormap
    cax = fig.add_axes([0.3,0.2,0.4,0.02])
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=my_cmap, norm=cNorm, spacing='proportional', orientation='horizontal')
    cb.set_label('DHDX log(k)', fontsize=10)
    if show_plot==True:
        fig.show()
    if save_plot==True:
        fig.savefig(outdir + outfile, bbox_inches=0, transparent=True)


def plot_fragment_chi_values(state, sig="model", outfile=None, show_plot=False, outdir="./output/"):

    if outdir is None:
        outdir = "./output/fragment_chi_plots/"
    else:
        outdir = outdir + "fragment_chi_plots/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    frags=state.frags
    maxchi=0
    for f in frags:
        if f.get_chi_value(sig) > maxchi:
            maxchi = f.chi

    nres=len(state.seq)
    fig, ax=plt.subplots(figsize=(12,6))
    max_overlap = 50
    color_map = plt.get_cmap('gnuplot')
    cNorm = matplotlib.colors.Normalize(vmin=0, vmax=maxchi)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=color_map)
    sorted_frags = sorted(frags, key=lambda x: x.start_res)
    y_val=1
    maxi=0
    end_res=[0]*max_overlap
    for f in sorted_frags:
        i=1
        while i < max_overlap:
            if f.start_res > end_res[i]:
                y_val=i
                end_res[i]=f.end_res
                break
            i=i+1
        if i> maxi:
            maxi=i
        colorVal=scalarMap.to_rgba(f.chi)
        ax.hlines(y_val, int(f.start_res)-0.5, int(f.end_res)+0.5, color=colorVal, lw=25)
        if nres < 100:
            ax.text(int(f.start_res)+1, y_val-0.3, f.seq)
            ax.text(int(f.start_res)-0.7, y_val, f.start_res)
            ax.text(int(f.end_res)+0.55, y_val, f.end_res)

    # Main Plot Labels
    ax.text(1, maxi+0.5, "Max chi value = " + str(maxchi))
    ax.set_title("Individual Fragment Fits to Model - " + state.state_name)
    ax.set_xlabel('Residue Number')
    ax.set_ylim([0,maxi+1])
    ax.set_xlim([0,nres+1])


    cax = fig.add_axes([0.95,0.2,0.02,0.6])
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=color_map, norm=cNorm, spacing='proportional')
    cb.set_label('Fragment Chi')

    if outfile is None:
        outfile = state.state_name+"_fragment_chi_fits.png"

    if outfile==None:
        plt.show()
    elif show_plot==False:
        plt.savefig(outdir+outfile, bbox_inches=0, format="png")
    else:
        plt.savefig(outdir+outfile, bbox_inches=0, format="png")
        plt.show()

def plot_fragment_avg_model_fits(state, sig, frags="all", outfile=None, show_plot=False):
    if frags=="all":
        frags=state.frags

    fig=plt.figure()
    ax=plt.gca()
    for f in frags:

        x=[]
        yavg=[]
        yerror=[]

        for t in f.timepoints:

            if t.get_model_avg() is not None and t.get_model_sd() is not None:
                x.append(t.time)
                xt=[int(t.time)]*len(t.replicates)
                yt=[float(r.deut) for r in t.replicates]
                plt.scatter(xt, yt)
        chi=f.get_chi_value(sig)
        for t in f.timepoints:
            if t.get_model_avg() is not None and t.get_model_sd() is not None:
                yavg.append(t.model_avg)
                yerror.append(t.model_sd)

        plt.errorbar(x, yavg, yerr=yerror)
        ax.set_xscale('log')

        plt.axis=(1,3600,0,100)
        plt.text(1,1,chi)
        fig.title=(str(f.seq)+"_"+str(chi))
        if outfile==None:
            plt.show()
        elif show_plot==False:
            plt.savefig(outfile, bbox_inches=0)      
        else:
            plt.show()
            plt.savefig(outfile, bbox_inches=0)

def plot_po_model_scores(po, show_plot=True, outfile=None, nscores=0):
    '''plots model vs. score for a sorted list'''
    
    scores = po.get_scores(sorted=True)
    if nscores==0:
        nscores =len(scores)   
    fig = plt.figure()
    x = range(nscores) 
    plt.xlabel("Model Rank")
    plt.ylabel("Model Score")
    plt.scatter(x,scores[:nscores])  

    #fig.title=(str(kinetic_model.state.state_name)+"_top_models")

    if outfile is not None:
        plt.savefig(outfile, bbox_inches=0)
    if show_plot==True:
        plt.show()

def plot_2state_fragment_avg_model_fits(state1, state2, sig, num_best_models=100, outdir=None, write_file=True, show_plot=False):
    ''' For two states (apo and 1 ligand, e.g.), plot the fits to model for all fragments.
        Output an individual plot for each fragment fit
    '''

    for s in [state1, state2]:
        bsm, scores=analysis.get_best_scoring_models(s.modelfile, s.scorefile, num_best_models=num_best_models, prefix=s.state_name, write_file=write_file)
        s.exp_model.import_model_deuteration_from_gridvals(s.frags, bsm)
        s.exp_model.import_models_from_gridvals(bsm)

    #takes a model and score file and writes a new model file with the best X scoring models.
    #This new file can then be imported into an HDXModel class for analysis

    if outdir is None:
        outdir = "./output/fragment_fit-to-data_"+ str(state1.state_name) + "-" + str(state2.state_name) +"/"
    else:
        outdir = outdir + "fragment_fit-to-data_"+ str(state1.state_name) + "-" + str(state2.state_name) +"/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for s1frag in state1.frags:
        # Get the state2 frag.  If this is not there, it will return ""
        s2frag = state2.get_frag(s1frag.seq, s1frag.start_res)

        outfile=str(s1frag.seq) + "_model_fits.png"

        fig=plt.figure()
        ax=plt.gca()
        x=[]
        yavg=[]
        yerror=[]

        s1avg=numpy.average(state1.exp_model.get_model_average()[s1frag.start_res+1:s1frag.end_res+1])

        for t in s1frag.timepoints:

            if t.get_model_avg() is not None and t.get_model_sd() is not None:
                x.append(t.time)
                xt=[int(t.time)]*len(t.replicates)
                yt=[float(r.deut) for r in t.replicates]
                plt.scatter(xt, yt, c='b')
        chi=s1frag.get_chi_value(sig)

        for t in s1frag.timepoints:
            if t.get_model_avg() is not None and t.get_model_sd() is not None:
                yavg.append(t.model_avg)
                yerror.append(t.model_sd)

        plt.errorbar(x, yavg, yerr=yerror, c='b')

        if s2frag != "":
            s2avg=numpy.average(state2.exp_model.get_model_average()[s2frag.start_res+1:s2frag.end_res])
            x2=[]
            yavg2=[]
            yerror2=[]

            for t in s2frag.timepoints:
                # Plot experimental data
                if t.get_model_avg() is not None and t.get_model_sd() is not None:
                    x2.append(t.time)
                    xt=[int(t.time)]*len(t.replicates)
                    yt=[float(r.deut) for r in t.replicates]
                    plt.scatter(xt, yt, c='r')
            chi2=s2frag.get_chi_value(sig)

            for t in s2frag.timepoints:
                #Plot model average and SD errorbars
                if t.get_model_avg() is not None and t.get_model_sd() is not None:
                    yavg2.append(t.model_avg)
                    yerror2.append(t.model_sd)

            plt.errorbar(x2, yavg2, yerr=yerror2, c='r')
            avg = sum(yavg2)/len(yavg2)

            fig.title=(str(s1frag.seq)+"_"+str(chi)+"_"+str(chi2))
        else:
            fig.title=(str(s1frag.seq)+"_"+str(chi))

        ax.set_xscale('log')
        ax.set_xlabel('Time')
        ax.set_ylabel('%D Incorporation')
        #ax.set_xlim=(1,3600)
        plt.axis=(1,3600,0,100)
        plt.text(2,avg+10,s1frag.seq)
        #plt.text(2,avg+5,"Apo chi:"+str(chi))
        #plt.text(2,avg,"Lilly Diff: "+str(lilly_diff))
        #plt.text(2,avg-5,"Sali Diff: "+str(sali_diff))
        #plt.text(1,avg-5,"Sali Exp Diff: "+str(exp_diff))

        if show_plot==False:
            plt.savefig(outdir+outfile, bbox_inches=0, format="png")      
        else:
            plt.show()
            plt.savefig(outdir+outfile, bbox_inches=0, format="png")
        plt.close()
        fig.clear()

def get_cdf(ax,data,pos, bp=False):
    '''
    create violin plots on an axis
    '''
    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.5)
    for d,p in zip(data,pos):
        k = gaussian_kde(d, bw_method="silverman") #calculates the kernel density
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        x = arange(m,M,(M-m)/100.) # support for violin
        v = k.evaluate(x) #violin profile (density curve)
        v = v/v.max()*w #scaling the violin to the available space
        ax.fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
        ax.fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
    if bp:
        ax.boxplot(data,notch=1,positions=pos,vert=1)

def plot_output_figure(oa):
    # Given an OutputAnalysis object, plot the Pf distribution functions
    pass

def plot_peptides(peptides, n_residues, color_chi=False, chi_max=50, cmap='viridis', ax=None, figwidth=10, outfile="peptides.png", show_fig=False, n_term_unobservable=2):
    '''
    Given a list of data.Peptide objects, plot the overlap of these
    @param peptides :: a list of data.Peptide objects
    @param n_residues :: the total number of residues in the system / molecule
    @param color_chi :: Color the peptides by their chi value
    @param chi_max :: the maximum chi to use for the color map
    @param cmap :: the color map to use
    @param ax :: A pyplot.Axes object if you want to plot the figure in a passed object rather than create a new one
    '''

    height_factor = 8

    max_overlap = tools.get_max_overlap(peptides, n_residues)

    figheight = max_overlap / height_factor * figwidth / 10

    # Figure size should be a function of n_residues and max_overlap
    # Use width of 10
    if ax is None:
        fig, ax = plt.subplots(figsize=(figwidth,figheight))

    if color_chi:
        color_map = plt.get_cmap(cmap)
        cNorm = matplotlib.colors.Normalize(vmin=0, vmax=chi_max)
        scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cmap)
        cmap = cmap
    else:
        cmap = matplotlib.cm.binary

    sorted_peptides = sorted(peptides, key=lambda x: (x.start_residue, len(x.sequence)))

    y_val = 1
    maxi = 0
    end_res=[0]*max_overlap
    patches = []
    colors = []
    for pep in sorted_peptides:
        i=1
        while i <= max_overlap:
            if pep.start_residue > end_res[i-1]:
                y_val=i
                end_res[i-1]=pep.start_residue+len(pep.sequence)-1
                break
            i=i+1
        if i > maxi:
            maxi = i
        if color_chi:
            colors.append(min([pep.chi/chi_max, 1.0]))
        else:
            colors.append(1.0)

        patches.append(Rectangle((int(pep.start_residue)-0.5+n_term_unobservable, y_val-0.3), width=len(pep.sequence)-n_term_unobservable, height=0.6))
        #ax.hlines(y_val, int(pep.start_residue), int(pep.start_residue+len(pep.sequence)), color=colorVal, lw=4*figwidth/10)

    if color_chi:
        pc = PatchCollection(patches, cmap=matplotlib.cm.viridis, edgecolor=None)
        pc.set_array(numpy.array(colors))
    else:
        pc = PatchCollection(patches, facecolor='black', edgecolor=None)
    

    ax.add_collection(pc)
    # Main Plot Labels
    #ax.text(1, maxi+0.5, "Max chi value = " + str(chi_max))
    ax.set_xlabel('Residue Number')
    ax.set_ylabel('Peptides')

    ax.set_ylim([0,max_overlap+0.5])
    ax.set_xlim([0,n_residues+1.0])
    ax.set_yticks([])

    #if color_chi:
    #    #cax = fig.add_axes([0.95,0.2,0.02,0.6])
    #    #cb = matplotlib.colorbar.ColorbarBase(cax, cmap=color_map, norm=cNorm, spacing='proportional')
    #    #cb.set_label('Fragment Chi')

    if outfile is not None:
        plt.tight_layout()
        plt.savefig(outfile, bbox_inches=0, format="png")

    if show_fig:
        plt.tight_layout()
        plt.show()


def plot_system_data_information_content(arr, protection_factors, figwidth=10, figheight=4, 
                                        ax = None, offset=0, outfile=None, show=False):
    '''
    Plots the system information content
    @param arr :: 2D numpy array of shape n_residues x n_protection_factors
    @param 
    '''
    if ax is None:
        fig, ax = plt.subplots(figsize=(figwidth,figheight))

    max_arr = numpy.max(arr)

    patches = []
    colors = []

    for res in range(len(arr)):
        xstart = offset + res + 0.5
        for pf in range(len(protection_factors)):
            if pf == 0:
                ystart = protection_factors[pf] - 0.5*(protection_factors[pf+1]-protection_factors[pf])
            else:
                ystart = protection_factors[pf] - 0.5*(protection_factors[pf]-protection_factors[pf-1])

            if pf == len(protection_factors)-1:
                height = protection_factors[pf]-protection_factors[pf-1]
            else:
                height = 0.5*(protection_factors[pf]-protection_factors[pf-1]) + 0.5*(protection_factors[pf+1]-protection_factors[pf])
            
            patches.append(Rectangle((xstart, ystart), width=1.0, height=height))
            colors.append(arr[res][pf]/max_arr)

    pc = PatchCollection(patches, cmap=matplotlib.cm.binary, edgecolor=None)
    pc.set_array(numpy.array(colors))

    ax.add_collection(pc)

    ax.set_xlabel('Residue Number')

    set_logpf_ytick_params(ax)

    #ax.set_yticklabels(["0","4","8"])
    
    #ax.set_ylim([math.floor(protection_factors[0]),10])#math.ceil(protection_factors[-1])])
    ax.set_xlim([0,len(arr)+1.0])
    xtick_skip = math.ceil(len(arr)/80)*10
    xticks = [1]
    xt = xtick_skip
    while xt < len(arr):
        xticks.append(xt)
        xt+=xtick_skip
    ax.set_xticks(xticks)

    ax.grid(b=True, which='major', axis = "y", color='black', linestyle=':', alpha = 0.5)
    ax.set_frame_on(False)

    if ax is None:
        plt.tight_layout()
    if show:
        plt.show()
    if outfile is not None:
        plt.savefig(outfile, dpi=300)


def plot_sectors(state, axes):
    # given one or more axes, plot a dashed line delineating the sectors
    # sectors are determined by the given peptides.
    sectors = state.get_sectors()
    if len(sectors) == 0:
        sectors = state.calculate_sectors()

    for s in sectors:
        sr = list(s.get_residues())
        for ax in axes:
            ylims = ax.get_ylim()
            ax.plot([min(sr)-0.5, min(sr)-0.5], ylims, ls=":", lw=0.5, color="darkgreen")
            ax.plot([max(sr)+0.5, max(sr)+0.5], ylims, ls=":", lw=0.5, color="darkgreen")

    #sr = list(sectors[-1].get_residues())
    #for ax in axes:
    #    ylims = ax.get_ylim()
    #    ax.plot([max(sr)+0.5,max(sr)+0.5], ylims, ls=":", lw=0.5, color="darkgreen")        

def plot_pf_distributions_and_peptide_chis(pof, width_scale = 1.0, outfile="pfs_and_chis.png", dpi=300):
    '''
    Given a POF object, plot the peptide overlap, colored by chi value
    and the peptide distribution profiles.

    TODO: Add in priors. Calculate posterior shrinkage.
    '''
    pep_chis = pof.calculate_peptide_chis()
    peptides = pof.get_all_peptides()

    for p in range(len(peptides)):
        peptides[p].chi = pep_chis[p]

    max_overlap = tools.get_max_overlap(peptides, len(pof.sequence))

    fig = plt.figure(figsize=(10*width_scale,max_overlap/4.0 + 2))
    gs = gridspec.GridSpec(2, len(pof.sequence), height_ratios=[max_overlap/20.0, 1]) 
    gs.update(wspace=0.0, hspace=0.1)
    ax0 = fig.add_subplot(gs[0, :])

    # Put residue lines on the peptides plot
    xt = [1] + list(range(20,len(pof.sequence), 20))
    for x in xt:
        ax0.plot([x,x],[-1,max_overlap+0.5], color="black", lw=0.5)

    #ax0.set_xticks(xt)

    # Figure out way to determine residue ticks
    xtickspacing = 20

    plot_peptides(pof.get_all_peptides(), len(pof.sequence), color_chi=True, chi_max=math.ceil(max(pep_chis)), outfile=None, ax=ax0)
    ax0.axis('off')

    ax = plot_residue_protection_factors([pof], gridspec=gs, fig=fig, return_ax=True, resnum_skip=20, resrange=(1,340))

    set_logpf_ytick_params(ax[0])

    plt.savefig(outfile, dpi=dpi)
    

def plot_overlap_and_information(state, protection_factors=None, sectors=True, outfile=None, dpi=300, show=False, figwidth_scale=1.0):
    # Given these datasets and protection factors, plot a nice figure with peptide overlap on top and
    # Residue Pf information content on the bottom
    

    if protection_factors is None:
        protection_factors = numpy.arange(-2,14,0.1)

    datasets = state.get_datasets()

    info = numpy.zeros((len(state.sequence), len(protection_factors)))
    peptides = []

    for ds in datasets:
        peptides += ds.get_peptides()
        info += ds.get_residue_information_content(protection_factors)

    max_overlap = tools.get_max_overlap(peptides, len(state.sequence))

    fig = plt.figure(figsize=(10*figwidth_scale,max_overlap/4.0 + 2))
    gs = gridspec.GridSpec(2, 1, height_ratios=[max_overlap/20.0, 1]) 
    ax0 = plt.subplot(gs[0])

    plot_peptides(peptides, len(state.sequence), outfile=None, ax=ax0)
    ax0.axis('off')

    ax1 = plt.subplot(gs[1], sharex=ax0)
    plot_system_data_information_content(info, protection_factors, ax=ax1)
    plt.setp(ax0.get_xticklabels(), visible=False)
    yticks = ax1.yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)

    if sectors:
        plot_sectors(state, [ax0, ax1])

    cb_height = 1/(max_overlap/20.0 + 1) - 0.2

    cax = fig.add_axes([0.92,0.15,0.02,cb_height])
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=matplotlib.cm.binary, norm=matplotlib.colors.Normalize(vmin=0, vmax=numpy.max(info)), spacing='proportional')
    cb.set_ticks([0,numpy.max(info)])
    cb.set_ticklabels(["0","+"])
    cb.set_label('Information Content')

    ax0.set_title(state.macromolecule.name+" "+state.name +" peptide map and residue information content")

    plt.gcf().subplots_adjust(bottom=0.15)
    plt.subplots_adjust(hspace=.0)

    if outfile is not None:
        plt.savefig(outfile, dpi=dpi)

    if show:
        plt.show()

    return fig


def calculate_shannon_bits(hist):
    # From a discrete probability distribution, represented as
    # a numpy array of probabilities, calculate the shannon information
    # gain over a uniform distribution

    if sum(hist)==0:
        return 0

    nbins = len(hist)
    base_info = numpy.log(nbins)
    hist_info = 0

    for p in hist:
        if p != 0:
            hist_info += p * numpy.log(1/p)

    return base_info - hist_info

def plot_residue_rate_distributions(model_files, rate_bins = None, resrange=None, plot_prior=True):
    # Input is a standard output model file and the rate bins
    # Should note the sectors as well, along with overlap.
    # Outputs a plot with a sing
    import csv

    colors = ["red", "blue", "yellow", "green", "orange", "darkred"]

    resnum_label_skip=10

    # Get data and place into list of lists.
    
    if type(model_files) != list:
        model_files = [model_files]

    d_list = []
    for i in model_files:
        d_list.append(numpy.loadtxt(i))

    nres = len(d_list[0][0])
    nmod = len(d_list[0])
    maxbin = int(numpy.max(d_list[0]))

    # How to calculate xlim? Keep the uniform prior the same proportion of the window
    # So proportional to 1/bins.  Say 3 or 4 times this value?
    xlim = 1.0/maxbin * 3

    plt.subplots_adjust(wspace=0)

    if resrange is None:
        resrange = range(1, nres+1)
    else:
        resrange = range(resrange[0], resrange[1]+1)

    bins = range(1,maxbin+1)

    # What is the optimal figsize?
    fig, ax = plt.subplots(1, len(resrange), sharey='row', figsize=(20,2))

    data = []
    if rate_bins is not None:
        if len(rate_bins) < maxbin:
            raise Exception("Number of inputted rate_bins is less than the maximum bin in output file.")
        x = bins
    else:

        x = bins

    # Calculate the histograms
    for n in resrange:
        d_hists=[]
        for d in d_list:
            nums = d[:,n]
            h = calculate_histogram(list(nums), numpy.array(range(maxbin+1))+0.5)
            hist = 1.0 * h[0] / nmod
            d_hists.append(hist)
        data.append(d_hists)
        # data is a list of lists.  Outer index is resnum, inner index is the dataset.

    # Figure out some way to determine the best ytick rate
    ytick_rate=2

    # Outer loop over all residues
    for nd in resrange:
        # n is the plot index; nd is the residue index
        n = nd - resrange[0]
        x_lists = data[n] 
        for i in range(len(x_lists)):
            xl = x_lists[i]
            arr = numpy.array(xl)

            ax[n].set_xticks([]) 

            if nd%resnum_label_skip == 0:
                ax[n].set_title(str(nd), fontsize=12)
                ax[n].set_xticks([0]) 
                ax[n].xaxis.set_ticks_position("top")
                ax[n].plot([0,0],[x[0]-0.2,x[-1]+0.2], color="black", lw=0.5)
            
            # Calculate bits of information
            bits = calculate_shannon_bits(arr)

            ax[n].set_xlim((-xlim,xlim))
            ax[n].set_ylim((-numpy.log(len(x))+1,x[-1]+0.2))

            if plot_prior:
                ax[n].fill_betweenx(x,0,1.0*numpy.ones(len(x))/len(x),facecolor='grey',alpha=0.5, lw=0)
                ax[n].fill_betweenx(x,0,-1.0*numpy.ones(len(x))/len(x),facecolor='grey',alpha=0.5, lw=0)
            
            if sum(arr) != 0:
                # Fill in the prior probability (uniform for now)
                ax[n].fill_betweenx(x,0,arr,facecolor=colors[i],alpha=0.5, lw=0)
                ax[n].fill_betweenx(x,0,-arr,facecolor=colors[i],alpha=0.5, lw=0)
                # Add in lower bar for information content
                ax[n].barh(y=-1*bits+1,width=2*xlim/len(x_lists), height=bits,left=-1*xlim+i*2*xlim/len(x_lists),color=colors[i], alpha=0.7, lw=0)
                for yval in range(1, int(numpy.max(d))+1,ytick_rate):
                    ax[n].axhline(y=yval,ls='-', lw=0.5)

                #ax[n].set_xticks([str(nd)])       
            ax[n].tick_params(axis='x', which='major', labelsize=0, color="grey")
            
            
            ax[n].set_frame_on(False)
            #ax[n].spines['top'].set_visible(False)
            #ax[n].spines['right'].set_visible(False)
            #ax[n].spines['bottom'].set_visible(False)
            #ax[n].spines['left'].set_visible(False)

    ax[0].set_ylabel("HX Rate Bin")
    ax[0].set_yticks(x[::ytick_rate])
    #ax[0].tick_params(axis='y', which='major', labelsize=8)

    plt.savefig("test_violins.png", dpi=300, format="png")
    plt.show()

def plot_priors(sfunc, n_residues, protection_factors=None, outfile="priors.png", dpi=300):
    if protection_factors is None:
        protection_factors = numpy.arange(-2,10,0.1)

    if not hasattr(sfunc, '__iter__'):
        sfunc = [sfunc]
    colors = ["red", "blue", "yellow", "green", "black"]

    prior_likelihoods = []
    for i in range(1,n_residues+1):
        like = numpy.empty(len(protection_factors))
        for pf in range(len(protection_factors)):
            like[pf] = sfunc[0].get_prior_likelihood(i, protection_factors[pf])
        prior_likelihoods.append(like/numpy.linalg.norm(like))

    if len(sfunc) == 1:
        fig,ax = plot_violin_distributions(prior_likelihoods, protection_factors, color="black")
    else:
        fig,ax = plot_violin_distributions(prior_likelihoods, protection_factors, color=colors[0])
        for j in range(1, min(len(sfunc), len(colors))):
            prior_likelihoods = []
            for i in range(1,n_residues+1):
                like = numpy.empty(len(protection_factors))
                for pf in range(len(protection_factors)):
                    like[pf] = sfunc[j].get_prior_likelihood(i, protection_factors[pf])
                prior_likelihoods.append(like/numpy.linalg.norm(like))
            fig, ax = plot_violin_distributions(prior_likelihoods, protection_factors, ax=ax, color=colors[j])

    plt.savefig(outfile, dpi=dpi)

def plot_violin_distributions(dists, yvals, figwidth=10, figheight=2, index_range=None, ax=None, label_skip=20, color="black", alpha=0.5):
    '''
    Given a 2D matrix of distributions, plot them as individual violin plots
    '''

    xlim = numpy.nanmax(dists)
    
    res_ixs = range(len(dists))
    if index_range is not None:
        if index_range[1] < len(dists):
            res_ixs = range(index_range[0], index_range[1]+1)
        else:
            raise Exception("index range must be contained in residue range", len(dists))

    if ax is None:
        fig, ax = plt.subplots(1, len(res_ixs), sharey='row', figsize=(figwidth,figheight))
    else:
        fig = ax[0].get_figure()

    for r in res_ixs:
        ax[r].set_xticks([])
        ax[r].set_yticks([0,4,8])
        ax[r].tick_params(axis='y', which='both',width=0.15, length=2)

        dist_vals = dists[r] 


        if r%label_skip == label_skip-1:
            #ax[r].set_title(str(r+1), fontsize=11)
            ax[r].set_xticks([0]) 
            ax[r].set_xticklabels([str(r+1)])
            ax[r].xaxis.set_ticks_position("both")
            #ax[r].plot([0,0],[yvals[0]-0.2,yvals[-1]+0.2], color="black", lw=0.2)
            pass
        # Calculate bits of information
        #bits = calculate_shannon_bits(arr)

        # Set x and y limits
        ax[r].set_xlim((-xlim,xlim))
        #ax[r].set_ylim((-numpy.log(len(yvals))+1,yvals[-1]+0.2))
        ax[r].set_ylim((numpy.floor(yvals[0]), numpy.ceil(yvals[-1])))
        # Add blocks for observable windows
        #window_threshold = 0.005 # what % of deuteration is considered observable?

        if not math.isnan(dist_vals[0]):
            ax[r].fill_betweenx(yvals,0,dist_vals,facecolor=color, alpha=alpha, lw=0)
            ax[r].fill_betweenx(yvals,0,-1*dist_vals,facecolor=color, alpha=alpha, lw=0)   
            # Add in lower bar for information content
            #ax[n].barh(y=-1*bits+minbin,width=2*xlim/len(x_lists), height=bits,left=-1*xlim+i*2*xlim/len(x_lists),color=colors[i], alpha=0.7, lw=0)
            #for yval in yvalsrange(minbin, maxbin, ytick_rate):
            #    ax[r].axhline(y=yval,ls='-', lw=0.5)
            pass
        else:
            ax[r].fill_betweenx(yvals,0,numpy.zeros(len(dists[0])),facecolor=color,alpha=alpha, lw=0)
            #ax[n].fill_betweenx(x,0,-arr,facecolor=colors[i],alpha=0.5, lw=0)  
            pass
            #ax[n].set_xticks([str(nd)])       
        #ax[r].tick_params(axis='x', which='major', labelsize=0, color="grey")
        
        ax[r].set_frame_on(False)

    ax[0].set_xticks([0])
    ax[0].set_xticklabels(["1"])
    ax[0].tick_params(axis='y', which='both',width=0.1, length=3)
    ax[0].set_ylabel('Log Protection Factor')
    ax[0].set_yticks([0,4,8])
    ax[0].set_yticklabels(["0","4", "8"])


    #ax[0].tick_params(axis='y', which='major', labelsize=8)

    return fig, ax

def plot_rhat_from_output_file(output_file, ax=None, show=False, outfile=None):

    oa = analysis.OutputAnalysis([output_file])
    psrf = numpy.array(oa.calculate_rhat())
    if ax is None:
        fig = plt.figure(figsize=(6,3))
        ax = fig.gca()
    plot_rhat(psrf, ax)

    ax.set_title(output_file+" residue convergence")

    if show:
        plt.show()

    if outfile is not None:
        plt.savefig(outfile, transparent=True)

    return ax


def plot_rhat(rhats, ax=None, color=None, show=False, outfile=None, dpi=300):

    colors = ["black", "red", "green", "blue", "pink"]

    n_r = len(rhats)

    if n_r > 5:
        raise Exception("Not a good idea to plot more than 5 Rhats at once")

    if n_r == 1:
        offsets = [0]
    else:
        offsets = [-0.1 + i*0.2/(n_r-1) for i in range(n_r)]

    if ax is None:
        fig = plt.figure(figsize=(6,3))
        ax = fig.gca()

    ax.set_ylim([0.9, 2.0])
    ax.set_yticks([1.0, 1.1, 1.3, 1.5, 2.0])
    ax.set_xlim([0,len(rhats[0])+1])
    xt = [1] + list(range(20,len(rhats[0]), 20))
    ax.set_xticks(xt)
    ax.set_ylabel(r'$\hat{R}$')
    ax.set_xlabel("Residue Number")
    ax.grid(b=True, which='major', axis = "y", color='black', linestyle=':', alpha = 0.5)

    if n_r == 1:
        for r in range(len(rhats[0])):
            ax.plot([r+1,r+1],[0,rhats[0][r]], lw=0.5, color=colors[0])
            ax.plot(r+1,rhats[0][r], marker='o', markersize=2, color=colors[0], lw=0.0)
    else:
        for rh in range(n_r):
            for r in range(len(rhats[rh])):
                ax.plot([r+1+offsets[rh],r+1+offsets[rh]],[0,rhats[rh][r]], lw=0.5, color=colors[rh])
                ax.plot(r+1+offsets[rh],rhats[rh][r], marker='o', markersize=2, color=colors[rh], lw=0.0)

    plt.tight_layout()

    if show:
        matplotlib.use("terminal")
        plt.show()
        matplotlib.use("agg")

    if outfile is not None:
        plt.savefig(outfile, transparent=True, dpi=dpi)

    return ax
        


def plot_residue_protection_factors(parse_output, rate_bins=None, 
                                    resrange=None, plot_prior=True, 
                                    resnum_skip=10, num_best_models=100,
                                    show=False, sort_sectors=False,
                                    true_vals=None, outputdir=None,
                                    gridspec=None, fig=None, return_ax=False):
    # Input is a standard output model file and the rate bins
    # Should note the sectors as well, along with overlap.
    # Outputs a plot with a sing
    import csv

    colors = ["red", "blue", "yellow", "green", "orange", "darkred"]

    resnum_label_skip = resnum_skip

    # Get data and place into list of lists.
    
    if type(parse_output) != list:
        parse_output = [parse_output]

    data_list = []
    po_states = ""

    for po in parse_output:
        po_states+=po.state_name+"_"
        bsm = po.get_best_scoring_models(num_best_models, return_pf=True, sort_sectors=sort_sectors)
        if len(bsm) > 10:
            data_list.append(po.get_best_scoring_models(num_best_models, return_pf=True, sort_sectors=sort_sectors))

    # data_list contains a list of a list of models.  
    #   The first index is the number of ParseOutput objects
    #   The second index is the model number
    #   The third index is 0 for the score and 1 for the model
    #   The fourth index is always zero, since the model is a list of a lists. Keeping this for potential multi-state
    #   The fifth index is the residue in the model. 
    nres = len(data_list[0][1][1][0])
    nmod = len(data_list[0])

    maxbin = 14#int(numpy.max(maxarr[~numpy.isnan(maxarr)]))
    #minarr = numpy.floor(data_list[0][0][1])
    minbin = -2#int(numpy.min(minarr[~numpy.isnan(minarr)]))

    pf_list = []

    # pf_list is an array of models
    for d in data_list:
        pfs = []
        # Loop over all models in this list     
        for i in d:
            # Append just the model itself. 
            pfs.append(numpy.array(i[1][0]))

        pf_list.append(numpy.array(pfs))

    # How to calculate xlim? Keep the uniform prior the same proportion of the window
    # So proportional to 1/bins.  Say 3 or 4 times this value?
    xlim = 1.0/maxbin * 4

    plt.subplots_adjust(wspace=0)

    if resrange is None:
        resrange = range(1, nres+1)
    else:
        resrange = range(resrange[0], resrange[1]+1)

    bins = range(1,maxbin+1)

    # What is the optimal figsize?
    if gridspec is None:
        fig, ax = plt.subplots(1, len(resrange), sharey='row', figsize=(10,2))
    else:
        ax = []
        for i in resrange:
            ax.append(fig.add_subplot(gridspec[1,i-1]))

    data = []
    if rate_bins is not None:
        if len(rate_bins) < maxbin:
            raise Exception("Number of inputted rate_bins is less than the maximum bin in output file.")
        x = bins
    else:
        x = numpy.linspace(minbin,maxbin,parse_output[0].grid_size)

    # Calculate the histograms over all rate bins for each residue in the desired range
    for n in resrange:
        d_hists=[]
        for d in pf_list:
            nums = d[:,n-1]
            if math.isnan(nums[0]):
                d_hists.append(numpy.zeros(parse_output[0].grid_size))
            else:
                h = calculate_histogram(list(nums), numpy.linspace(minbin-0.0001,maxbin+0.0001,parse_output[0].grid_size+1))
                hist = 1.0 * h[0] / nmod
                d_hists.append(hist)
        data.append(d_hists)
        # data is a list of lists.  Outer index is resnum, inner index is the dataset.

    # Figure out some way to determine the best ytick rate
    ytick_rate=20

    # Outer loop over all residues
    for nd in resrange:
        # n is the plot index; nd is the residue index
        n = nd - resrange[0]
        x_lists = data[n]  
        # Loop over all datasets 
        for i in range(min(len(x_lists),5)):
            xl = x_lists[i] 
            arr = numpy.array(xl) # arr is the histogram
            
            ax[n].set_xticks([]) 

            if nd%resnum_label_skip == 0 or nd==1:
                ax[n].set_title(str(nd), fontsize=10)
                ax[n].set_xticks([]) 
                ax[n].xaxis.set_ticks_position("top")
                ax[n].plot([0,0],[x[0]-0.2,x[-1]+0.2], color="black", lw=0.5)
            
            # Calculate bits of information
            bits = calculate_shannon_bits(arr)

            # Set x and y limits
            ax[n].set_xlim((-xlim,xlim))
            ax[n].set_ylim((-2,10))
            ax[n].set_yticks([]) 

            # Add blocks for observable windows
            window_threshold = 0.05 # what % of deuteration is considered observable?

            if plot_prior:
                # Fill in the prior probability (uniform for now)
                ax[n].fill_betweenx(x,0,1.0*numpy.ones(len(x))/len(x),facecolor='grey',alpha=0.5, lw=0)
                ax[n].fill_betweenx(x,0,-1.0*numpy.ones(len(x))/len(x),facecolor='grey',alpha=0.5, lw=0)

            if not math.isnan(arr[0]): 
                ax[n].fill_betweenx(x,0,arr,facecolor=colors[i],alpha=0.5, lw=0)
                ax[n].fill_betweenx(x,0,-arr,facecolor=colors[i],alpha=0.5, lw=0)   
                # Add in lower bar for information content
                #ax[n].barh(y=-1*bits+minbin,width=2*xlim/len(x_lists), height=bits,left=-1*xlim+i*2*xlim/len(x_lists),color=colors[i], alpha=0.7, lw=0)
                for yval in range(minbin, maxbin, ytick_rate):
                    ax[n].axhline(y=yval,ls='-', lw=0.5)
            else:
                ax[n].fill_betweenx(x,0,numpy.zeros(parse_output[0].grid_size-0),facecolor=colors[i],alpha=0.5, lw=0)
                #ax[n].fill_betweenx(x,0,-arr,facecolor=colors[i],alpha=0.5, lw=0)  

                #ax[n].set_xticks([str(nd)])       
            ax[n].tick_params(axis='x', which='major', labelsize=0, color="grey")
            
            ax[n].set_frame_on(False)

        if true_vals is not None:
            if nd-1 in true_vals.keys():
                ax[n].plot((-1, 1), (true_vals[nd-1], true_vals[nd-1]), "k-",  color='limegreen', lw=2)
            # How to sort these??

    if return_ax:
        return ax

    ax[0].set_ylabel("Log(Pf)")
    ax[0].set_yticks([0,4,8])
    #ax[0].tick_params(axis='y', which='major', labelsize=8)

    if outputdir is not None:
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        if sort_sectors:
            plt.savefig(outputdir+"/protection_factor_distributions_sorted.png", dpi=600, format="png", transparent=True)
        else:
            plt.savefig(outputdir+"/protection_factor_distributions.png", dpi=600, format="png", transparent=True)
    else:
        plt.savefig("protection_factor_distributions-"+po_states[:-1]+".png", dpi=600, format="png", transparent=True)
    if show:
        plt.show()

def plot_sampling_convergence_stats(pvals, cvs, percents, cutoff_list, sampling_precision=None, show=False, outfile=None, dpi=300):
    fig = plt.figure()
    ax = fig.gca()

    ax.plot(cutoff_list, pvals, marker='.', color="orange")
    ax.plot(cutoff_list, cvs, marker='X', color="green")
    ax.plot(cutoff_list, numpy.array(percents)/100, marker='^', color="purple")
    ax.set_xlabel("Sampling Precision")

    if sampling_precision:
        ax.plot([sampling_precision, sampling_precision], [0,1], marker=None, lw=1.0, color="black", ls=":")

    if outfile is not None:
        plt.savefig(outfile, dpi=dpi, transparent=True)
    if show:
        plt.show()

def plot_incorporation_curve_fits(po, num_models, write_plots=True, single_plot=False, output_directory=None, log_time=False):
    '''
    For each dataset:
      For each peptide:
         1) Calculate the mean, SD and chi of deuterium incorporation 
            over all models at each timepoint
         2) Plot the mean/SD along with experimental values
         3) Write plot to the output directory
    '''
    
    state_name_prefix = po.state_name
    print(" --- Plotting incorporation curve fits for state", state_name_prefix)


    if output_directory is None:
        output_directory = "state_"+state_name_prefix+"/incorporation_plots/"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Get the best scoring models as protection factors
    # Returns as a list of tuples [(score, [pfs])]
    bsms = po.get_best_scoring_models(N=num_models, return_pf=True)

    # Get the list of datasets used to compute the models
    datasets = po.get_datasets()

    # Loop over all datasets
    for d in datasets:

        # Get a list of all timepoints in the dataset sorted from low to high
        sorted_times = numpy.sort(list(d.get_all_times()))

        # Get the experimental data as a dictionary
        exp_deuteration_by_peptide_time = d.get_dataset_as_dictionary(log_time)

        # Open and initialize peptide chi file, which lists all peptides 
        # and the chi^2 value for each timepoint
        chi_out = open(output_directory+"/peptide_chi_values.dat", "w")
        init_string = ""
        for tp in sorted_times:
            init_string += str(tp)+"s, "
        chi_out.write("Pep_seq, start_res, chi^2 | "+init_string[:-2]+"\n")

        #print("  Peptide | Chi^2")
        for pep in d.get_peptides():

            # pep_id is the data_dict identifier for each peptide
            pep_id = pep.sequence+"_"+str(pep.start_residue)

            tot_chi = 0
            chi_string = ""

            # Get model incorporations over all good scoring models
            model_deuteration_by_time = {}
            extra_times = []
            if log_time:
                extra_log_times = numpy.arange(1, 6, 0.2)
                extra_times = [int(10**e)for e in extra_log_times]

            all_times = numpy.sort(list(set(extra_times+list(sorted_times))))

            for time in all_times:
                if log_time:
                    t = numpy.log(time)
                else:
                    t = time
                model_deuteration_by_time[t] = []

            for gsm in bsms:
                # Extract GSM protection factors
                protection_factors = gsm[1][0]

                for time in all_times:
                    incorp = (tools.get_timepoint_deuteration(pep, time, protection_factors) / pep.get_number_of_observable_amides())*100
                    if log_time:
                        t = numpy.log(time)
                    else:
                        t = time
                    model_deuteration_by_time[t].append(incorp)

            n_tps = len(sorted_times)
            for time in sorted_times:
                tp = pep.get_timepoint_by_time(time)
                if tp is None:
                    n_tps-=1
                    continue
                if log_time:
                    time = numpy.log(tp.time)
                diff = -1*(numpy.average(exp_deuteration_by_peptide_time[pep_id][time])-numpy.average(model_deuteration_by_time[time]))
                chi = diff**2/numpy.average(exp_deuteration_by_peptide_time[pep_id][time])
                tot_chi+=chi
                chi_string += str(chi)+", "
            if n_tps > 0:
                tot_chi = tot_chi / float(n_tps) 
            else:
                tot_chi = -1
            
            chi_out.write(pep.sequence +", " +str(pep.start_residue) +", " +str(tot_chi) + " | "+chi_string[:-2])

            fig, ax = plot_incorporation_curve_spread(model_deuteration_by_time, exp_deuteration_by_peptide_time[pep_id])
            #print("  "+str(pep.start_residue)+"_"+pep.sequence+"\t|\t %.2f2" % tot_chi)
            ax.set_title(str(pep.start_residue) + "-" + pep.sequence)
    
            if log_time:
                ttime = 1
                ax.text(1, 95, "Chi^2 = %.2f2" % tot_chi)
            else:
                ttime = sorted_times[1]
            ax.text(ttime, 95, "Chi^2 = %.2f2" % tot_chi)

            if write_plots:
                try:
                    os.mkdir(output_directory)
                except:
                    pass
                try:
                    os.mkdir(output_directory  + state_name_prefix)
                except:
                    pass
            
                fig.savefig(output_directory  + str(pep.start_residue) + "_" + pep.sequence + "_NM" + str(num_models) +".png", transparent=True)
            plt.close(fig)

def plot_incorporation_curve(deuteration_by_time, exp, ax=None, color='blue', plot=False):
    # Given a set of experimental and model deuteration values
    # Plot each of these and return the fig, ax objects.
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()


    ax.set_xlim([0,max(deuteration_by_time.keys())*1.1])
    ax.set_ylim([0,100])
    ax.set_ylabel("%Deuterium incorporation")
    ax.set_xlabel("time (s)")

    #ax.get_xaxis().set_ticks(range(0,maxx,xtics))
    ax.get_yaxis().set_ticks([0,20,40,60,80,100])

    for time in exp.keys():
        for value in exp[time]:
            # Plot each experimental timepoint
            ax.scatter(time, value, c=color)

    for time in deuteration_by_time.keys():
        avg = numpy.average(deuteration_by_time[time])
        sd = numpy.std(deuteration_by_time[time])
        ax.plot(time, avg, c='red', lw = 0.2)
        ax.errorbar(time, avg, yerr=sd, c='red')

    if plot:
        plt.show()

    return fig, ax

def plot_incorporation_curve_spread(deuteration_by_time, exp, ax=None, color='blue', plot=False):
    # Given a set of experimental and model deuteration values
    # Plot each of these and return the fig, ax objects.
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()


    ax.set_xlim([0,max(deuteration_by_time.keys())*1.1])
    ax.set_ylim([0,100])
    ax.set_ylabel("%Deuterium incorporation")
    ax.set_xlabel("time (s)")

    #ax.get_xaxis().set_ticks(range(0,maxx,xtics))
    ax.get_yaxis().set_ticks([0,20,40,60,80,100])

    # Make 5 lines for the quantiles, 0, 0.25, 0.5, 0.75, 1.0
    qs = [0.0, 0.25, 0.5, 0.75, 1.0]
    times = sorted(list(deuteration_by_time.keys()))
    q_data = {}
    for q in qs:
        q_data[q] = []
        for time in times:
            q_data[q].append(numpy.quantile(deuteration_by_time[time], [q])[0])

        if q == 0.0 or q == 1.0:
            ax.plot(times, q_data[q], c='red', lw = 0.5, linestyle='dashed')
        #elif q == 0.25 or q==0.75:
        #    ax.plot(times, q_data[q], c='red', lw = 0.4)
        elif q==0.5:
            ax.plot(times, q_data[q], c='red', lw = 1.0)

    for time in exp.keys():
        ax.scatter(time, numpy.mean(deuteration_by_time[time]), c='red', marker='o', s=12)

    plt.fill_between(times, q_data[qs[1]], q_data[qs[3]], facecolor='red', alpha=0.2)

    for time in exp.keys():
        for value in exp[time]:
            # Plot each experimental timepoint
            ax.scatter(time, value, c=color, s=10, marker="D")

    if plot:
        plt.show()

    return fig, ax

def plot_peptide_avg_model_fits(dataset, model_pfs, num_models=100, outfile=None, show_plot=False):
    # Given a dataset and a set of models, plot the fit to the experimental data

    fig = plt.figure()
    ax = plt.gca()

    # First, calculate the average and SD of Deuterium incorporation at each peptide,
    # timepoint

    model_deuts = {}

    for pep in dataset.get_peptides():
        pep_deuts = {}
        for tp in pep.get_timepoints():
            tp_deuts = []
            for pfs in model_pfs:
                tp_deuts.append(tools.get_timepoint_deuteration(peptide, tp.time, pfs))
            pep_deuts[tp.time] = tp_deuts

        model_deuts[pep.sequence] = pep_deuts

    for pep in dataset.get_peptides():

        pep_deut = model_deuts[pep.sequence]

        x=[]
        yavg=[]
        yerror=[]

        for tp in pep.get_timepoints():

            x.append(t.time)
            xt=[int(t.time)]*len(tp.get_replicates())
            yt=[float(r.deut) for r in t.get_replicates()]
            plt.scatter(xt, yt)

            yavg.append(numpy.average(pep_deut[tp.time]))
            yerror.append(numpy.stdev(pep_deut[tp.time]))

        plt.errorbar(x, yavg, yerr=yerror)
        ax.set_xscale('log')
        ax.set_ylabel("%D Incorporation")
        ax.set_xlabel("Time (seconds)")

        plt.axis=(1,3600,0,100)
        #plt.text(1,1,chi)
        fig.title=(dataset.name +"_"+pep.sequence)#+"_"+str(chi))
        if outfile==None:
            plt.show()
        elif show_plot==False:
            plt.savefig(outfile, bbox_inches=0, transparent=True)      
        else:
            plt.show()
            plt.savefig(outfile, bbox_inches=0, transparent=True)
