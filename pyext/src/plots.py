"""Various plotting functions for HDX output
"""

from __future__ import print_function
import os
import system_setup
import hdx_models
import numpy
import math
import analysis
from pylab import *
import matplotlib.pyplot as plt
import matplotlib

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

    print(xtics)
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
            #print("NNN", avg_apo[n], avg_lig[n], grid_val_lo, "|", grid_val_lo + tol, grid_val_hi -tol, flag)

            
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
        fig.savefig(outdir + outfile, bbox_inches=0)


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
        plt.show()
        plt.savefig(outdir+outfile, bbox_inches=0, format="png")

def plot_fragment_avg_model_fits(state, sig, frags="all", outfile=None, show_plot=False):
    if frags=="all":
        frags=state.frags

    fig=plt.figure()
    ax=plt.gca()
    for f in frags:

        x=[]
        yavg=[]
        yerror=[]
        #print f.seq, x
        for t in f.timepoints:
            #print t.model
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
                #print f.seq, t.time, t.model_avg, t.model_sd, len(t.models)
        #plt.show()
        plt.errorbar(x, yavg, yerr=yerror)
        ax.set_xscale('log')
        #ax.set_xlim=(1,3600)
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

def plot_model_scores(kinetic_model, show_plot=True, outfile=None):
    '''plots model vs. score for a sorted list'''
    try:
        scores = kinetic_model.model_scores
    except:
        scores = numpy.sort(kinetic_model.calc_model_scores())

    fig = plt.figure()
    x=range(len(scores)) 
    plt.xlabel("Model Rank")
    plt.ylabel("Model Score")
    plt.scatter(x,scores)  

    fig.title=(str(kinetic_model.state.state_name)+"_top_models")

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
        #print f.seq, x
        s1avg=numpy.average(state1.exp_model.get_model_average()[s1frag.start_res+1:s1frag.end_res+1])

        for t in s1frag.timepoints:
            #print t.model
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
                #print(s1frag.seq, t.time, t.model_avg, t.model_sd, len(t.models))
        plt.errorbar(x, yavg, yerr=yerror, c='b')

        if s2frag != "":
            s2avg=numpy.average(state2.exp_model.get_model_average()[s2frag.start_res+1:s2frag.end_res])
            x2=[]
            yavg2=[]
            yerror2=[]
            #print(s2frag, len(s2frag.timepoints))
            for t in s2frag.timepoints:
                # Plot experimental data
                #print(t.models)
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
                    #print f.seq, t.time, t.model_avg, t.model_sd, len(t.models)
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
