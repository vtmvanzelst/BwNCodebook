# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 17:10:29 2023

@author: zelst
"""

import math
import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import join
import pandas as pd
import pickle
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes
import shutil
import subprocess
import xarray as xr
import sys


def XB_load_results(dir_nc, vegopt = True):
    ds     = xr.open_dataset(dir_nc)
    
    x       = ds['globalx'].values.squeeze()  # xcoords
    zb      = ds['zb_mean'].values.squeeze()  # bed level
    zs      = ds['zs_mean'].values.squeeze()  # water level
    zsvar   = ds['zs_var'].values.squeeze()   # water level variance
    Hrms_IG = 2*np.sqrt(2)*np.sqrt(abs(zsvar)) 
    
    Hrms   = ds['H_mean'].values.squeeze()   # rms wave height
    if vegopt == True:
        veg    = ds['vegtype'].values.squeeze()[0,:]  # vegtype
    
    zs_grid = np.zeros_like(zb)
    zs_grid[:] = zs
    
    df_results = pd.DataFrame.from_dict({'x':x,
                                         'bed':zb,
                                         'wl':zs_grid,
                                         'hrms': Hrms,
                                         'hrms_ig':Hrms_IG})
    if vegopt == True:
        df_results['veg'] = veg
    
    # clip results once levee is reached
    if (df_results['wl']<=df_results['bed']).sum() >0:
        first_idx  = df_results[df_results['wl']<=df_results['bed']].index[0]
        df_results = df_results.iloc[0:first_idx+1,:]

    return df_results

def XB_plot_val(df, var = 'hrms', vegopt = True):
    fig, ax = plt.subplots()
    
    ax.plot(df.x, df.bed, color = 'saddlebrown')  # bed
    
    if var == 'hrms':
        
        ax.plot(df.x, df.wl + df.hrms, color = 'dodgerblue')
        ax.plot(df.x, df.wl , color = 'lightblue', lw = 1)
    else:
        ax.plot(df.x, df[var], color = 'dodgerblue')

    if vegopt == True:
        ax.plot(df.x.loc[df.veg==1], df.bed.loc[df.veg==1], color = 'limegreen')
        # add vertical line
        ax.axvline(df.x.loc[df.x.loc[df.veg==1].index[0]], color = 'limegreen', alpha = 0.5, ls = '--', lw = 0.5)
        ax.axvline(df.x.loc[df.x.loc[df.veg==1].index[-1]], color = 'limegreen', alpha = 0.5, ls = '--', lw = 0.5)
        
    ax.spines[['right', 'top']].set_visible(False)

    ax.set_ylabel('z [m]')
    ax.set_xlabel('Distance [m]')



def XB_plot_dif_val(df_veg, df_noveg, var = 'hrms', save_path = ''):
    plt.figure(figsize=(8,6))
    ax1 = plt.subplot2grid(shape=(3, 4), loc=(0, 0), colspan=3, rowspan = 2)
    ax2 = plt.subplot2grid(shape=(3, 4), loc=(0, 3), colspan=1, rowspan = 2)
    ax3 = plt.subplot2grid(shape=(3, 4), loc=(2, 0), colspan = 3, rowspan =1)                      
    # ax4 = plt.subplot2grid(shape=(3, 4), loc=(2, 2))                       
    # ax5 = plt.subplot2grid(shape=(3, 4), loc=(2, 3))
    axes = [ax1, ax2, ax3]
    
    
    # ax1 (top left plot)
    ax1.plot(df_veg.x, df_veg.bed, color = 'saddlebrown')  # bed
    ax1.plot(df_veg.x.loc[df_veg.veg==1], df_veg.bed.loc[df_veg.veg==1], color = 'limegreen')

    if var == 'hrms':
        ax1.plot(df_noveg.x, df_noveg.wl + df_noveg.hrms, color = 'dodgerblue')
        ax1.plot(df_noveg.x, df_noveg.wl , color = 'lightblue', lw = 1)
        
        ax1.plot(df_veg.x, df_veg.wl + df_veg.hrms, color = 'tomato', ls = '--', lw = 1.5)
        ax1.plot(df_veg.x, df_veg.wl , color = 'salmon', lw = 1, ls = '')        
        
        
    else:
        ax1.plot(df_noveg.x, df_noveg[var], color = 'dodgerblue')    
        ax1.plot(df_veg.x  , df_veg[var], color = 'tomato', ls = '--')

    # add vertical line
    ax1.axvline(df_veg.x.loc[df_veg.x.loc[df_veg.veg==1].index[0]], color = 'limegreen', alpha = 0.5, ls = '--', lw = 0.5)
    ax1.axvline(df_veg.x.loc[df_veg.x.loc[df_veg.veg==1].index[-1]], color = 'limegreen', alpha = 0.5, ls = '--', lw = 0.5)

    
    ax1.set_ylabel('z [m]')
    ax1.set_ylim(-1,)
    ax1.set_xticklabels([])
    
    # ax2  (top right plot)
    idx_veg_end = df_veg.loc[df_veg['veg'] == 1].index[-1]
    hrms_end_veg    = df_veg['hrms'].loc[idx_veg_end] 
    hrms_end_noveg  = df_noveg['hrms'].loc[idx_veg_end]
    
    ax2.bar([0,1],[hrms_end_veg, hrms_end_noveg], color = ['limegreen','saddlebrown'])
    ax2.set_ylabel('H$_{rms}$ [m]'); 
    ax2.set_xticks([0,1]); ax2.set_xticklabels(['veg','noveg'])
    ax2.set_xlabel('Wave height \nat end of \nvegetation zone')
    
    # ax3 (bottom plot)
    ax3.set_xlabel('Distance [m]')
    ax3.set_ylabel('H$_{rms}$ [m]')
    ax3.plot(df_veg.x, df_veg.hrms, color = 'limegreen')
    ax3.plot(df_noveg.x, df_noveg.hrms, color = 'saddlebrown')
    

    for ax in [ax1, ax2, ax3]:
        ax.spines[['right', 'top']].set_visible(False)
        
    savefig(save_path, dpi = 300,bbox_inches='tight')


def XBparams_sb(XB_templatefiledir,XB_outputfiledir, nxb, zs0):                                     #maken van param file for surfbeat [vvz]
    XBloadpathin= os.path.join(XB_templatefiledir, 'paramssb_0.txt');
    XBloadpathout= os.path.join(XB_outputfiledir, 'params.txt');
    fin = open(XBloadpathin, 'r')
    fout = open(XBloadpathout, 'w')
    endno = {}
    end = 0
    for i,line in enumerate(fin):
        end += len(line)
        endno[end] = i   
    fin.close()
    fin = open(XBloadpathin, 'r')
    for i in np.arange(endno[end]+1):
        line=fin.readline()
        print (line)
        if i==17:
            fout.write('nx           = ' + str(nxb) + '\n')
        elif i==35:
            fout.write('zs0          = ' + str(zs0) + '\n')
        else:             
            fout.write(line)
    fin.close()
    fout.close()

def XBjonswap(XB_templatefiledir, XB_outputfiledir,Hm0,fp):                                        #maken van param file Jonswap spectrum
    XBloadpathin= os.path.join(XB_templatefiledir, 'jonswap_0.txt');
    XBloadpathout= os.path.join(XB_outputfiledir, 'jonswap.txt');
    fin = open(XBloadpathin, 'r')
    fout = open(XBloadpathout, 'w')
    endno = {}
    end = 0
    for i,line in enumerate(fin):
        end += len(line)
        endno[end] = i   
    fin.close()
    fin = open(XBloadpathin, 'r')
    for i in np.arange(endno[end]+1):
        line=fin.readline()
        print (line)
        if i==0:
            fout.write('Hm0        =  ' + str(Hm0) + '\n')
        elif i==1:
            fout.write('fp         =  ' + str(fp) + '\n')
        else:             
            fout.write(line)
    fin.close()
    fout.close()



def XBwriteveg(XB_sim_dir, d_veg):
    XB_veggiefile  = os.path.join(XB_sim_dir, 'veggiefile.txt');         # hold veg name
    XB_vegtypefile = os.path.join(XB_sim_dir, d_veg['name_v'] + '.txt'); # hold veg parameters

    # write veggiefile
    fout = open(XB_veggiefile, 'w')
    fout.write(d_veg['name_v'] +  '.txt\n')
    fout.close()

    # write vegtypefile
    fout = open(XB_vegtypefile, 'w')
    fout.write('nsec = 1\n')
    fout.write('ah = ' + str(d_veg['hv']) + '\n')
    fout.write('bv = ' + str(d_veg['bv'])+ '\n')
    fout.write('N = ' + str(d_veg['nv'])+ '\n')
    fout.write('Cd = ' + str(d_veg['cd'])+ '\n')
    fout.close()



def write_xb_hydrodynamics(xgrid, main_dir, wd, zs0, Hm0, fp):
    nxb = len(xgrid)-1 # number of xgrid cells
    XB_templatefiledir = join(main_dir,'02_XB_sims/00_dummy_input_xbfiles')
    XB_outputfiledir   = wd
    
    # write params file
    XBparams_sb(XB_templatefiledir,XB_outputfiledir, nxb, zs0)
    # write wave related input
    XBjonswap(XB_templatefiledir, XB_outputfiledir,Hm0,fp)