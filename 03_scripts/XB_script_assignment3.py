# -*- coding: utf-8 -*-
"""

Script to set-up, run and post-process 1D XB simulation as part of the BwN coarse 2023.


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


# local funcs
from support_funcs.XB_func import XB_plot_val, XB_plot_dif_val, XBparams_sb, XB_load_results, write_xb_hydrodynamics, XBwriteveg

def to_pickle(fn_pickle, df):
	file = open(fn_pickle,'wb')
	pickle.dump(df,file)
	file.close()
	
def from_pickle(fn_pickle):
    file = open(fn_pickle, 'rb')
    df = pickle.load(file)
    file.close()      
    return df


def check_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def check_dir_plot(directory):
    if not os.path.exists(directory):
        print('ERROR: output is not avilable in the provided folder..')
        opt = False, 
    else:
        opt = True
    return opt



def main(main_dir:str          = 'path\to\BwNCodebook',
         simulation_name:str   = 'name',
         selected_transect:int = 1,
         zs0: float            = 3.5,
         Hm0:float             = 1.25,
         fp:float              = 1/4.5,
         d_veg:dict            = {},
         opt_prepare_and_run_xb= True,
         opt_postprocessing    = True):
    
    """
     Description:
     ------------
     This function loads data and writes XBeach input data to file for one of the three provided transects
    
     Parameters:
     ------------
     main_dir: String
         directory of the BwNCodeBook checkout
     simulation_name : String
         user defined name of the simulation (a folder will be created if not exist)
     selected_transect: Integer
         number of the transect
     
     zs0 : Float
         design storm water level in meters (vertical reference frame according to the *.dep-file)
     Hm0: integer
         significant wave height
    
     fp: Float
         peak frequency (= 1 / Tp)
    
     Returns:
     ------------
     """
    
    # directories
    wd       = join(main_dir, '02_XB_sims/' + simulation_name)
    dir_noveg = join(wd,'noveg'); 
    dir_veg   = join(wd,'veg'); 
    dirs      = [dir_noveg, dir_veg]


    # PREPARE AND RUN XB
    if opt_prepare_and_run_xb == True:   
        # in the working directory we create two folders (veg and noveg)
        check_dir(wd);check_dir(dir_noveg);check_dir(dir_veg)
        # load pickle file for selected transect. The file contains bed level information, and information on vegetation extent
        df_transect = from_pickle(join(main_dir, '01_transects' + '/transect_' + str(selected_transect)) + '.pkl')
        
        #     Write grid and depth to XBeach text files
        for wd_sel in dirs:
            xgrid  = df_transect.dist_total.values; np.savetxt(join(wd_sel,'xgrid.grd'), xgrid)
            dep    = df_transect.elevation.values; np.savetxt(join(wd_sel,'bed.dep'), dep)
        
            # hydrodynamics
            write_xb_hydrodynamics(xgrid, main_dir, wd_sel, zs0, Hm0, fp)
        
        
        #writing vegetation presence to XBeach text file
        # veg
        df_transect.vegetation.loc[df_transect.vegetation == -1] = 0
        np.savetxt(join(dir_veg,'vegpatch.txt'), df_transect.vegetation.values)
      
        # noveg
        df_transect['noveg'] = 0
        np.savetxt(join(dir_noveg,'vegpatch.txt'), df_transect.noveg.values)

      
        # -- B.  Vegetation in XBeach 
        #    The XBeach vegetation module requires different veg input (files):    
        #        1. the veggiemap file -> indicates where vegetation is present in the computational domain(0 = no veg, 1 = vegtype1, 2 = vegtype2, ..n)
        #        2. the veggie file    -> indicates where WHAT type of vegetation is present by providing the vegtype.txt files that corresponds to the integers provided in the veggiepatch file. 
        #        3. the vegtypes file(s) -> indicate the vegetation schematization. (nsec = vertical layer nr, N = density, d = representative diamter, h = height, Cd = drag coef.)
        #        4. 'vegetation' = 1 keyword in the params file (set vegetation = 0 to turn off the vegetation module)
        #
        #     For detailed information visit: https://xbeach.readthedocs.io/en/latest/xbeach_manual.html#vegetation-input
        
        # writing vegetation characteristics to file
        for wd_sel in dirs:
            XBwriteveg(wd_sel, d_veg)
      
        # execute simulation
        for wd_sel in dirs:
            shutil.copyfile(join(main_dir,'02_XB_sims/00_dummy_input_xbfiles/run_xb.bat'),join(wd_sel,'run_xb.bat'))             # copy bat file to simulation folder
        subprocess.run(join(dir_veg,'run_xb.bat'), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)     # Start the XB simulation by running: 'run_xb.bat'
        subprocess.run(join(dir_noveg,'run_xb.bat'), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) 
        
  
    # POST-PROCESSING
    if opt_postprocessing == True:
        opt = check_dir_plot(join(dir_veg, 'xboutput.nc'))
        
        if opt == True:        
            df_results_veg   = XB_load_results(join(dir_veg,'xboutput.nc'), vegopt = True)
            df_results_noveg = XB_load_results(join(dir_noveg,'xboutput.nc'), vegopt = False)
        
            # plotting (simple function to get you started..)
            XB_plot_val(df_results_veg, vegopt = True)
            XB_plot_val(df_results_noveg, vegopt = False)
            XB_plot_dif_val(df_results_veg, df_results_noveg, save_path = join(wd, 'veg_versus_noveg.png'))
            
    
#%%
    
    
if __name__ == '__main__':

    # USER INPUT 
    # =========================================================================
    # directory input
    main_dir             = r'c:\Users\zelst\OneDrive - Stichting Deltares\Documents\GitHub\BwNCodebook'
    simulation_name      = 'Transect3_test'
    
    # transect (1,2 or 3)
    selected_transect    = 3
    
    # hydrodynamic input
    zs0                  = 3.5        # design storm water level in meters (vertical reference frame according to the *.dep-file)
    Hm0                  = 1.25       # significant wave height
    fp                   = 1/4.5      # peak frequency (= 1 / Tp)
    
    # vegetation input
    name_v  = 'saltmarsh1'            # vegetation name (for your own administration) 
    hv      = 0.3                     # vegetation height  [m]
    bv      = 0.005                   # representative diameter [m]
    nv      = 100                     # vegetation density [1/m2]
    cd      = 0.5                     # drag coeficient [-]
    
    # options to run and/or postprocess
    opt_prepare_and_run_xb = True
    opt_postprocessing     = True
    # =========================================================================
    # END USER INPUT
    
    
    
    # no changes required below
    d_veg = {'name_v':name_v, 'hv':hv,
             'bv':bv,'nv':nv,'cd':cd}
    
    main(main_dir          = main_dir, 
         simulation_name   = simulation_name, 
         selected_transect = selected_transect,
         zs0               = zs0,
         Hm0               = Hm0,
         fp                = fp,
         d_veg             = d_veg, 
         opt_prepare_and_run_xb = opt_prepare_and_run_xb,
         opt_postprocessing     = opt_postprocessing)













