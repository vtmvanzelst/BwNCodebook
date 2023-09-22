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
# import subprocess  # only needed if you would like to start the simulations via a python script. 



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
    ax.plot(df.x, df[var], color = 'dodgerblue')

    if vegopt == True:
        ax.plot(df.x.loc[df.veg==1], df.bed.loc[df.veg==1], color = 'limegreen')
        # add vertical line
        ax.axvline(df.x.loc[df.x.loc[df.veg==1].index[0]], color = 'limegreen', alpha = 0.5, ls = '--', lw = 0.5)
        ax.axvline(df.x.loc[df.x.loc[df.veg==1].index[-1]], color = 'limegreen', alpha = 0.5, ls = '--', lw = 0.5)
        
    ax.spines[['right', 'top']].set_visible(False)

    # ax.set_ylabel()
    # ax.set_xlabel()


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




#%%


def main(main_dir:str          = 'path\to\BwNCodebook',
         simulation_name:str   = 'name',
         selected_transect:int = 1,
         zs0: float            = 3.5,
         Hm0:float             = 1.25,
         fp:float              = 1/4.5):
    
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
    
    # working directory and create folder if not exist
    wd       = join(main_dir, '02_XB_sims/' + simulation_name)
    check_dir(wd)

    # load pickle file for selected transect. The file contains bed level information
    df_transect = from_pickle(join(main_dir, '01_transects' + '/transect_' + str(selected_transect)) + '.pkl')
    
    #     Write grid and depth to XBeach text files
    xgrid  = df_transect.dist_total.values; np.savetxt(join(wd,'xgrid.grd'), xgrid)
    dep    = df_transect.elevation.values; np.savetxt(join(wd,'bed.dep'), dep)
    
  
    
  
    
  
    
if __name__ == '__main__':
    main(main_dir          = main_dir, 
         simulation_name   = simulation_name, 
         selected_transect = selected_transect,
         zs0               = waterlevel,
         Hm0               = Hm0,
         fp                = fp
         )
    
  
    
    # USER INPUT BELOW
    # directory input
    main_dir             = r'c:\Users\zelst\OneDrive - Stichting Deltares\Documents\GitHub\BwNCodebook'
    simulation_name      = 'test01'
    
    # transect (1,2 or 3)
    selected_transect    = 1
    
    # hydrodynamic input
    zs0                  = 3.5        # design storm water level in meters (vertical reference frame according to the *.dep-file)
    Hm0                  = 1.25       # significant wave height
    fp                   = 1/4.5      # peak frequency (= 1 / Tp)
    
    





#=============================================================================
# 2. Create input for XBeach
#_____________________________________________________________________________
# -- A. grid and topography


#_____________________________________________________________________________
# -- B.  Vegetation in XBeach 
#    The XBeach vegetation module requires different veg input (files):    
#        1. the veggiemap file -> indicates where vegetation is present in the computational domain(0 = no veg, 1 = vegtype1, 2 = vegtype2, ..n)
#        2. the veggie file    -> indicates where WHAT type of vegetation is present by providing the vegtype.txt files that corresponds to the integers provided in the veggiepatch file. 
#        3. the vegtypes file(s) -> indicate the vegetation schematization. (nsec = vertical layer nr, N = density, d = representative diamter, h = height, Cd = drag coef.)
#        4. 'vegetation' = 1 keyword in the params file (set vegetation = 0 to turn off the vegetation module)
#
#     For detailed information visit: https://xbeach.readthedocs.io/en/latest/xbeach_manual.html#vegetation-input

# Example of writing vegetation presence to XBeach text file
df_transect.vegetation.loc[df_transect.vegetation == -1] = 0
np.savetxt(join(wd,'vegpatch.txt'), df_transect.vegetation.values)
    

# copy veggiefile and vegtype_dummy to working dir
shutil.copyfile(join(main_dir,'02_XB_sims/00_dummy_input_xbfiles/veggiefile.txt'),join(wd,'veggiefile.txt'))
shutil.copyfile(join(main_dir,'02_XB_sims/00_dummy_input_xbfiles/vegtype_dummy.txt'),join(wd,'vegtype_dummy.txt'))
# Manually make change to these files to provide the model with your own vegetation parameterization!


#_____________________________________________________________________________
# -- C. Hydrodynamic boundary conditions
# define input and output dir locations
nxb = len(xgrid)-1 # number of xgrid cells
XB_templatefiledir = join(main_dir,'02_XB_sims/00_dummy_input_xbfiles')
XB_outputfiledir   = wd

# write params file
XBparams_sb(XB_templatefiledir,XB_outputfiledir, nxb, zs0)
# write wave related input
XBjonswap(XB_templatefiledir, XB_outputfiledir,Hm0,fp)




#=============================================================================
# 3. Start simulation
# copy bat file to simulation folder
fn_bat = join(wd,'run_xb.bat')
shutil.copyfile(join(main_dir,'02_XB_sims/00_dummy_input_xbfiles/run_xb.bat'),fn_bat)


#%%
# Start the XB simulation by running: 'run_xb.bat'
log_xbeach = subprocess.run(fn_bat, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)




#%%

#=============================================================================
# 4. Load results
dir_nc     = join(wd,'xboutput.nc')
df_results = XB_load_results(dir_nc, vegopt = True)



#=============================================================================
# 5. Post-processing
# simple function to get you started..
XB_plot_val(df_results)








