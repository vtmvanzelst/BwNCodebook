![](https://github.com/vtmvanzelst/BwNCodebook/blob/main/imgs/saltmarsh.PNG)

# BwN Codebook
 Supporting software for Delft University of Technology Building with Nature course.


## Usage

For the design assignment we are making use of the hydrodynamic numerical model XBeach in combination with some Python scripts. 
These scripts provide some basic functionality. Feel challenged to adapt/develop these scripts to retrieve better system understanding. 
Below we provide instructions on:
1. Download and installation of Python Enviroment
2. Downland and use of XBeach
3. Provided Python scripts


## Python Environment and downloading assignment files
Make a clone of this repository using GitHub Desktop and create a Python environment using a *.yaml-file to prevent dependencies issues.

### 1. GitHub Desktop client

By these steps, the files that are hosted at GitHub are "pulled" to your machine:
* 1.1 Download GitHub Desktop via (https://desktop.github.com/) and install on your computer.
* 1.2 Browse to (https://github.com/vtmvanzelst/BwNCodebook), click on the green "Code"
   		button and select "Open with GitHub Desktop", or paste the URL into the GitHub client.

### 2. Mamba package manager

If you're not familiar with managing Python environments, please have a look at this
[introduction](https://earth-env-data-science.github.io/lectures/environment/python_environments.html?highlight=conda)
first. The bottom line is that it is good practice to manage your software environments
to avoid dependency conflicts.

Continue with installing a package manager (here we use Mambaforge):

#### Windows
1. Download the [Mambaforge](https://mamba.readthedocs.io/en/latest/installation.html) executable file for Windows from [Miniforge GitHub page](https://github.com/conda-forge/miniforge#mambaforge). On that page there are also binaries for Mac and Linux; and for `conda` package managers, so make sure you download
the `mambaforge` executable file for Windows. Install the executable by clicking on it; you can stay with the default settings by just clicking next through the installation client.

#### Mac and Linux
1. We recommend to install Mambaforge on Linux and Mac using a terminal. On Mac, you can
	   open a terminal by searching for "terminal" or "iterm". On Linux the hotkey to open a
	   terminal is "cntrl + shift + t". The commands to
	   install the package manager are copied from their documentation and can be run by
	   copying the commands below over to your terminal and pressing enter:  
	   ```bash
	   curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
	   bash Mambaforge-$(uname)-$(uname -m).sh
	   ```
2. Accept the user agreements, and allow the installation script to edit your profile
	   file. The profile file (`~/.bashrc` on Linux or possibly `~/.zshrc` on Mac) is the
	   first script which is being executed when you open a new terminal. The installation
	   script will add a few lines to that file to make the `mamba` command available every
	   time open a new terminal. 
3. Close the terminal.


### 3. Create Python Environment
Now that we have our package manager we will create the software
environments. 

1. Now that mambaforge is available on your machine, open a terminal. On Windows you
   should open the Miniforge prompt, which you can find by searching for it in the Start
   window. On Mac you can open a terminal by searching for "terminal" or "iterm". For
   Linux it's "cntrl + shift + t". 
2. You can check if mamba was installed by running the following command in the terminal: 
   ```bash 
   mamba --version
   ```` 
   It should output something like: 

   ```console
   ~ base ? mamba --version
   mamba 1.1.0
   conda 22.9.0
   ```
3. Now that mambaforge is installed, navigate in the terminal to the directory
   where you cloned the GitHub CoastalCodeBook repository. You can navigate the terminal
   using `cd`, which stands for "change directory". 
   - **Windows**: if you are on Windows and you installed the GitHub client using their default settings you can
   simply run `cd %userprofile%\Documents\GitHub\BwNCodeBook`. 
   - **Linux/Mac**: change to the directory where you cloned the GitHub repository. This
     will be something like `cd ~/path/to/github/repository`. 
4. The BwNCodebook root directory contains a yaml file:, (that describe the software
   dependencies). 

   ```
   Create the BwN environment by running: 
   ``` 
   mamba env create -f environment-coastal.yml
   ``` 
   If you want to give other names to
   the environments you can do so by adding a `-n` or `--name` flag, i.e., 
   ```bash
   mamba env create -n type_here_your_own_name -f environment-jupyterlab.yml
   ``` 


## XBeach
XBeach is a two-dimensional model for wave propagation, long waves and mean flow, sediment transport and morphological change on the coastal nearshore area (Roelvink et al., 2009). 
It includes a vegetation module that solves short-wave and long-wave -vegetation interaction (Van Rooijen et al., 2015; Van Rooijen et al., 2016). 

* XBeach software can be downloaded at: 	(https://download.deltares.nl/xbeach)

* XBeach manual:				(https://xbeach.readthedocs.io/en/latest/xbeach_manual.html)



## Python scripts




## References
1. Roelvink, D. et al. Modelling storm impacts on beaches, dunes and barrier islands. Coast. Eng. 56, 1133–1152 (2009).
2. van Rooijen, A. A. et al. Modeling of wave attenuation by vegetation with XBeach. E-proceedings 36th IAHR World Congr. 7 (2015).
3. Van Rooijen, A. A. et al. Modeling the effect of wave-vegetation interaction on wave setup. J. Geophys. Res. Ocean. 121, 4341–4359 (2016).
