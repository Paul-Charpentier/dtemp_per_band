# Read Me
In this repository is a small python script that compute band per band dTemp time series

## Dependencies 
Not so much, just make sure you have all python packages that are on the 10 first lines but they are very common you might probably have them already and otherwise just install them using `pip`

## Get Started 

Uneless you know what you are doing, the only lines you might modify are lines 17 to 30. 

First define the path to the lbl reposotory :
```
cwd = "/Volumes/dataDRS/spirou/lbl"
```
This directory should be organized as follow : 
```
. 
├── lblrdb/
├── lblrv/
│   └── TARGET_TEMPLATE/ 
│       ├── 1234567o_pp_e2dsff_tcorr_AB_TARGET_TEMPLATE_lbl.fits  #all your individual lblrv files
│       ├── ...
│       └── ...
|   └── ..._.../
│       ├── ...
│       └── ...
├── models/
├── recon/
├── template/
...
```
Then define the directory where you want the output files to be stored.
Put it wherever you want. 
```
outdir = change/me # Directory of output files
```

Then define which dtemp vectors that will be computed. 
Here I chose to compute them all. 
```
temps = [3000,3500,4000,4500,5000,5500,6000]
```

Then define the wavelength range of your bands, you can edit it or even add custom ones. 
```
bands = {
    "Y": [980,1080],
    "J": [1100,1400],
    "H": [1450,1800],
    "K": [2000,2500]
}
```

Then enter the list of targets you want this script to run on, and the list of corresponding template (usualy the same name as the target) : 
```
targets = ['AUMIC', 'EV_LAC', 'V1298TAU']
templates = targets
```

And that's it ! Now just run `python3 dtemp_per_band.py`. 
It should save all the computed data into `XXXX_DTEMP_bands.h5` files inside the output directory. 

## How to open and use the computed data ? 
Output files are `.h5` files. 
You can open them using : 
```
with h5py.File("PATH/TO/OUTPUT/FILES/AUMIC_DTEMP_bands.h5","r") as h5:

    dtemp = h5["dtemp"][:]
    sdtemp = h5["sdtemp"][:]
    time = h5["time"][:]
    bands = h5["bands"][:].astype(str)
    temps = h5["temperatures"][:]
```

Then you can use masks to use only the data you need. 
for example, if you want to plot the dTemp data from the H band computed using the 3500 K vector, use
```
ib = np.where(bands == "H")[0][0]
it = np.where(temps == 3500)[0][0]

import matplotlib.pyplot as plt
plt.errorbar(time, dtemp[:,ib,it], yerr=sdtemp[:,ib,it], fmt='ro')
plt.show()
```
