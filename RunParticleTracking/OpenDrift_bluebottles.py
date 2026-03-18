#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 15:22:57 2024

@author: Gaby Mayorga-Adame
"""

import pandas as pd
import numpy as np
import xarray as xr
from datetime import datetime
from datetime import timedelta
import matplotlib.pyplot as plt
from glob import glob

from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_global_landmask, reader_netCDF_CF_generic #, reader_netCDF_CF_unstructured # for corners



# NCI Paths 
BRAN_path = '/srv/scratch/oceanopen/BRAN2020_files/'
#BRAN_path = '/scratch/fu5/ca5240/BRAN_katana/'
#BRAN_path = '/scratch/fu5/ca5240/BRANsurf/'
#BRAN_path = '/g/data/gb6/BRAN/BRAN2020/daily/'
#BRAN_path = '/scratch/fu5/ca5240/BlueBlottles/'
#ERA5_path = '/g/data/rt52/era5/single-levels/reanalysis/' #10u/2023/
#ERA5_path = '/scratch/fu5/ca5240/ERA/'
ERA5_path = '/srv/scratch/oceanopen/ERA/'

#Forcing files
years = ['2012','2013']
ERA_ufiles = sorted([f for year in years for f in glob(ERA5_path + '10v_era5_oper_sfc_' + year + '*.nc')])
ERA_vfiles = sorted([f for year in years for f in glob(ERA5_path + '10u_era5_oper_sfc_' + year + '*.nc')])

#BRAN_ufiles = sorted([f for year in years for f in glob(BRAN_path + 'surf_ocean_u_' + year + '_*_opendrift.nc')])
#BRAN_vfiles = sorted([f for year in years for f in glob(BRAN_path + 'surf_ocean_v_' + year + '_*_opendrift.nc')])

# Reader Ocean Current
reader_curr_BRAN_u = reader_netCDF_CF_generic.Reader(BRAN_path + 'ocean_u_merged_' + '2010_2019_surf2m5_opendrift.nc')
reader_curr_BRAN_v = reader_netCDF_CF_generic.Reader(BRAN_path + 'ocean_v_merged_' + '2010_2019_surf2m5_opendrift.nc')

# Reader Wind
reader_wind_ERA5_u10 = reader_netCDF_CF_generic.Reader(ERA_ufiles)
reader_wind_ERA5_v10 = reader_netCDF_CF_generic.Reader(ERA_vfiles)


# Check dates
print('ERA5 start time:' + str(reader_wind_ERA5_u10.start_time))
print('ERA5 end time:' + str(reader_wind_ERA5_u10.end_time))

print('BRAN start time:' + str(reader_curr_BRAN_u.start_time))
print('BRAN end time:' + str(reader_curr_BRAN_u.end_time))

# random % of wind drift
mu, sigma = 0.017, 0.0013 # mean and standard deviation
np.random.seed(0)
s = np.random.normal(mu, sigma, 1000)

# Run right handed / left-handed
from BBDrift_AS2024 import BBDrift

# using the Bluebottle parametrisation
o = BBDrift(loglevel=20)
# Adding readers
o.add_reader([reader_curr_BRAN_u, reader_curr_BRAN_v,reader_wind_ERA5_u10,reader_wind_ERA5_v10])

# Config
o.set_config('environment:fallback:x_sea_water_velocity', None) # so that the particles are deactivated when outide of ROMS domain
#o.set_config('environment:fallback:x_wind', None)  # so that the particles are deactivated when outide of BARRA2 domain
o.set_config('drift:horizontal_diffusivity', 10) # not for Luca: adding random walk diffusion for unresolved variability
o.set_config('general:coastline_action', 'previous') # keep going after beaching

# SEED
# Orientation =1 is right-handed, -1 is left-handed
# To change the drift angle to a constant, e.g. 40 degree, add: Angle_course_fromW = np.zeros(N_particles)+30
N_particles = 1000
o.seed_elements(lon = -6.25268, lat= 36.45809, radius = 1000, number = N_particles,
                time = datetime(2013, 4, 2, 0, 0),
                Wind_perc = s,
                Orientation = -np.ones(N_particles))
o.description = 'left-handed backtrack'

# RUN
o.run(time_step=-3600, duration=timedelta(days=365),time_step_output=86400, outfile='2s_1yr_left.nc')
