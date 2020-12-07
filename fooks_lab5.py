#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import scipy
import rasterio
import numpy as np
import pandas as pd


# In[46]:


#Read the DEM as a numpy array
dem_raster = rasterio.open('bigElk_dem.tif')
dem_array = dem_raster.read(1)
cell_size = dem_raster.transform[0]

#Calculate slope and aspect
def slopeAspect(dem, cs):
    from math import pi
    from scipy import ndimage
    kernel = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]])
    dzdx = ndimage.convolve(dem, kernel, mode='mirror') / (8 * cs)
    dzdy = ndimage.convolve(dem, kernel.T, mode='mirror') / (8 * cs)
    global slp
    slp = np.arctan((dzdx ** 2 + dzdy ** 2) ** 0.5) * 180 / pi
    ang = np.arctan2(-dzdy, dzdx) * 180 / pi
    global aspect
    aspect = np.where(ang > 90, 450 - ang, 90 - ang)
    return slp, aspect
slp, asp = slopeAspect(dem_array, cell_size)

#Reclassify aspect grid into 8 cardinal directions
def reclassAspect(npArray):
    return np.where((npArray > 22.5) & (npArray <= 67.5), 2,
    np.where((npArray > 67.5) & (npArray <= 112.5), 3,
    np.where((npArray > 112.5) & (npArray <= 157.5), 4,
    np.where((npArray > 157.5) & (npArray <= 202.5), 5,
    np.where((npArray > 202.5) & (npArray <= 247.5), 6,
    np.where((npArray > 247.5) & (npArray <= 292.5), 7,
    np.where((npArray > 292.5) & (npArray <= 337.5), 8, 1)))))))
aspect = reclassAspect(asp)

#Reclassify slope grid into 10 classes
def reclassByHisto(npArray, bins):
    histo = np.histogram(npArray, bins)[1]
    rClss = np.zeros_like(npArray)
    for i in range(bins):
        rClss = np.where((npArray >= histo[i]) & (npArray <= histo[i + 1]),
                         i + 1, rClss)
    return rClss
slope = reclassByHisto(slp, 10)


# In[47]:


#Calculate NDVI for all years of analysis from the Landsat images and Calculate Recovery Ratio for each pixel for each year.
fire_raster = rasterio.open('fire_perimeter.tif')
fire_array = fire_raster.read(1)

subdir = 'L5_big_elk/'
files = os.listdir(subdir)
files = [file for file in files if file.endswith('.tif')]
b3list = [subdir+file for file in files if file.endswith('B3.tif')]
b4list = [subdir+file for file in files if file.endswith('B4.tif')]

recovery_ratios = []
for i in range(10):
    b3 = b3_raster.read(1)
    b4 = b4_raster.read(1)
    b4_raster = rasterio.open(b4list[i])
    b3_raster = rasterio.open(b3list[i])
    ndvi = (b4-b3) / (b4+b3)
    ndvi2 = ndvi[fire_array == 2].mean()
    rr = ndvi / ndvi2
    print(rr[fire_array == 1].mean())
    recovery_ratios.append(rr)


# In[48]:


#Calculate the trend of the RR for each pixel across the years of analysis
#Print the mean RR for each year
ndvi_slope = np.zeros_like(recovery_ratios[0])
xs = range(10)
for px in range(recovery_ratios[0].size):
    ys = [p[px] for p in recovery_ratios]
    ndvi_slope[px] = scipy.polyfit(xs, ys, 1)[0]
ndvi_slope = ndvi_slope.reshape(b3.shape)
print('mean rr:', ndvi_slope[fire_array == 1].mean())
ndvi_slope[fire_array != 1] = np.nan


# In[49]:


#Write a generic function that calculates Zonal Statistics as Table.
def zstats_astable(zone_raster, value_raster):
    values = np.unique(zone_raster)
    table = pd.DataFrame(columns = ['zone', 'mean', 'stddev', 'min', 'max','count',])
    for zone in table['zone']:
        inzone = value_raster[zone_raster == zone]
        table.at[table['zone'] == zone, 'mean'] = inZone.mean()
        table.at[table['zone'] == zone, 'stddev'] = inZone.std()
        table.at[table['zone'] == zone, 'min'] = inZone.min()
        table.at[table['zone'] == zone, 'max'] = inZone.max() 
        table.at[table['zone'] == zone, 'count'] = (zone_raster == zone).sum()
    return table


# In[50]:


#Calculate zonal stats of the coefficient of recovery for each slope and each aspect.
aspect_burn = np.where(fire_array == 1, aspect, np.nan)
slope_burn = np.where(fire_array == 1, slope, np.nan)
zstats_astable(ndvi_slope, aspect_burn).to_csv('aspect.csv')
zstats_astable(ndvi_slope, slope_burn).to_csv('slope.csv')


# In[ ]:


#Export the final coefficient of recovery for the  burned  areaas  a  GeoTiff
fire_profile = fire_raster.profile
fire_profile['nodata'] = -99
fire_profile['dtype'] = np.float32
with rasterio.open('cof_recovery.tif', 'w', **fire_profile) as ds:
    ds.write_band(1, ndvi_slope)


# In[ ]:


#Final conclusions
print("Lower elevations and east-facing aspects appear to have the best recvoery. ")


# In[ ]:




