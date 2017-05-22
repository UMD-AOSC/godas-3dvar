#!/usr/bin/env python3
import argparse
import netCDF4 as nc
import numpy as np
from scipy.spatial import KDTree
import math

parser = argparse.ArgumentParser()
args=parser.parse_args()

args.grid_file="/export/cpc-lw-tsluka/travis.sluka/proj/godas-3dvar-mom6/test/INPUT/ocean_geometry.nc"



ncd = nc.Dataset(args.grid_file, 'r')
lons = ncd.variables['geolon'][:]
lats = ncd.variables['geolat'][:]
lon_1d = ncd.variables['lonq'][:]
lat_1d = ncd.variables['latq'][:]
lons *= math.pi/180.0
lats *= math.pi/180.0
mask = ncd.variables['wet'][:]
print(mask.shape)
ncd.close()

mask_f = mask.flatten()



R=6371000.0
def ll2xyz(lat,lon):
    ar = np.zeros( (3,) )
    ar[0] = R * math.cos(lat) * math.cos(lon)
    ar[1] = R * math.cos(lat) * math.sin(lon)
    ar[2] = R * math.sin(lat)
    return ar

print("Creating KDTree...")
lons_l = lons.flatten()[mask_f == 0.0]
lats_l = lats.flatten()[mask_f == 0.0]
coords=np.zeros((len(lons_l),3))
for i in range(len(lons_l)):
    coords[i,:] = ll2xyz(lats_l[i], lons_l[i])
kd = KDTree(coords)



print("calculating distances to coast...")
coast_dist = np.zeros(mask.shape)
for y in range(mask.shape[0]):
    print("{:05.2f}% complete".format(y/mask.shape[0]*100.0))
    for x in range(mask.shape[1]):
        if mask[y,x] == 0.0:
            coast_dist[y,x] = 0.0
        else:            
            c = ll2xyz(lats[y,x],lons[y,x])
            res = kd.query(c , 1)
            coast_dist[y,x] = res[0]



ncd = nc.Dataset('coast_dist.nc', 'w')
lat_dim = ncd.createDimension("lat", mask.shape[0])
lon_dim = ncd.createDimension("lon", mask.shape[1])

var = ncd.createVariable("lat", "f4", ("lat",))
var[:] = lat_1d
var.units="degrees_N"

var = ncd.createVariable("lon", "f4", ("lon",))
var[:] = lon_1d
var.units="degrees_E"

var = ncd.createVariable("coast_dist", "f4", ("lat","lon"))
var[:] = coast_dist
ncd.close()

