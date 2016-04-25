# Script to perform spatial analysis of correlationa and regression for mean temperature data and covariates pulled from shapefils

# import libraries
import numpy as np
import pandas as pd
from osgeo import ogr, osr
from shapely.geometry import Point
import shapely.geometry
import shapely.wkt
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader
import gdal
import spatialfunctions 

# import metadata
meta = pd.DataFrame(pd.read_csv('./data/falldownload/TempSensorFinal_results-4.csv', sep = ','))
# select sensors used for analysis
ebaltsensorsi = np.where(meta['location:Longitude']>= -76.6125831)#-76.61)# -76.6072591 )
ebaltsensors = meta.sensornumber.iloc[ebaltsensorsi]
selected = ebaltsensors

# Extract spatial variates from various files 
lon = meta.drop(64, axis=0)['location:Longitude'][selected].values
lat = meta.drop(64, axis=0)['location:Latitude'][selected].values

feature_file = 'data/Parks_Dissolved_reproj.shp'
distance_to_park = spatialfunctions.compute_distance_to_feature(lon,lat,feature_file, feature_name = 'none')

feature_file = 'data/Water-2/geo_export_a5e512a2-2fa2-47e8-ac32-78f85e798a3c.shp'
distance_to_water = spatialfunctions.compute_distance_to_feature(lon,lat,feature_file, feature_name = 'Inner Harbor')

rasterfile = 'data/BaltimoreDEM'
elevation = spatialfunctions.extract_raster_values(lon,lat, rasterfile)

rasterfile ='../satellite/LC80150332015229LGN00_B10.TIF' 
band10 = spatialfunctions.extract_raster_values(lon,lat, rasterfile)
LST = spatialfunctions.band10_toLST(band10)

rasterfile = 'data/TreeBaltimore_CanopyCover.img'
tree_cover = spatialfunctions.extract_raster_values(lon,lat, rasterfile)

path = '../satellite/data/BaltimoreLandsatSummer2015/L8 OLI_TIRS/LC80150332015229LGN00/'
B = np.zeros((5,lon.shape[0]))
for i in (1,2,3,4,5): 
    file = path+ 'LC80150332015229LGN00_B'+'%s'%i +'.TIF'
    B[i-1, :] = spatialfunctions.extract_raster_values(lon,lat,file)

alb = spatialfunctions.albedo(B[0,:], B[1,:], B[2,:], B[3,:], B[4,:], )
