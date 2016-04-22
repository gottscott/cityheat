#s functions
from osgeo import ogr, osr
import os
from shapely.geometry import Point
import shapely.geometry
import shapely.wkt
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader
import gdal

# reproject 'Parks_Dissolved.shp'
def mapmean(tempDF, meta, name = '', option = 0):
    import cartopy.crs as ccrs
    from cartopy.io.img_tiles import MapQuestOSM
    #fig  = plt.figure(figsize=(30, 30))
    x = meta['location:Longitude'].values
    y = meta['location:Latitude'].values
    c = tempDF[meta.sensornumber].mean()
    marker_size = 100
    fig = plt.figure(figsize=[15,15])
    imagery = MapQuestOSM()
    ax = plt.axes(projection=imagery.crs)
    
    ax.set_extent(( meta['location:Longitude'].min()-.005, 
                   meta['location:Longitude'].max()+.005 , 
                   meta['location:Latitude'].min()-.005,
                   meta['location:Latitude'].max()+.005))
    ax.add_image(imagery, 14)

    cmap = matplotlib.cm.OrRd
    bounds = np.linspace(round((c.mean()-3)),round((c.mean()+3)),13)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    plotHandle = ax.scatter(x,y,c = c, s = 350, transform=ccrs.Geodetic(), 
                 cmap = cmap,
                 norm = norm)
    if option ==0 : 
        cbar1 = plt.colorbar(plotHandle, label = 'Temperature in $^\circ $C')
    else : 
        cbar1 = plt.colorbar(plotHandle, label = option)
    
    lon = x[np.nanargmax(c)]
    lat = y[np.nanargmax(c)]
    at_x, at_y = ax.projection.transform_point(lon, lat,
                                               src_crs=ccrs.Geodetic())
    plt.annotate(
        '%2.1f'%np.nanmax(c.values), xy=(at_x, at_y), #xytext=(30, 20), textcoords='offset points',
        color='black', backgroundcolor='none', size=22,
        )

    lon = x[np.nanargmin(c)]
    lat = y[np.nanargmin(c)]
    at_x, at_y = ax.projection.transform_point(lon, lat,
                                               src_crs=ccrs.Geodetic())
    plt.annotate(
        '%2.1f'%np.nanmin(c.values), xy=(at_x, at_y), #xytext=(30, 20), textcoords='offset points',
        color='black', size = 22, backgroundcolor='none')

    plt.annotate(
        '$\mu = $ %2.1f'%np.nanmean(c.values), (0.01,0.01), xycoords ='axes fraction', #xytext=(30, 20), textcoords='offset points',
        color='black', size = 22, backgroundcolor='none')
    
    plt.title('Mean Temperature %s'%name)
    filename = './plots/meantempmap%s.eps'%name
    plt.savefig(filename, format = 'eps', dpi = 600)
    
def reproject_shapefile(fname,outfilename,outProjection = "WGS84"): 
    #fname = 'Parks_Dissolved.shp'
    driver = ogr.GetDriverByName('ESRI Shapefile')

    # input SpatialReference
    datasource = ogr.GetDriverByName('ESRI Shapefile').Open(fname)
    inlayer  = datasource.GetLayer()
    inSpatialRef = inlayer.GetSpatialRef()

    # output SpatialReference
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.SetWellKnownGeogCS(outProjection)

    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    # get the input layer
    inDataSet = driver.Open(fname)
    inLayer = inDataSet.GetLayer()

    # create the output layer
    outputShapefile = outfilename + '.shp'
    if os.path.exists(outputShapefile):
        driver.DeleteDataSource(outputShapefile)
    outDataSet = driver.CreateDataSource(outputShapefile)
    outLayer = outDataSet.CreateLayer("basemap_wgs84", geom_type=ogr.wkbMultiPolygon)

    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()

    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # destroy the features and get the next input feature
        outFeature.Destroy()
        inFeature.Destroy()
        inFeature = inLayer.GetNextFeature()
        
    spatialRef = osr.SpatialReference()
    spatialRef.SetWellKnownGeogCS(outProjection)
    spatialRef.MorphToESRI()
    file = open(outfilename+'.prj', 'w')
    file.write(spatialRef.ExportToWkt())
    file.close()

    # close the shapefiles
    inDataSet.Destroy()
    outDataSet.Destroy

def extract_raster_values(X,Y, rasterfile  ): 
    # Transform lat/lon values to the raster projection 
    
	sourceEPSG = 4326
	source = osr.SpatialReference()
	source.ImportFromEPSG(sourceEPSG)

	# Read in raster data to 	
    file = rasterfile #'exportImage'
    layer = gdal.Open(file)
    gt =layer.GetGeoTransform()
    bands = layer.RasterCount

    inProj = osr.SpatialReference()
    inProj.ImportFromWkt(layer.GetProjection())
	#target = osr.SpatialReference()
	#target =layer.GetGeoTransform()
	#target.ImportFromEPSG(calculationProjection)

	transform = osr.CoordinateTransformation(source, target)

    elevation = np.zeros(X.shape[0])

    i = 0 
    for x,y in zip(X,Y): 
        if ~np.isnan(x) & ~np.isnan(y): 
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(x,y)
            point.Transform(transform)
            #print point.ExportToWkt()

            x = point.GetPoints()[0][0]
            y = point.GetPoints()[0][1]

            rasterx = int((x - gt[0]) / gt[1])
            rastery = int((y - gt[3]) / gt[5])

            elevation[i] = layer.GetRasterBand(1).ReadAsArray(rasterx,rastery, 1, 1) 
            #print layer.GetRasterBand(1).ReadAsArray(rasterx,rastery, 1, 1) 
        else:
            print 'missing data at ', i
        i = i+1
    return elevation

def compute_distance_to_feature(X,Y,feature_file, feature_name = 'none', calculationProjection = 6347):
# compute distance from an array of lons/lats to a feature
# if multiple features in shapefile, specify feature_name
#feature_file = 'data/Parks_Dissolved_reproj.shp'
#feature_name = 'none'
#calculationProjection = 6347
#X = meta.drop(64, axis=0)['location:Longitude'][selected].values
#Y = meta.drop(64, axis=0)['location:Latitude'][selected].values
    
    #feature_file = 'data/Parks_Dissolved_reproj.shp'
    # Read in shapefile for the feature
    shapefile = ogr.Open(feature_file)
    layer = shapefile.GetLayer(0)
    
    # Select the correct feature
    if feature_name == 'none': 
        feature = layer.GetFeature(0)
        geometry = feature.GetGeometryRef()
    else: 
        for i in range(layer.GetFeatureCount()):
            feature = layer.GetFeature(i)
            name = feature.GetField('name')
            if name == feature_name : 
                geometry = feature.GetGeometryRef()
    
    ### Transform
    # reproject lat/lon values to the raster
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(4326) # lat lon

    outSpatialRef = osr.SpatialReference()
    #outSpatialRef = layer.GetSpatialRef()
    outSpatialRef.ImportFromEPSG(calculationProjection)
    
    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    
    # check that the shapefile is in the correct projection
    if geometry.GetSpatialReference() != outSpatialRef : 
        # reproject 
        coordTrans2 = osr.CoordinateTransformation(geometry.GetSpatialReference(), outSpatialRef)
        geometry.Transform(coordTrans2)
    
    shape = shapely.wkt.loads(geometry.ExportToWkt())
    
    i = 0
    distance_to_park = np.zeros(X.shape)
    for (x,y) in zip(X,Y): 
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x,y)
        point.Transform(coordTrans)
        distance_to_park[i] = Point(point.GetPoints()[0]).distance(shape)
        i = i+1
    #distance_to_park = distance_to_park*2*np.pi/360*6371000.    
    return distance_to_park

def map_data(data, lat,lon, shapefiles = 'none'): 
    # gis_layers a list of filenames
    
    # first, check that lat, lon, and data all same size
    if lat.shape != lon.shape : 
        print 'Lat and lon data shape do not match' 
    elif lat.shape != data.shape : 
        print 'Lat/lon & data shape do not match'
        
    # set up figures and axes    
    plt.figure(figsize=[15,15])
    if shapefiles == 'none': 
        imagery = MapQuestOSM()
        ax = plt.axes(projection=imagery.crs)
        ax.add_image(imagery, 14)
    else: 
        ax = plt.axes(projection=ccrs.PlateCarree()) # note to self:put in UTM 18/19

    ax.set_extent((lon.min()-.05, 
               lon.max()+.05 , 
               lat.min(),
               lat.max()+.05))
    
    # if there are shapefiles, read in and display
    if shapefiles != 'none' : 
        for file in shapefiles: 
            shape_feature = ShapelyFeature(Reader(file).geometries(),
                                ccrs.PlateCarree(), facecolor='none', edgecolor='black')
            ax.add_feature(shape_feature)
    
    # define colormap
    cmap = matplotlib.cm.OrRd
    bounds = np.linspace(round((data.mean()-2*data.std())),round((data.mean()+2*data.std())),13)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    # plot
    plotHandle = ax.scatter(lon,lat,c= data, s = 150, transform=ccrs.Geodetic(), 
                 cmap = cmap,
                 norm = norm)
    
    return ax, plotHandle

def band10_toLST(band10): 
    # convert landsat band 10 to land surface temperature
    # convert from digital numbers to (spectral) radiance
    m_l = 3.3420E-04 # gain
    a_l = 0.10000 # bias
    radiance = m_l * band10 + a_l
    # convert from radiance to temperature by inverting the planck function
    k1 = 774.8853 
    k2 = 1321.0789
    temp = k2/ (np.log(k1/radiance +1)) -273.15
    temp[temp < -100] = 'nan'
    
    return temp

def albedo(B1,B2,B3,B4,B5) : 
# ((0.356*B1) + (0.130*B2) + (0.373*B3) + (0.085*B4) + (0.072*B5) -0.018) / 1.016
    return ((0.356*B1) + (0.130*B2) + (0.373*B3) + (0.085*B4) + (0.072*B5) -0.018)/1.016

#reproject_shapefile('data/Parks_Dissolved.shp',outfilename='data/Parks_Dissolved_reproj')
