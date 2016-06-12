# Plotting utilities for ibutton observation data 
# Anna Scott
# September 18, 2015 

import numpy as np 
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd
import matplotlib 
import matplotlib.pylab as pylab
#from mpl_toolkits.basemap import Basemap

def histPlot(tempDF, title='Histogram of Temperature'): 
    fig  = plt.figure(figsize=(12, 6))
    data = tempDF.values[~np.isnan(tempDF.values)] 
    n,bins, patches= plt.hist(data,20)
    #    n, bins, patches = plt.hist(tempDF.values[:].reshape([tempDF.shape[0]*tempDF.shape[1]]),10) # plots each sensor as a separate bar
    plt.axvline(data.mean(), 
                linestyle='dashed', 
                color = pd.tools.plotting._get_standard_colors(3)[2], 
                linewidth=2, 
                label = 'mean, %2.1f $^\circ$ C'%np.nanmean(data),
                alpha = 1.0)
    plt.axvline(data.mean()+data.std(),
                color = pd.tools.plotting._get_standard_colors(3)[1], 
                linestyle='dashed', 
                linewidth=2, 
                label = '$\pm \sigma$, %2.1f $^\circ$ C'%np.nanstd(data)
                )

    plt.axvline(data.mean()-data.std(), 
                color = pd.tools.plotting._get_standard_colors(3)[1], 
                linestyle='dashed', 
                linewidth=2, 
                )
    plt.xlabel('Temperature in $^\circ$ C')
    plt.ylabel('Count')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) #, ncol = 3)
    plt.title(title)
    filename = './plots/histogram' + title.replace(" ", "") + ".eps"
    plt.savefig(filename, format = 'eps', dpi = 600, )
    return  n, bins

# function draws plots of diurnal data given a pandas dataframe that has an 'hour' column added to it 
# with optional 2nd category to sort by 

def timeseriesplots(tempDF,meta, sorttype, sorttype2=0, option2=0):
        import cartopy.crs as ccrs
	from cartopy.io.img_tiles import MapQuestOSM
	# define dictionaries of options and titles
        options = {'landcoverclass': ['impervious', 'grass', 'soil'],#'dirt'],
                   'sunorshade': ['sun', 'partial', 'shade'],
                   'attachment': ['metal', 'wood', 'tree']}

        titles = {'landcoverclass': 'Land Cover Class',
                   'sunorshade': 'Shadiness',
                   'attachment': 'Attachment'}
        fig  = plt.figure(figsize=(30, 12))
        plt.subplot(1,2,2)
        imagery = MapQuestOSM()
	ax = plt.axes(projection=imagery.crs)
    
	ax.set_extent(( meta['location:Longitude'].min()-.005, 
                   meta['location:Longitude'].max()+.005 , 
                   meta['location:Latitude'].min()-.005,
                   meta['location:Latitude'].max()+.005))
	ax.add_image(imagery, 14) 
        # loop over every curve to draw 
        i = 0

        
        # loop over every curve to draw 
        for option in options[sorttype]:
            # if there is only one category to find 
            if sorttype2 == 0 : 
                index = meta['sensornumber'].values[np.where(meta[sorttype]==option)].astype(int)
                title = 'Temperature by %s'%titles[sorttype] 
            # if there are two categories to find
            else : 
                index = meta['sensornumber'].values[np.where( (meta[sorttype]==option) & (meta[sorttype2]==option2) )].astype(int)
                #print index.shape # use for debugging
                title = 'Temperature by %s for %s'%(titles[sorttype], option2 )
            # make sure that there was data found satisfying the criteria 
            if index.shape[0] >0 :
                lab = '%s, %s/%s sensors (%2.1f %%)'%(option, index.shape[0], meta.sensornumber.shape[0], index.shape[0]/float(meta.sensornumber.shape[0])*100 )                
                plt.subplot(1,2,1)
                plt.plot(tempDF.index, tempDF[index].mean(axis=1), linewidth = 3, label = lab)
                
                plt.subplot(1,2,2)
                #print 'now plotting map for %s'%option
                x = meta['location:Longitude'][index].values
                y = meta['location:Latitude'][index].values
                marker_size = 150
                for j in index:
                        plt.scatter(meta['location:Longitude'][j], meta['location:Latitude'][j], 
                                  s = marker_size, 
                          c = pd.tools.plotting._get_standard_colors(3)[i], 
                          label = option if j == index[0] else "", )
                #plt.hold(True)
                i = i+1


            else : 
                print 'skipping plot %s'%title
            #annotate maximum and minimum
        plt.subplot(1,2,1)
        plt.title(title)
        plt.ylabel('Temperature in $^\circ $C')
        lgd = plt.legend()#bbox_to_anchor = (.5, -.28),
                    #loc=8,
                   #borderaxespad=0.)
        plt.ylim([-6,6])
        matplotlib.rcParams.update({'font.size': 22})
        
        plt.subplot(1,2,2)
        plt.legend(#[options[sorttype]], 
                   bbox_to_anchor = (.5, -.28),
                    loc=8,
                   borderaxespad=0.)        
        matplotlib.rcParams.update({'font.size': 22})
        filename = './plots/timeseries%s%s.eps'%(option,option2)
        plt.savefig(filename, format = 'eps', dpi = 600, 
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.show()

def mapmean(tempDF, meta, name = '', option = 0): 
    import cartopy.crs as ccrs
    from cartopy.io.img_tiles import MapQuestOSM
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    #fig  = plt.figure(figsize=(30, 30))
    x = meta['location:Longitude'].values
    y = meta['location:Latitude'].values
    c = tempDF[meta.index].mean()
    marker_size = 350 
    imagery = MapQuestOSM()
    fig = plt.figure(figsize=[15,15])
    ax = plt.axes(projection=imagery.crs)
    
    ax.set_extent(( meta['location:Longitude'].min()-.005, 
                   meta['location:Longitude'].max()+.005 , 
                   meta['location:Latitude'].min()-.005,
                   meta['location:Latitude'].max()+.005))
    ax.add_image(imagery, 14)

    cmap = matplotlib.cm.OrRd
    bounds = np.linspace(round((c.mean()-3)),round((c.mean()+3)),13)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    plotHandle = ax.scatter(x,y,c = c, s = marker_size, transform=ccrs.Geodetic(), 
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
        '$\mu = $ %2.1f, $\sigma = $ %2.1f'%(np.nanmean(c.values), np.nanstd(c.values)), (0.01,0.01), xycoords ='axes fraction', #xytext=(30, 20), textcoords='offset points',
        color='black', size = 22, backgroundcolor='none')
    
    plt.title('Mean Temperature %s'%name)
    filename = './plots/meantempmap%s.eps'%name
    plt.savefig(filename, format = 'eps', dpi = 600)

# function draws plots of diurnal data given a pandas dataframe that has an 'hour' column added to it 
# with optional 2nd category to sort by 
def diurnalplots(tempDF, meta, sorttype, sorttype2=0, option2=0):
        # define dictionaries of options and titles
        import cartopy.crs as ccrs
	from cartopy.io.img_tiles import MapQuestOSM
        options = {'landcoverclass': ['impervious', 'grass', 'soil'],#'dirt'],
                   'sunorshade': ['sun', 'partial', 'shade'],
                   'attachment': ['metal', 'wood', 'tree']}

        titles = {'landcoverclass': 'Land Cover Class',
                   'sunorshade': 'Shadiness',
                   'attachment': 'Attachment'}

        fig  = plt.figure(figsize=(30, 12))
        plt.subplot(1,2,2)
        matplotlib.rcParams.update({'font.size': 22})
        
        imagery = MapQuestOSM()
        ax = plt.axes(projection=imagery.crs)
    
        ax.set_extent(( meta['location:Longitude'].min()-.005, 
                   meta['location:Longitude'].max()+.005 , 
                   meta['location:Latitude'].min()-.005,
                   meta['location:Latitude'].max()+.005))
        ax.add_image(imagery, 14)
        # loop over every curve to draw 
        i = 0
        for option in options[sorttype]:
            # if there is only one category to find 
            if sorttype2 == 0 : 
                index = meta['sensornumber'].values[np.where(meta[sorttype]==option)].astype(int)
                title = 'Diurnal Cycle by %s'%titles[sorttype] 
            # if there are two categories to find
            else : 
                index = meta['sensornumber'].values[np.where( (meta[sorttype]==option) & (meta[sorttype2]==option2) )].astype(int)
                #print index.shape # use for debugging
                title = 'Diurnal Cycle by %s for %s'%(titles[sorttype], option2 )
                
            # make sure that there was data found satisfying the criteria 
            if index.shape[0] >0 :
                lab = '%s, %s/%s sensors (%2.1f %%)'%(option, index.shape[0], meta.sensornumber.shape[0], index.shape[0]/float(meta.sensornumber.shape[0])*100 )
                plt.subplot(1,2,1)
 		tempDF.groupby(tempDF.index.hour).mean()[index].mean(axis=1).plot(linewidth = 3,              
                	yerr= tempDF.groupby(tempDF.index.hour).mean().loc[index].std(axis=1),
                        label = lab)

                plt.text(15.5, tempDF.groupby(tempDF.index.hour).mean()[index].mean(axis=1).max(),
                '%2.1f'%tempDF.groupby(tempDF.index.hour).mean()[index].mean(axis=1).max(),
                 bbox=dict(facecolor='white', edgecolor = '#636363')
                 )
                plt.text(5.5, tempDF.groupby(tempDF.index.hour).mean()[index].mean(axis=1).min(),
                 '%2.1f'%tempDF.groupby(tempDF.index.hour).mean()[index].mean(axis=1).min(),
                 bbox=dict(facecolor='white', edgecolor = '#636363')
                 )
                
                plt.subplot(1,2,2)
                #print 'now plotting map for %s'%option
                x = meta['location:Longitude'][index].values
                y = meta['location:Latitude'][index].values
                marker_size = 150
                for j in index:
                        ax.scatter(meta['location:Longitude'][j], meta['location:Latitude'][j], 
                                  s = marker_size, 
                          c = pd.tools.plotting._get_standard_colors(3)[i], 
                          label = option if j == index[0] else "", 
                          )
                #plt.hold(True)
                i = i+1

            else : 
                print 'skipping plot %s'%title
            #annotate maximum and minimum
        plt.subplot(1,2,1)
        plt.title(title)
        plt.ylabel('Temperature in $^\circ $C')
        plt.xlabel('Hour')
        plt.legend(bbox_to_anchor = (.5, -.28),
                    loc=8,
                   borderaxespad=0.)
        plt.ylim([10,35])
        
        plt.subplot(1,2,2)
        lgd = plt.legend(#[options[sorttype]], 
                   bbox_to_anchor = (.5, -.28),
                    loc=8,
                   borderaxespad=0.)        
        #matplotlib.rcParams.update({'font.size': 22})
        filename = './plots/diurnal%s%s%s.eps'%(sorttype,sorttype2, option2)
        plt.savefig(filename, format = 'pdf', dpi = 600,) 


# function draws plots of diurnal data given a pandas dataframe that has an 'hour' column added to it 
# with optional 2nd category to sort by 

def diurnalplotsgeneral(tempDF,meta, parks, filename):

        fig  = plt.figure(figsize=(30, 12))
        plt.subplot(1,2,2)
        
        m = Basemap(llcrnrlon=meta['location:Longitude'].min()-.005,
            llcrnrlat=meta['location:Latitude'].min()-.0005,
            urcrnrlon=meta['location:Longitude'].max()+.005,
            urcrnrlat=meta['location:Latitude'].max()+.0005,
            projection='mill',
            #projection = 'merc',
            resolution ='h',
            #area_thresh=1000.
            epsg=3857
            )
        wms_server = "http://osm.woc.noaa.gov/mapcache" 
        #m.wmsimage(wms_server, layers = ["osm"], verbose = False)
        # loop over every curve to draw 
        i = 0
        for key in parks.keys():
            # if there is only one category to find 
            index = parks[key]                
            # make sure that there was data found satisfying the criteria 
            if index.shape[0] >0 :
                lab = '%s, %s/%s sensors (%2.1f %%)'%(key, index.shape[0], meta.sensornumber.shape[0], index.shape[0]/float(meta.sensornumber.shape[0])*100 )
                plt.subplot(1,2,1)
                tempDF.groupby(tempDF.index.hour).mean()[index].mean(axis=1).plot(linewidth = 3,
                                                                                                 yerr= tempDF.groupby(tempDF.index.hour).mean()[index].std(axis=1),
                                                                                                 label = lab)

                plt.text(15.5, tempDF.groupby(tempDF.index.hour).mean()[index].mean(axis=1).max(),
                 '%2.1f'%tempDF.groupby(tempDF.index.hour).mean()[index].mean(axis=1).max(),
                 bbox=dict(facecolor='white', edgecolor = '#636363')
                 )
                plt.text(5.5, tempDF.groupby(tempDF.index.hour).mean()[index].mean(axis=1).min(),
                 '%2.1f'%tempDF.groupby(tempDF.index.hour).mean()[index].mean(axis=1).min(),
                 bbox=dict(facecolor='white', edgecolor = '#636363')
                 )
                
                plt.subplot(1,2,2)
                #print 'now plotting map for %s'%option
                x = meta['location:Longitude'][index].values
                y = meta['location:Latitude'][index].values
                marker_size = 150
                for j in index:
                        m.scatter(meta['location:Longitude'][j], meta['location:Latitude'][j], 
                                  s = marker_size, 
                          c = pd.tools.plotting._get_standard_colors(3)[i], 
                          label = key if j == index[0] else "", 
                          latlon = True)
                #plt.hold(True)
                i = i+1

            else : 
                print 'skipping plot %s'%title
            #annotate maximum and minimum
        plt.subplot(1,2,1)
        title = 'Diurnal Cycle'
        plt.title(title)
        plt.ylabel('Temperature in $^\circ $C')
        plt.xlabel('Hour')
        plt.legend(bbox_to_anchor = (.5, -.28),
                    loc=8,
                   borderaxespad=0.)
        plt.ylim([18,32])
        
        plt.subplot(1,2,2)
        lgd = plt.legend(#[options[sorttype]], 
                   bbox_to_anchor = (.5, -.28),
                    loc=8,
                   borderaxespad=0.)        
        matplotlib.rcParams.update({'font.size': 22})
#        filename = './plots/diurnalparks%s.eps'%(parks.keys())
        plt.savefig(filename, format = 'eps', dpi = 600, 
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.show()
        
