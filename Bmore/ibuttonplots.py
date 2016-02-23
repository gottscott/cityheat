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
from mpl_toolkits.basemap import Basemap

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

    plt.axvline(data.std()+data.mean(), 
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
        # define dictionaries of options and titles
        options = {'landcoverclass': ['impervious', 'grass', 'soil'],#'dirt'],
                   'sunorshade': ['sun', 'partial', 'shade'],
                   'attachment': ['metal', 'wood', 'tree']}

        titles = {'landcoverclass': 'Land Cover Class',
                   'sunorshade': 'Shadiness',
                   'attachment': 'Attachment'}
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
                        m.scatter(meta['location:Longitude'][j], meta['location:Latitude'][j], 
                                  s = marker_size, 
                          c = pd.tools.plotting._get_standard_colors(3)[i], 
                          label = option if j == index[0] else "", 
                          latlon = True)
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

def mapmean(tempDF, meta, name = ''): 
    #fig  = plt.figure(figsize=(30, 30))
    x = meta['location:Longitude'].values
    y = meta['location:Latitude'].values
    c = tempDF.mean(axis=0).values # The colors will show the mean temp at that location
    marker_size = 150
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])

    m = Basemap(llcrnrlon=meta['location:Longitude'].min()-.005,
                llcrnrlat=meta['location:Latitude'].min()-.005,
                urcrnrlon=meta['location:Longitude'].max()+.005,
                urcrnrlat=meta['location:Latitude'].max()+.005,
                projection='mill',
                #projection = 'merc',
                resolution ='h',
                #area_thresh=1000.
                epsg=3857
                )
    wms_server = "http://osm.woc.noaa.gov/mapcache" 
#    m.wmsimage(wms_server, layers = ["osm"], verbose = True)
    #define the color map 
    cmap = matplotlib.cm.RdBu_r
    #bounds = np.linspace(round((np.nanmean(c)-3)*9/5.+32),round((np.nanmean(c)+3)*9/5.+32),13)
#    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
#    
#    fahr = m.scatter(x,y, s = marker_size, c =(c*9./5. + 32.),
#                     cmap = cmap, 
#                     norm = norm, 
#                     #cmap = matplotlib.cm.RdBu_r, 
#                        latlon = True, 
#                        lw = 2, 
#                        edgecolor = 'gray', 
#                        #vmin = (c.mean()-3)*9/5., 
#                        #vmax = (c.mean()+3)*9/5.,
#                        )
    cmap = matplotlib.cm.YlOrRd #RdBu_r
    bounds = np.linspace(round(np.nanmean(c))-3,round(np.nanmean(c))+3,7)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    
    celsius = m.scatter(x,y, s = marker_size, 
                        c = c, 
                        cmap = cmap, 
                        norm = norm,
                        #cmap = matplotlib.cm.RdBu_r, 
                        latlon = True, 
                        lw = 2, 
                        edgecolor = 'gray', 
                        vmin = np.nanmean(c)-3, 
                        vmax = (np.nanmean(c))+3,
                        )
    cbar1 = m.colorbar(celsius, location = 'right', label = 'Temperature in $^\circ $C')
    #cbar2 = m.colorbar(fahr, location = 'bottom', label = 'Temperature in $^\circ $F')
    #ax1 = fig.add_axes([.8, .15, 0.05, .7] )
    #cbar2 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm = norm, orientation = 'vertical')
    
    lon = x[np.nanargmin(c)]
    lat = y[np.nanargmin(c)]
    a,b = m(lon,lat)
    plt.text(a,b, '%2.1f'%np.nanmin(c),
             ha = 'center', 
             va = 'bottom',
             fontsize=16, 
             #fontweight='bold',
             color = 'k')
    
    lon = x[np.nanargmax(c)]
    lat = y[np.nanargmax(c)]
    a,b = m(lon,lat)
    plt.text(a,b, '%2.1f'%np.nanmax(c), 
             ha = 'center', 
             va = 'bottom',
             fontsize=18, 
             #fontweight='bold',
             color = 'k')
    
    plt.title('Mean Temperature %s'%name)

    filename = './plots/meantempmap%s.eps'%name
    plt.savefig(filename, format = 'eps', dpi = 600)

# function draws plots of diurnal data given a pandas dataframe that has an 'hour' column added to it 
# with optional 2nd category to sort by 
def diurnalplots(tempDF, meta, sorttype, sorttype2=0, option2=0):
        # define dictionaries of options and titles
        options = {'landcoverclass': ['impervious', 'grass', 'soil'],#'dirt'],
                   'sunorshade': ['sun', 'partial', 'shade'],
                   'attachment': ['metal', 'wood', 'tree']}

        titles = {'landcoverclass': 'Land Cover Class',
                   'sunorshade': 'Shadiness',
                   'attachment': 'Attachment'}

        fig  = plt.figure(figsize=(30, 12))
        plt.subplot(1,2,2)
        matplotlib.rcParams.update({'font.size': 22})
        
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
                	yerr= tempDF.groupby(tempDF.index.hour).mean()[index].var(axis=1),
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
                          label = option if j == index[0] else "", 
                          latlon = True)
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
        plt.savefig(filename, format = 'eps', dpi = 600, 
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.tight_layout
        plt.show()

# function draws plots of diurnal data given a pandas dataframe that has an 'hour' column added to it 
# with optional 2nd category to sort by 

def diurnalplotsgeneral(diurnalDF,meta, parks, filename):

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
                                                                                                 yerr= tempDF.groupby(tempDF.index.hour).mean()[index].var(axis=1),
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
        
