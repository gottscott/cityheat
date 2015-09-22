# import ibutton data located in ./data/ and put it in a pandas dataframe
# In[4]:

import numpy as np 
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd
import matplotlib 
import matplotlib.pylab as pylab
from mpl_toolkits.basemap import Basemap
pylab.rcParams['figure.figsize'] = 12, 8
pd.options.display.mpl_style = 'default'
get_ipython().magic(u'matplotlib inline')

def importdata(files,meta) : 
	# read in all files in current directory, assuming that they have the file struction nt.csv or nh.csv, eg 1t.csv, 13t.csv, 40h.csv, etc 
	#files = glob.glob('./data/*[tT].csv') # read in data files 

	#meta = pd.DataFrame(pd.read_csv('./data/bmoremetadata.csv', sep = ',')) # read in metadata 
	meta = meta.set_index(meta['sensornumber']) # set the row names for metadata to be the sensor number

	for index in meta.index: # loop over the sensor numbers in metadata and check if they are numbers
	    try: 
		if str.isdigit(index)== 0: 
		    meta = meta.drop(meta.loc[index]) # delete the rows without sensor numbers
	    except (TypeError):
		pass
		
	meta = meta.sort(axis=0) # sort the data in order of the sensor number 
	meta = meta.set_index(meta['sensornumber'].astype(int))
	# read in all the csv files into a Panda's dataframe, skipping the headers of the csv files
	frames = []
	sensornumbers = []
	date_spec = {'Date': [ 0]}
	for file in files: 
	    try : 
		frames.append(pd.read_csv(file, sep = ',', skiprows = 19, parse_dates = date_spec, keep_date_col=True))
		sensornumbers.append(int(os.path.splitext(os.path.basename(file))[0][0:-1]) )
	    except ValueError: 
		print "oops... something went wrong"
	    
	data = pd.concat(frames, axis =1)

	# Clean up data and build a new dataset
	# next, normalize the dates: the sensors probably weren't set to all go off at the same time
	# assume they all are within an hour of each other
	maxStart = data['Date'].values[0,:].max() # the start time of the last sensor to start
	minEndIndex = pd.isnull(data['Date']).any(1).nonzero()[0].min() # find the earliest row instance of Nans
	maxEnd = data['Date'].values[minEndIndex-1,:].max() #picks the latest time in that row
	rng = pd.date_range(maxStart, maxEnd, freq='H') # make a date range from maxStart to maxEnd, with hourly frequency

	# produces a data structure with time from variable rng as the rows and sensornumbers as the column heading
	# the column 
	tempDF = pd.DataFrame(data['Value'].values[0:minEndIndex,:], rng[:], sensornumbers).sort(axis=1)

	# Clean all flipped data
	#tempDF1 = tempDF # make a copy 
	#flipped1 = np.where(tempDF['June 12 2015 05'].values < 23)[1]
	#nonflipped = np.setdiff1d(meta.sensornumber.values, meta.sensornumber.values[flipped1])
	#tempDF = tempDF[nonflipped] #eliminate all the nonflipped
	#tempDF = tempDF.drop(40,1)

	#anomalyDF = anomalyDF[nonflipped]
	#anomalyDF = anomalyDF.drop(40,1)
	#clean the metadata so that missing sensors don't mess with analysis
	meta = meta.loc[np.intersect1d(tempDF.columns.values, meta.sensornumber.values)]

	#### Now compute the anomaly statistics
	clim = tempDF.mean(axis=1) # the temperature 'climatology'
	anomaly = tempDF - np.tile(clim[:,np.newaxis], (1, tempDF.shape[1])) # the anomaly data
	anomalyDF = pd.DataFrame(anomaly, rng[:], tempDF.columns).sort(axis=1)

	meta = meta.loc[np.intersect1d(tempDF.columns.values, meta.sensornumber.values)]


	return tempDF, anomalyDF, meta



