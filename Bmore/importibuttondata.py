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
#	meta = pd.DataFrame(pd.read_csv(metafile, sep = ','))
	meta = meta.set_index(meta['sensornumber']) # set the row names for metadata to be the sensor number
	date_spec = {'Date': [ 0]}

	for index in meta.index: # loop over the sensor numbers in metadata and check if they are numbers
	    try : 
		str.isdigit(index)
		if str.isdigit(index)== 0: 
		    meta = meta.drop(meta.loc[index]) # delete the rows without sensor numbers
	    except TypeError: 
		pass
	meta = meta.sort(axis=0) # sort the data in order of the sensor number 
	meta = meta.set_index(meta['sensornumber'].astype(int))

	# read in all the csv files into a Panda's dataframe, skipping the headers of the csv files
	frames = []
	sensornumbers = []
	for file in files:
	    try : 
		colnumber = int(os.path.splitext(os.path.basename(file))[0][0:-1])
		frame = pd.read_csv(file, sep = ',', 
            		skiprows = 19, 
            		parse_dates = date_spec, 
            		keep_date_col=True)
		frame = frame.set_index('Date').rename(columns = {'Value': colnumber})
	
		try: 
			startdate = pd.to_datetime(meta['time'][colnumber]) #when the sensor was installed according to metadata
			ind = np.argmin([abs(frame.index -startdate)]) # Find the closest recording time to when put out
			frames.append(frame[ind+1:]) # only save the data from the hour after the sensor was installed	
		except KeyError: 
			frames.append(frame)
		#ind = np.argmin([abs(frame.index -startdate)]) # Find the closest recording time to when put out
		#frames.append(frame[ind+1:]) # only save the data from the hour after the sensor was installed	
		#frames.append(frame)
		#frames.append(pd.read_csv(file, sep = ',', skiprows = 19, parse_dates = date_spec, keep_date_col=True).set_index('Date').rename(columns = {'Value': int(os.path.splitext(os.path.basename(file))[0][0:-1])}))
	    	
	    except ValueError: 
		print "oops... something went wrong"
	   
	tempDF = pd.concat(frames, axis =1).resample('H')
	#tempDF = pd.concat(frames, axis =1).resample('H').dropna()
	meta = meta.loc[np.intersect1d(tempDF.columns.values, meta.sensornumber.values)]

	clim = tempDF.mean(axis=1) # the temperature 'climatology'
	anomaly = tempDF - np.tile(clim[:,np.newaxis], (1, tempDF.shape[1])) # the anomaly data
	anomalyDF = pd.DataFrame(anomaly, tempDF.index, tempDF.columns).sort(axis=1)

	meta = meta.loc[np.intersect1d(tempDF.columns.values, meta.sensornumber.values)]
	return tempDF, anomalyDF, meta



