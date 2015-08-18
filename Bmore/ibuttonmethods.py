# Methods for analyzing and plotting ibutton data 

# make diurnal plots of data defined in the dataframe diurnalDF
# that has an hour as its index
# for each type of data (landcoverclass, sunorshade) 
# as defined in the metadata 
import matplotlib.pyplot as plt

def diurnalplots(diurnalDF, sorttype):
	import matplotlib.pyplot as plt	
	options = {'landcoverclass': ['impervious', 'dirt', 'grass', ], 
		   'sunorshade': ['sun', 'partial', 'shade'], 
		   'attachment': ['metal', 'deadwood', 'tree']}

	titles = {'landcoverclass': 'Land Cover Class', 
		   'sunorshade': 'Shadiness', 
		   'attachment': 'Attachment'}

	fig  = plt.figure(figsize=(15, 12))


	for option in options[sorttype]: 
	    index = meta['sensornumber'].values[np.where(meta[sorttype]==option)].astype(int)
	    lab = '%s, %s/%s sensors (%2.1f %%)'%(option, index.shape[0], meta.sensornumber.shape[0], index.shape[0]/float(meta.sensornumber.shape[0])*100 )
	    diurnalDF.set_index("hour").groupby(level=0).mean()[index].mean(axis=1).plot(linewidth = 3, 
												 yerr= diurnalDF.set_index("hour").groupby(level=0).mean()[index].var(axis=1), 
												 label = lab)
	    #annotate maximum and minimum
	    plt.text(15.5, diurnalDF.set_index("hour").groupby(level=0).mean()[index].mean(axis=1).max(), 
		 '%2.1f'%diurnalDF.set_index("hour").groupby(level=0).mean()[index].mean(axis=1).max(), 
		 bbox=dict(facecolor='white', edgecolor = '#636363')
		 )
	    plt.text(5.5, diurnalDF.set_index("hour").groupby(level=0).mean()[index].mean(axis=1).min(), 
		 '%2.1f'%diurnalDF.set_index("hour").groupby(level=0).mean()[index].mean(axis=1).min(), 
		 bbox=dict(facecolor='white', edgecolor = '#636363')
		 )

	plt.title('Diurnal Cycle by %s'%titles[sorttype])
	plt.ylabel('Temperature in $^\circ $C')
	plt.xlabel('Hour')
	plt.legend(bbox_to_anchor = (.5, -.28),
		    loc=8,
		   borderaxespad=0.)
	plt.ylim([18,32])
	matplotlib.rcParams.update({'font.size': 22})


