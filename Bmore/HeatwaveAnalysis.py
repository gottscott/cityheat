
# coding: utf-8

# In[2]:

import numpy as np
import matplotlib.pyplot as plt
import iris
import glob
import time
import netCDF4


# In[3]:




###### 

# In[6]:




# In[121]:

# Heat waves as defined by H01-H04, mean daily temperature > median threshold for consecutive tdays
files =sorted(glob.glob('/data3/blaustein/waugh/NCEP/daily/air*.nc'))
cube = iris.load_cube(files[0], iris.Constraint(Level = 1000))
climatologies = {} 
thresholds = np.array([81, 90, 95, 97.5, 98, 99])
percentilethresholds = np.zeros([cube.shape[2], cube.shape[1],thresholds.size])
toc = time.time()
for x in range(0,cube.shape[2]) : 
    lon = cube.coord('longitude').points[x]
    #for y in range(0,37) : #cube.shape[1]) : # for northern hemisphere
    for y in range(36,73) : #cube.shape[1]) : # for southern hemisphere
        lat = cube.coord('latitude').points[y] 
        cubelist = iris.load(files, iris.Constraint(Level = 1000) & iris.Constraint(longitude =lon) & iris.Constraint(latitude = lat))
        temp = np.array([])
        timearray = np.array([])
                
        for c in cubelist : 
            temp = np.append(temp, c.data)
            timearray = np.append(timearray, c.coord('time').points)

        percentilethresholds[x,y,:] = np.array(np.percentile(temp, list(thresholds))) # compute the 90th, 95, 95th and 99th percentile for daily mean temperature
        for i in [1,2,4,5]:
            hotIndex = np.array(np.where(temp >= percentilethresholds[x,y,i]))[0]
            grad = hotIndex[1:]-hotIndex[0:-1]
            heatEventIndex = np.array(np.where(grad>6)) 
            eventStart = timearray[heatEventIndex.astype(int)] 
            eventEnd = timearray[heatEventIndex.astype(int)+1] 
            eventDuration = eventEnd - eventStart 
            ind = np.where(eventDuration > 48)[0]
            if ind.size > 0: 
                    climatologies[thresholds[i], x, y] = eventStart[ind], eventEnd[ind], eventDuration[ind] 
    np.save('HeatWaveClimatologyMedianThreshold',climatologies)                
tic = time.time()
print tic-toc


# In[150]:

np.save('percentilethresholds',percentilethresholds)


# In[86]:




# In[122]:

# Heat waves as defined by H05-H07, min or max daily temperature > median threshold for consecutive days
files =sorted(glob.glob('/data3/blaustein/waugh/NCEP/subdaily/air*.nc'))
cube = iris.load_cube(files[0], iris.Constraint(Level = 1000))
mindailyclimatologies = {} 
maxdailyclimatologies = {} 

toc = time.time()

for x in range(0,cube.shape[2]) : 
    lon = cube.coord('longitude').points[x]
    #for y in range(0,37) : #cube.shape[1]) : # for northern hemisphere
    for y in range(36,73) : #cube.shape[1]) : # for southern hemisphere
        lat = cube.coord('latitude').points[y] 
        cubelist = iris.load(files, iris.Constraint(Level = 1000) & iris.Constraint(longitude =lon) & iris.Constraint(latitude = lat))
        temp = np.array([])
        timearray = np.array([])
                
        for c in cubelist : 
            temp = np.append(temp, c.data)
            timearray = np.append(timearray, c.coord('time').points)
            
        mintemp = temp.reshape(temp.shape[0]/4, 4).min(1)
        maxtemp = temp.reshape(temp.shape[0]/4, 4).max(1)
        
        #minimum temp analysis
        i95 = np.where(thresholds ==95)
        hotIndex = np.array(np.where(mintemp >= percentilethresholds[x,y,i95]))[0] #where min daily temp exceeds 95% threshold
        grad = hotIndex[1:]-hotIndex[0:-1]
        heatEventIndex = np.array(np.where(grad>6)) 
        eventStart = timearray[heatEventIndex.astype(int)] 
        eventEnd = timearray[heatEventIndex.astype(int)+1] 
        eventDuration = eventEnd - eventStart 
        ind = np.where(eventDuration > 48)[0]
        if ind.size > 0: 
            mindailyclimatologies[x, y] = eventStart[ind], eventEnd[ind], eventDuration[ind] 

        #maximum temp analysis
        hotIndex = np.array(np.where(maxtemp >= percentilethresholds[x,y,i95]))[0] #where max daily temp exceeds 95% threshold
        grad = hotIndex[1:]-hotIndex[0:-1]
        heatEventIndex = np.array(np.where(grad>6)) 
        eventStart = timearray[heatEventIndex.astype(int)] 
        eventEnd = timearray[heatEventIndex.astype(int)+1] 
        eventDuration = eventEnd - eventStart 
        ind = np.where(eventDuration > 48)[0]
        if ind.size > 0: 
            maxdailyclimatologies[x, y] = eventStart[ind], eventEnd[ind], eventDuration[ind] 
    np.save('HeatWaveClimatologyMinThreshold',mindailyclimatologies)
    np.save('HeatWaveClimatologyMaxThreshold',maxdailyclimatologies) 
tic = time.time()
print tic-toc


# In[11]:




# In[149]:

# maximum daily heat index 
files1 =sorted(glob.glob('/data3/blaustein/waugh/NCEP/subdaily/air*.nc'))
files2 =sorted(glob.glob('/data3/blaustein/waugh/NCEP/subdaily/rhum*.nc'))
cube = iris.load_cube(files1[0], iris.Constraint(Level = 1000))
hiclimatologies = {}
# coefficients for heat index calculation 
c1 = -42.379
c2 = 2.04901523
c3 = 10.14333127
c4 = - 0.22475541
c5 = - 0.00683783
c6 = - 0.05481717
c7 = 0.00122874
c8 = 0.00085282
c9 = - 0.00000199

#heat index thresholds, 
hithresholds = np.array([80, 90, 105, 130])

for x in range(0,cube.shape[2]) : 
    lon = cube.coord('longitude').points[x]
    #for y in range(0,37) : #cube.shape[1]) : # for northern hemisphere
    for y in range(36,73) : #cube.shape[1]) : # for southern hemisphere 
        lat = cube.coord('latitude').points[y] 
        tempcubelist = iris.load(files1, iris.Constraint(Level = 1000) & iris.Constraint(longitude =lon) & iris.Constraint(latitude = lat))
        rhcubelist = iris.load(files2, iris.Constraint(longitude =lon) & iris.Constraint(latitude = lat))

        temp = np.array([])
        rh = np.array([])
        timearray = np.array([])
                
        for i in range(0,len(tempcubelist)-1) : 
            rh = np.append(rh, rhcubelist[i].data)
            temp = np.append(temp, tempcubelist[i].data)
            timearray = np.append(timearray, rhcubelist[i].coord('time').points)
        
        # convert temp to fahrenheit, weap for the state of science in america
        temp = (temp-273)*9/5. +32
        hi = c1 +c2*temp +c3*rh + c4*temp*rh +c5*temp**2 +c6*rh**2+c7*temp**2*rh+c8*temp*rh**2+c9*temp**2*rh**2
        
        #maximum temp analysis
        for maxdailyhi in hithresholds : 
            hotIndex = np.array(np.where(hi >= maxdailyhi))[0] #where max daily temp exceeds 95% threshold
            grad = hotIndex[1:]-hotIndex[0:-1]
            heatEventIndex = np.array(np.where(grad>6)) 
            eventStart = timearray[heatEventIndex.astype(int)] 
            eventEnd = timearray[heatEventIndex.astype(int)+1] 
            eventDuration = eventEnd - eventStart 
            hiclimatologies[maxdailyhi, x, y] = eventStart[ind], eventEnd[ind], eventDuration[ind] 

    np.save('HeatWaveClimatologyHI',hiclimatologies) 
tic = time.time()
print tic-toc


# In[11]:




# In[11]:




# In[ ]:



