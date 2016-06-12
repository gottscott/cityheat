
# coding: utf-8

# In[12]:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib
#%matplotlib inline
import re
from skimage import io, exposure
from matplotlib import pyplot as plt


# In[3]:

# open images
#band10 = mpimg.imread('LC80150332015229LGN00_B10.TIF')
#band11 = mpimg.imread('LC80150332015229LGN00_B11.TIF')


# In[32]:

def image_adjust_brightness(img, limit_left, limit_right, title):
    """
    Adjust image brightness and plot the image
    Input:
    img - 2D array of uint16 type
    limit_left - integer
    limit_right - integer
    title - string
    """
    img_ha = exposure.rescale_intensity(img, (limit_left, limit_right))
    
    fig = plt.figure(figsize=(10, 10))
    fig.set_facecolor('white')
    plt.imshow(img_ha, cmap='gray')
    plt.title(title)
    plt.show()
    
def landsat8LST(band10file, title, vmin = 20, vmax = 37): 
    band10 = mpimg.imread(band10file)
    # convert from digital numbers to (spectral) radiance
    m_l = 3.3420E-04 # gain
    a_l = 0.10000 # bias
    radiance = m_l * band10 + a_l
    # convert from radiance to temperature by inverting the planck function
    k1 = 774.8853 
    k2 = 1321.0789
    temp = k2/ (np.log(k1/radiance +1)) -273.15
    mpimg.imsave('LST%s.tif'%title, temp, format = 'TIF')
    temp[temp < -100] = 'nan'
    matplotlib.rcParams.update({'font.size': 60})
    fig  = plt.figure(figsize=(75, 75))
    imgplot = plt.imshow(temp, 
                         cmap = 'gray', #matplotlib.cm.RdBu_r, #matplotlib.cm.OrRd,
                         vmin = vmin, #20,
                         vmax = vmax, #37,
                         )
    plt.title(title)
    plt.colorbar()
    plt.savefig('LSTcolorbar%s.eps'%title)
    return temp, imgplot

if __name__ == '__main__' : 
    import fnmatch 
    import os 
    files = []
    for root, dirnames, filenames in os.walk('.'): 
        for filename in fnmatch.filter(filenames, '*B10.TIF'): 
            files.append(os.path.join(root, filename))
            
    for file in files: 
        title = os.path.splitext(os.path.basename(file))[0]
        landsat8LST(file, title)

