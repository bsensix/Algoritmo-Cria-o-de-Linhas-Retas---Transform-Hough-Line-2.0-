#!/usr/bin/env python
# coding: utf-8

# ##  Bibliotecas:

# In[1]:


#python and matplotlib libraries
import numpy as np
import math 
import matplotlib.pyplot as plt
from matplotlib import cm, gridspec

#geospatial libraries
from osgeo import gdal, ogr, osr
import rasterio
import geopandas as gpd 
import fiona

#machine learning libraries
from skimage import feature
from skimage import morphology
from skimage.transform import hough_line, hough_line_peaks
from skimage.morphology import medial_axis, skeletonize ,thin


# ##  Importar Raster

# In[2]:


caminho = r'C:\Users\breno\Desktop\TESTE_SRICPT\SCRIPT_LINHAS\1230426\t1.gpkg'
imgRaster = rasterio.open(caminho)
im = imgRaster.read(1)


# ##  Algoritmo de Esqueletização

# In[3]:


skeleton = morphology.skeletonize(im, method='zhang')


# ##  Transform Hough Line 

# In[4]:


selImage = skeleton

# Classic straight-line Hough transform


n_linhas = 6000
limite = 0.4

tested_angles = np.linspace(-np.pi / 2, np.pi / 2, int(n_linhas), endpoint=False)
h, theta, d = hough_line(selImage, theta=tested_angles)


# ##  Transformar para SHP

# In[5]:


#define representative length of interpreted lines in pixels
#selDiag = int(np.sqrt(selImage.shape[0]**2+selImage.shape[1]**2))
selDiag = selImage.shape[1]
progRange = range(selDiag)
totalLines = []
angleList = []
for _, angle, dist in zip(*hough_line_peaks(h, theta, d,threshold= limite * h.max())):
    (x0, y0) = dist * np.array([np.cos(angle), np.sin(angle)])
    if angle in [np.pi/2, -np.pi/2]:
        cols = [prog for prog in progRange]
        rows = [y0 for prog in progRange]
    elif angle == 0:
        cols = [x0 for prog in progRange]
        rows = [prog for prog in progRange]
    else:
        c0 = y0 + x0 / np.tan(angle)
        cols = [prog for prog in progRange]
        rows = [col*np.tan(angle + np.pi/2) + c0 for col in cols]

    partialLine = []
    for col, row in zip(cols, rows):
        partialLine.append(imgRaster.xy(row,col))
    totalLines.append(partialLine)
    if math.degrees(angle+np.pi/2) > 90:
        angleList.append(180 - math.degrees(angle+np.pi/2))
    else:
        angleList.append(math.degrees(angle+np.pi/2))

#show on sample of the points and lines
print(totalLines[0][:5])


# ##  Salvando em SHP

# In[6]:


#define schema
schema = {
    'properties':{'angle':'float:16'},
    'geometry':'LineString'
}

#out shapefile
outShp = fiona.open('croprows.shp',mode='w',driver='ESRI Shapefile',schema=schema,crs=imgRaster.crs)
for index, line in enumerate(totalLines):
    feature = {
        'geometry':{'type':'LineString','coordinates':line},
        'properties':{'angle':angleList[index]}
         }
    outShp.write(feature)
outShp.close()


# ##  Vetorização do Raster de Cana 

# In[7]:


dataset = gdal.Open(caminho)
band = dataset.GetRasterBand(1)

proj = dataset.GetProjection()
shp_proj = osr.SpatialReference()
shp_proj.ImportFromWkt(proj)

output_file = 'vetor.shp'
call_drive = ogr.GetDriverByName('ESRI Shapefile')
create_shp = call_drive.CreateDataSource(output_file)
shp_layer = create_shp.CreateLayer('layername', srs = shp_proj)
new_field = ogr.FieldDefn(str('ID'),ogr.OFTInteger)
shp_layer.CreateField(new_field)

gdal.Polygonize(band,band,shp_layer,0,[],callback=None)
create_shp.Destroy()
raster = None


# In[ ]:





# In[ ]:




