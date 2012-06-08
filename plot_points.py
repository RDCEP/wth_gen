#!/usr/bin/env python
from pylab import figure, plot, show
from read_2d_array import read_2d_array
from string import split
from mpl_toolkits.basemap import NetCDFFile

f = NetCDFFile('/gpfs/pads/projects/see/data/raw/isi_test/precip/precip_2001.nc')
lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
f.close()

print lat
print lon

lat_list = []
lon_list = []
f = open('grid_id_list.txt')
for line in f:
  line_array = split(line)
  lat_list.append(lat[int(line_array[1])-1])
  lon_list.append(lon[int(line_array[2])-1])

figure()
plot(lon_list,lat_list,'bo')
show()


