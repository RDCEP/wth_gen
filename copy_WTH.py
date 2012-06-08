#!/usr/bin/env python
from string import split
from numpy import floor
from os import system, path
from mpl_toolkits.basemap import NetCDFFile
from pylab import find

in_dir = '/gpfs/pads/projects/see/CCSM_Library/nc_test_1/'
out_dir = '/gpfs/pads/projects/see/CCSM_Library/nc_test_1_copy/afri/'

f = NetCDFFile('/gpfs/pads/projects/see/data/raw/isi_test/precip/precip_2001.nc')
lat_lores_array = f.variables['lat'][:]
lon_lores_array = f.variables['lon'][:]
f.close()

fname = '/autonfs/home/dmcinern/new_wth_gen/continents/grid_total_mask_afri.csv'
f = open(fname)
for line in f:

  line_array = split(line,',')
  grid_id_hires = int(line_array[0]) 
  lat_hires = float(line_array[1])
  lon_hires =  float(line_array[2]) 
  grid_id_lores = float(line_array[3])
  lat_lores = float(line_array[4])
  lon_lores =  float(line_array[5])

  lat_lores_index = find(lat_lores_array == lat_lores)[0] + 1
  lon_lores_index = find(lon_lores_array + 180. == lon_lores )[0] + 1

  lat_lon_id_lores = '%03d_%03d' % (lat_lores_index,lon_lores_index)

  old_sub_dir = in_dir + '/' + lat_lon_id_lores
  new_sub_dir = out_dir + '/' + str(grid_id_hires) 

  if not path.exists(old_sub_dir):
    print 'path ' + old_sub_dir + ' does not exist'

  else:

    system('mkdir ' + new_sub_dir)

    fin = open(old_sub_dir + '/GENERIC1.WTH')
    fout = open(new_sub_dir + '/GENERIC1.WTH','w')
    line_num = 0
    for line in fin:
      if line_num == 0: 
        str_old = 'cell ' + lat_lon_id_lores
        str_new = 'cell %i' % grid_id_hires
        fout.write('hi')
        fout.write( line.replace( str_old, str_new ) )
      elif line_num == 2:    
        str_old = 'CI%9.4f%9.4f' % (lat_lores,lon_lores-180.)
        str_new = 'CI%9.4f%9.4f' % (lat_hires,lon_hires-180.)
        fout.write( line.replace( str_old, str_new ) )
      else:
        fout.write( line )
      line_num += 1
    fin.close()
    fout.close()






