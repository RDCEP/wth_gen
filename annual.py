#!/usr/bin/env python
import sys
from os import system, listdir
from fnmatch import fnmatch
import calendar
from mpl_toolkits.basemap import NetCDFFile

outdir = sys.argv[1]

for fileName in listdir('./'):
  if fnmatch(fileName,'*.nc'):
    print 
    print fileName
    fileName_split = fileName.split('_')
    varName = fileName_split[0]
    years = fileName_split[-1][:-3]

    if '-' not in years:
      command = 'cp ' + fileName + ' ' + outdir
      print command
#      system(command)

    else:  
      year_0 = int(years[:4])
      year_1 = int(years[-4:])
      print year_0, year_1

      time = NetCDFFile(fileName).variables['time'][:]
      day_0 = time[0]

      for year in range(year_0,year_1+1):
        print year
        if calendar.isleap(year): 
          n_days = 366
        else:
          n_days = 365
        print n_days
        day_1 = day_0 + n_days - 1
        print day_0, day_1
        new_fileName = fileName.replace(years,str(year))
        command = 'ncks -d lat,-60.,67. -d time,' + str(day_0) + ',' + str(day_1) + ' ' + fileName + ' ' + outdir + '/' + new_fileName
        print command
#        system(command)

        day_0 = day_1 + 1
           
       


