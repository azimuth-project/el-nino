## Experiments in El Nino analysis and prediction

This software is connected to the Azimuth Code Project **[Experiments in El Nino analysis and prediction](http://www.azimuthproject.org/azimuth/show/Experiments%20in%20El%20Ni%C3%B1o%20analysis%20and%20prediction)**. 

#### R / netcdf-convertor.R.

This converts netCDF files containing surface temperature data from NOAA (eg  air.sig995.1951.nc) into a format more easily readable by other programs in other languages. You will need to edit it to your requirements:

* Change the working directory 
* Change the ranges of latitude and longitude
* Change the range of years

As supplied, it converts 3 years for 4 grid points covering Scotland. I’ve put the Ludescher *et al* Pacific co-ordinates in comments. More instructions are in the script. 

Then start R, and copy and paste the whole file into the R console. (There are other ways of running R scripts but this is simplest for novices).

#### R / grj / ludescher.R

Aimed at replicating Ludescher et al, 2013. As of 26 June 2014, it is close, but not identical.  For an explanation see  [Part 4](http://johncarlosbaez.wordpress.com/2014/07/08/el-nino-project-part-4/) of the El Ni&ntilde;o Project series.


#### R / average-link-strength.txt

This file has the average link strength S as computed by `ludescher.R` at 10-day intervals, starting from day 730 and going until day 12040, where day 1 is the first of January 1948.  For an explanation see  [Part 4](http://johncarlosbaez.wordpress.com/2014/07/08/el-nino-project-part-4/) of the El Ni&ntilde;o Project series.

#### R / average-link-strength-1948-2013.txt

The second column in this file lists the average link strengths S as computed by `ludescher.R` at 10-day intervals, starting from day 730, and going until day 24090, where day 1 is the first of January 1948.  The first column numbers these items from 1 to 2337.  For an explanation see  [Part 4](http://johncarlosbaez.wordpress.com/2014/07/08/el-nino-project-part-4/) of the El Ni&ntilde;o Project series.

#### R / average-link-strength-daily.txt

The second column in this file lists the average link strengths S as computed by `ludescher.R` at daily intervals, starting from day 730, and going until day 24090, where day 1 is the first of January 1948. The first column numbers these items from 730 to 24090. For an explanation see [Part 4](http://johncarlosbaez.wordpress.com/2014/07/08/el-nino-project-part-4/) of the El Niño Project series.

#### R / grj / covariances-basin-vs-rest.R

Makes maps of the Pacific, one per quarter from 1951 to 1979, showing covariances of grid points with the "Ludescher et al basin"

```
