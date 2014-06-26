## Experiments in detection and prediction of El Nino/La Nina events.

#### R/netcdf-convertor.R.

This converts netCDF files containing surface temperature data from NOAA (eg  air.sig995.1951.nc) into a format more easily readable by other programs in other languages. You will need to edit it to your requirements:

* Change the working directory 
* Change the ranges of latitude and longitude
* Change the range of years

As supplied, it converts 3 years for 4 grid points covering Scotland. I’ve put the Ludescher *et al* Pacific co-ordinates in comments. More instructions are in the script. 

Then start R, and copy and paste the whole file into the R console. (There are other ways of running R scripts but this is simplest for novices).

#### R/grj/ludescher.R

Aimed at replicating Ludescher et al, 2013. As of 26 June 2014, it is close, but not identical.

```
