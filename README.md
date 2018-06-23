# Pointless Spatial Modeling

## Summary

This repository hosts the code and data necessary to run the simulation and data analysis in the paper "Pointless continuous spatial surface reconstruction" by Wilson and Wakefield.

The algorithm to run the fully Bayesian Scotland example can take approximately a week to run; therefore, we include output that can be directly read in.

## Population data

Population data for Kenya and Scotland were obtained from the Center for International Earth Science Information Network (CIESIN). Information on the datasets can be found [here](http://sedac.ciesin.columbia.edu/binaries/web/sedac/collections/gpw-v4/gpw-v4-documentation.pdf) and can be downloaded [here](http://sedac.ciesin.columbia.edu/data/sets/browse). We used 2000 as the target year. To run the code, the data needs to be downloaded and can then be cropped and projected. For example:

```R
pop.dat <- raster("Data/Population/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2000.tif")
pop.dat <- setMinMax(pop.dat)
scot.extent <- extent(-8.7, -0.5, 54, 61)
pop.dat.scot <- crop(pop.dat, scot.extent)
pop.dat.scot <- projectRaster(pop.dat.scot, crs=CRS("+init=epsg:27700"))
```

where `gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2000.tif` is the unzipped file downloaded from CIESIN.

Gridded Population of the World, Version 4 (GPWv4): Population Count Adjusted to Match 2015 Revision of UN WPP Country Totals Center for International Earth Science Information Network - CIESIN - Columbia University. 2016. Gridded Population of the World, Version 4 (GPWv4): Population Count Adjusted to Match 2015 Revision of UN WPP Country Totals. Palisades, NY: NASA Socioeconomic Data and Applications Center (SEDAC). http://dx.doi.org/10.7927/H4SF2T42. Accessed 22 November 2016.

## Kenya survey coordinates

As an example, we use 399 survey coordinates (which are jittered) from the 2008 Kenya DHS and simulate an additional cluster location (not from the DHS). The locations, with a simulated number of households surveyed is provided in `Data/Kenya/kenya_dataNEW.RData`.

## Other information

- Approximate household size in Kenya was taken from [here](https://dhsprogram.com/pubs/pdf/fr308/fr308.pdf).

