# Urban Metric System
> "An urban metric system (UMS) is here implemented in order to henceforth solve the
> numerous problems stemming from of the various concepts of urban areas. It is based on space‐
> the absence of a mathematical definition of the boundaries economy, especially, land‐rent theory,
> the concepts of attractive and repulsive forces, as well as vector fields. [...] City. Based on
> the proposed UMS, a synthetic index of urban sprawl is defined and computed for the studied
> metropolitan
> areas."

> -Tellier et al (2018)

###### v0.1
Version 0.1 is the version of the Urban Metric System program as
realized by Jeremy Gelb according to [this article](./Tellier_et_al-2018-Regional_Science_Policy__26_Practice.pdf)


###### installation:

####### Unix

GDAL (osgeo)
https://trac.osgeo.org/gdal/wiki/GdalOgrInPython

https://stackoverflow.com/questions/43587960/gdal-installation-error-using-pip

https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal

https://ljvmiranda921.github.io/notebook/2019/04/13/install-gdal/#using-your-package-manager


pip install --global-option=build_ext --global-option="-IC:\\repos\\UrbanMetricSystem\\2.2\\gdal\\port\\" GDAL==`gdal-config --version`

####### windows
In Anaconda:
    -install libgdal
    -install gdal
In command prompt:
    -install pip packages using "pip install -r requirements.txt"


