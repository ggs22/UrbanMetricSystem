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
realized by Jeremy Gelb according to [this article](./Tellier_et_al-2018-Regional_Science_Policy__26_Practice.pdf).

##### Install dependencies:
This project uses Python 3.8
##### Ubuntu
###### 1 Install GDAL binaries and headers

    sudo apt install libgdal-dev gdal-bin

###### 2 Set environment variables

    export CPLUS_INCLUDE_PATH=/usr/include/gdal
    export C_INCLUDE_PATH=/usr/include/gdal

###### 3 (Recommended) create a virtual Python environment

    python3 -m venv ~/.venv_ums
    source ~/.venv_ums/bin/activate

###### 4 Install osgeo (python binding for both GDAL and OGR)
    pip3 install --global-option=build_ext --global-option="-I/usr/include/gdal" GDAL==`gdal-config --version`

###### 5 install Python packages batch

    # from root directory of this project
    pip3 install -r requirements.txt
