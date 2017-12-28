fmajorAstroUtils
======================

Fmajor's tools for astronomy data reduction.

* [imheader](doc/imheader.rst): A command line tool to check the headers or do filtering for fits file using header information.
    * You can quickly overview the structure of a fits (like [fv](https://heasarc.gsfc.nasa.gov/ftools/fv/), but in terminal)
    * You can filter fits by their header (e.g. only select fits with exptime>100)
    
* [ds10](doc/ds10.rst): A python wrapper of [ds9](http://ds9.si.edu/site/Home.html) (use [pyds9](https://github.com/ericmandel/pyds9)).    

    * It can open 2D image, 1D array or binary table by a single command "ds10 \*.fits". It will open the image using ds9, open a 1D array using matplotlib (also with a PyQt control window) and open binary table using ipython interactive mode
    
    * You can use the same arguments as imheader to filter the files and the extensions
    
    * You can modify the plot data, plot style and figure xlimits and ylimits using the PyQt control panel for 1D plot
    
    * You can easily copy regions between different frames
    
    * You can directly plot the projections in matplotlib (e.g. you want to check the projection of two 2D spectrums and see the difference)
    
    * You can tell ds9 want to do after it loads all the files (e.g. change the pan, zoom, rotate, limits...)

You will like these tools if you have numerious fits image to check.

Dependences
-----------
* python3
* pyds9: if the ds10 tell you that it can to load fits due to XPA error, maybe you should install the latest version of pyds9 from github, like

    ``pip3 install git+https://github.com/ericmandel/pyds9.git#egg=pyds9 --upgrade``
* ds9: make sure you can open ds9 from your terminal.
* pyqt (optional): you will have a contral panel when open 1D array if you install ``pyqt``. You should install ``qt`` first.


Install
-------
```
git clone https://github.com/Fmajor/fmajorAstroTools.git
cd fmajorAstroTools
make
```
