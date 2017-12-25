fmajorAstroUtils
======================

Developing and debuging now, please wait until this Line disappear
--------------------------------------------------

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
* pyqt: you should install it from github by

    ``pip install [--user] git+https://github.com/ericmandel/pyds9.git#egg=pyds9``
* ds9: make sure you can open ds9 from your terminal.
* pyqt4 (optional): you will have a contral panel when open 1D array if you install ``pyqt4``.

In macOS, you can install it by
```
    brew install sip
    brew tap cartr/qt4
    brew tap-pin cartr/qt4
    brew install qt@4
    brew install cartr/qt4/pyqt@4 --with-python3
```

Install
-------
```
python3 setup.py install
```
