Source Extractor Wrapper for Python
===================================

The tiny `sewpy` module let's you run SExtractor as if it would all be native python... 

```python 
import sewpy
sew = sewpy.SEW(params=["X_IMAGE", "Y_IMAGE", "FLUX_RADIUS(3)", "FLAGS"])
out = sew("myimage.fits", config={"DETECT_MINAREA":10, "PHOT_FLUXFRAC":"0.3, 0.5, 0.8"})
print out["table"] # this is an astropy table.
```

... but also allows for a more sophisticated use, for instance if you want to use existing SExtractor input files,
or reveal the output files. This module

- is based on `astropy` (instead of astroasciidata)
- uses standard `logging` (no prints)
- uses `tempfile` to hide all input and output files (if you don't care about them)
- has some convenience functionality to use SExtractor's `ASSOC` process

Installation
------------

Clone or download this repository, `cd` into it, and 

```
python setup.py install --user
```

For a system-wide install, remove the `--user`. More information about how to install python modules with distutil can be found [here](https://docs.python.org/2/install/index.html#install-index).

Requirements are [astropy](http://www.astropy.org/) 0.4.2 or later, and obviously [SExtractor](http://www.astromatic.net/software/sextractor) (`sewpy` does not require a specific version).

Documentation
-------------

The API is documented (using Sphinx) on the github project page: http://megalut.github.io/sewpy/


