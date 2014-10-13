# Source Extractor Wrapper for Python

The tiny `sewpy` module let's you run [SExtractor](http://www.astromatic.net/software/sextractor) as if it would all be native python... 

```python 
import sewpy
sew = sewpy.SEW(params=["X_IMAGE", "Y_IMAGE", "FLUX_RADIUS(3)", "FLAGS"])
out = sew("myimage.fits", config={"DETECT_MINAREA":10, "PHOT_FLUXFRAC":"0.3, 0.5, 0.8"})
print out["table"] # this is an astropy table.
```

... but also allows for a more sophisticated use, for instance if you want to use existing SExtractor input files,
or reveal the output files. The module

- is based on `astropy` (instead of astroasciidata)
- uses standard `logging` (no prints)
- uses `tempfile` to hide all input and output files (if you don't care about them)
- has some convenience functionality to use SExtractor's `ASSOC` process

To learn more about how to use ``sewpy``, the best is to have a look at the scripts in the ``examples`` directory. You can also view the full API documentation, hosted at [http://sewpy.readthedocs.org](sewpy.readthedocs.org).


## Installation

Clone or download this repository, `cd` into it, and 

```
python setup.py install --user
```

For a system-wide install, remove the `--user`. More information about how to install python modules with distutil can be found [here](https://docs.python.org/2/install/index.html#install-index).

Requirements are [astropy](http://www.astropy.org/) 0.4.2 or later, and obviously [SExtractor](http://www.astromatic.net/software/sextractor) (`sewpy` does not require a specific version).

## Documentation

### General philosophy (to be rewritten)

- If you do not specify a workdir, I'll write all my required internal files somewhere in /tmp,
  and you don't have to bother about this (this is done using the tempfile module).
- "params" (a list) refers to the features that you want SExtractor to measure
  (e.g., settings you find in "default.param").
- "config" (a dict) refers to the settings (e.g., stuff you find in "default.sex").
- A SExtractor instance ("se" in the example above) can well be reused for different images
  that you want to analyse with the same params but a different config.
  Indeed you usually don't want to change params from image to image, but you might have to change
  the config (e.g., gain, seeing, ...).
- In the params list, you have to specify all the parameters that you want to be measured.
- On the other hand, in the config dict, you only have to give those settings that deviate from
  the *default*! We take as *default* the output of "sextractor -d" (if not told otherwise).
- When repeatedly calling run(), we avoid writing the SExtractor input files to disk over and over again.
  Instead, param is written only once, and config settings are passed as command line arguments to
  the SExtractor executable, superseding the default config.
  So to change config from image to image, simply edit se.config between calls of run().
- There is special helper functionality for using ASSOC (see note below)


### The ASSOC helper

The ASSOC helper assists you in measuring galaxies from an existing input catalog,
instead of just making a new catalog of all sources. In summary, you pass an existing input
catalog to run(), and you'll get this same catalog as output, but with the new columns
corresponding to the SExtractor params appended.

To use the ASSOC helper:
- Add `VECTOR_ASSOC(3)` to your params (at the beginning, not at the end, of the params list).
- Add for instance `{"ASSOC_RADIUS":10.0, "ASSOC_TYPE":"NEAREST"}` to your config.
  These values are the defaults used if you don't specify anything.
- Give the relevant arguments (`assoc_cat, assoc_xname, assoc_yname`) when calling.
		   
The output will contain an astropy table, with the same rows as `assoc_cat`, but 
to which the new SExtractor columns will be appended.
Those SExtractor columns might be **masked** columns (leading to a masked table),
as some of your sources might not have been found by SExtractor.
Note that the attribute mytable.masked tells you if an astropy table "mytable" is masked.
To make it even more foolproof, I systematically add a boolean column named
`prefix + "assoc_flag"`. True means that the source was found.

### Generating the API documentation

See `sphinx/source/index.rst`

		
## Recent improvements (latest on top):

- better verbosity about masked output of ASSOC procedure
- ASSOC helper implemented
- run() now returns a dict containing several objects, such as the output astropy table, catfilepath, workdir, and logfilepath.
- now also works with vector parameters such as MAG_APER(4)
- possibility to "nice" SExtractor
- a log file is written for every run() if not told otherwise
- filenames change according to FITS image file name where required
- but you can also pass an "imgname" argument to run, and this will be used instead.
- params and config files are written only once, as discussed
- appropriate warnings and behaviour when a workdir already exists, or when you rerun on the same file
- possibility to use existing param / config / conv / nnw files
- run() returns either the catalog, or the filepath to the catalog


## To do:

- move "config" to run ?
- check that all masked columns of ASSOC do indeed share the same mask.
- implement _check_config()
- better detection of SExtractor failures
- implement raising Exceptions when SExtractor fails
- implement CHECK IMAGE "helper" ?
- give access to several conv and nnw settings (if needed)


## Credits

The code of the present module contains elements inspired by
- previous MegaLUT, alipy, and cosmouline implementations
- pysex from Nicolas Cantale
- sextractor.py by Laurent Le Guillou (from the "forgotten" COSMOGRAIL pipeline)
- pysextractor by Nicolas Gruel


