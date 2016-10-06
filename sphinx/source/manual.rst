User manual
===========

General philosophy (to be rewritten)
------------------------------------

* If you do not specify a workdir, I'll write all my required internal files somewhere in /tmp,
  and you don't have to bother about this (this is done using the tempfile module).
* ``params`` (a list) refers to the features that you want SExtractor to measure
  (e.g., settings you find in "default.param").
* ``config`` (a dict) refers to the settings (e.g., stuff you find in "default.sex").
* A SEW instance can well be reused for different images
  that you want to analyse with the same params but a different config.
  Indeed you usually don't want to change params from image to image, but you might have to change
  the config (e.g., gain, seeing, ...).
* In the params list, you have to specify all the parameters that you want to be measured.
* On the other hand, in the config dict, you only have to give those settings that deviate from
  the *default*! We take as *default* the output of "sextractor -d" (if not told otherwise).
* When repeatedly calling run(), we avoid writing the SExtractor input files to disk over and over again.
  Instead, param is written only once, and config settings are passed as command line arguments to
  the SExtractor executable, superseding the default config.
  So to change config from image to image, simply edit se.config between calls of run().
* There is special helper functionality for using ASSOC (see note below)

Logging
-------

Sewpy uses the logging module.
To see a detailed log of what is going on, insert this into your script::

	import logging
	logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)


The ASSOC helper
----------------

The ASSOC helper assists you in measuring galaxies from an existing input catalog,
instead of just making a new catalog of all sources. In summary, you pass an existing input
catalog, and you'll get this same catalog as output, but with the new columns
corresponding to the SExtractor params appended.

To use the ASSOC helper:

* Add ``VECTOR_ASSOC(3)`` to your params (The following is valid only for Astropy v1.1.2 and prior: at the beginning, not at the end, of the params list).
* Add for instance ``{"ASSOC_RADIUS":10.0, "ASSOC_TYPE":"NEAREST"}`` to your config.
  These values are the defaults used if you don't specify anything.
* Give the relevant arguments (``assoc_cat, assoc_xname, assoc_yname``) when calling.
		   
The output will contain an astropy table, with the same rows as ``assoc_cat``, but 
to which the new SExtractor columns will be appended.
Those SExtractor columns might be **masked** columns (leading to a masked table),
as some of your sources might not have been found by SExtractor.
Note that the attribute mytable.masked tells you if an astropy table "mytable" is masked.
To make it even more foolproof, I systematically add a boolean column named
``prefix + "assoc_flag"``. True means that the source was found.


		

