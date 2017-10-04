Source Extractor Wrapper for Python
===================================

The tiny ``sewpy`` module lets you run `SExtractor <http://www.astromatic.net/software/sextractor>`_ as if it would all be native python...

.. code-block:: python 
	
	import sewpy
	sew = sewpy.SEW(params=["X_IMAGE", "Y_IMAGE", "FLUX_RADIUS(3)", "FLAGS"],
		config={"DETECT_MINAREA":10, "PHOT_FLUXFRAC":"0.3, 0.5, 0.8"})
	out = sew("myimage.fits")
	print out["table"] # this is an astropy table.

... but also allows for a more sophisticated use, for instance if you want to use existing SExtractor input files.

Why `yet <https://pypi.python.org/pypi/pysex/>`_ `another <https://gitorious.org/pysextractor>`_ SExtractor wrapper ? Because we needed one that:

* is based on `astropy <http://www.astropy.org>`_ (in particular ``astropy.table``),
* works with both python 2 and 3,
* uses standard ``logging`` (and also logs SExtractor's stdout & stderr to file),
* uses ``tempfile`` to hide all input and output files, *except if you want to see them*
* has some convenience functionality to use SExtractor's ``ASSOC`` process (give me an input catalog, and I append columns with SExtractor measurements to it).

The **demos** in the ``examples`` directory can be run without installing sewpy, and provide a quick overview. 

Installation (in short)
-----------------------

.. code-block:: bash
	
	python setup.py install --user
	

Documentation
-------------

To learn more about how to install and use ``sewpy``, proceed to its `documentation <http://sewpy.readthedocs.org>`_.





