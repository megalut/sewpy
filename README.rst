Source Extractor Wrapper for Python
===================================

The tiny ``sewpy`` module let's you run `SExtractor <http://www.astromatic.net/software/sextractor>`_ as if it would all be native python...

.. code-block:: python 
	
	import sewpy
	sew = sewpy.SEW(params=["X_IMAGE", "Y_IMAGE", "FLUX_RADIUS(3)", "FLAGS"],
		config={"DETECT_MINAREA":10, "PHOT_FLUXFRAC":"0.3, 0.5, 0.8"})
	out = sew("myimage.fits")
	print out["table"] # this is an astropy table.

... but also allows for a more sophisticated use, for instance if you want to use existing SExtractor input files, or reveal the output files. The ``sewpy`` module

* is based on `astropy <http://www.astropy.org>`_
* uses standard ``logging``
* uses ``tempfile`` to hide all input and output files (except if you care about them)
* has some convenience functionality to use SExtractor's ``ASSOC`` process

Documentation
-------------

To learn more about how to install and use ``sewpy``, proceed to its `documentation <http://sewpy.readthedocs.org>`_.**





