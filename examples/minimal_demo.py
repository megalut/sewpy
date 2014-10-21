# The following two lines let this script find sewpy, even if sewpy has not yet been installed.
import sys, os
sys.path.insert(0, os.path.abspath('../'))

# Once sewpy is installed, you can start your script from here...

# The next two lines are temporary and will be moved from the minimal demo to the second minimal demo...
import logging
logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)

import sewpy

sew = sewpy.SEW(
		params=["X_IMAGE", "Y_IMAGE", "FLUX_APER(3)", "FLAGS"],
		config={"DETECT_MINAREA":10, "PHOT_APERTURES":"5, 10, 20"},
		sexpath="sex"
	)
# By default, this assumes that SExtractor can be called as "sex"
# If this is not the case, or if the executable is not in your path,
# specify the path by adding the argument sexpath="/path/to/sextractor"
# to the above instantiation.
	
out = sew("image.fits")

print out["table"] # This is an astropy table.

