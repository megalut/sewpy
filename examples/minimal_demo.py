# The following two lines let this script find sewpy, even if sewpy has not yet been installed.
import sys, os
sys.path.insert(0, os.path.abspath('../'))

# Once sewpy is installed, you can start your script from here...

import sewpy

sew = sewpy.SEW(
		params=["X_IMAGE", "Y_IMAGE", "FLUX_RADIUS(3)", "FLAGS"],
		config={"DETECT_MINAREA":10, "PHOT_FLUXFRAC":"0.3, 0.5, 0.8"}
	)

out = sew("image.fits")

print out["table"] # This is an astropy table.

