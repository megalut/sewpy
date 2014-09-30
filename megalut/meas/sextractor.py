"""
Yet another attempt to build a reusable and transparent SExtractor wrapper, this time

- using only astropy (no more astroasciidata or pyfits)
- using logging
- using tempfile (if needed)
- with support for ASSOC
- not only with MegaLUT in mind

Several SExtractor wrappers for python can be found online, and the code for this module mixes elements inspired by:

- previous MegaLUT, alipy, and cosmouline implementations
- pysex from Nicolas Cantale
- sextractor.py by Laurent Le Guillou (from the "forgotten" COSMOGRAIL pipeline)
- pysextractor by Nicolas Gruel


Minimal example::

 from sextractor import SExtractor
 se = SExtractor()
 cat = se.run("myimage.fits")
 print cat

The philosophy is the following:

- A SExtractor instance ("se" in the example above) is meant to be reused for different images
  that you want to analyse with roughly the same settings.
- These "permanent" settings can be given via the constructor.
- When repeatedly calling run(), we avoid rewriting the SExtractor configuration files over and over again.
  Instead, "non-permanent" settings (say the gain of each image, or the output catalog name) are passed as command line arguments.
- Options should let you keep the catalog files, output images etc etc.


To do:

- clarify the above philosophy, find out how we want to use it
- fix Exceptions
- better handling of what files to keep and where
- Do we need anything more to use ASSOC
- option to log sextractor bla bla
- delete / cleanup afterwards

"""


import os
import astropy
import subprocess
import tempfile
import re
from datetime import datetime

import logging
logger = logging.getLogger(__name__)


defaultparams = ["XWIN_IMAGE", "YWIN_IMAGE", "AWIN_IMAGE", "BWIN_IMAGE", "THETAWIN_IMAGE", "BACKGROUND", "FLUX_AUTO"]
defaultconfig = {}


class SExtractorError(Exception):
    pass


class SExtractor():
	"""
	Holds together all the configuration and settings.
	"""
	
	def __init__(self, workdir = None, sexpath="sex", params=None, config=None):
		"""
		All arguments have default values and are optional.
		
		:param workdir: where I'll write my *internal* files. If None, I create a temporary directory myself, usually in /tmp.
		:param sexpath: path to the sextractor executable (e.g., "sex" or "sextractor", if in your PATH)
		:param params: the parameters you want SExtractor to measure (i.e., what you usually write in the "default.param" file).
		:param config: special settings (overwrite the defaults from the usual "default.sex" file).
		
		"""
	
		if workdir is not None:
			self.workdir = workdir
			self.tmp = False
		else:
			self.workdir = tempfile.mkdtemp(prefix='sextractor_dot_py_workdir_')
			self.tmp = True
		
		self.sexpath = sexpath
		
		if params == None:
			self.params = defaultparams
		else:
			self.params = params
			
		if config == None:
			self.config = defaultconfig
		else:
			self.config = config
		
		self.catfilename = "test.cat" # That's the official SExtractor default.
		
	
	def get_version(self):
		"""
		To find the SExtractor version, we call it without arguments and parse the stdout.
		"""
		p = subprocess.Popen([self.sexpath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate()
		version_match = re.search("[Vv]ersion ([0-9\.])+", err)
		assert version_match is not False
		version = str(version_match.group()[8:])
		assert len(version) != 0
		return version
		
	
	def __str__(self):
		return "SExtractor in %s" % (self.workdir)
		
	
	def get_params_filepath(self):
		return os.path.join(self.workdir, "params.txt")
		
	def get_config_filepath(self):
		return os.path.join(self.workdir, "config.txt")

	def get_conv_filepath(self):
		return os.path.join(self.workdir, "conv.txt")
	
	def get_cat_filepath(self):
		return os.path.join(self.workdir, self.catfilename)
	
	def _write_params(self):
		"""
		Writes the parameters to the file
		"""
		f = open(self.get_params_filepath(), 'w')
		f.write("\n".join(self.params))
		f.write("\n")
		f.close()
		logger.debug("Wrote %s" % (self.get_params_filepath()))


	def _write_default_config(self):
		p = subprocess.Popen([self.sexpath, "-d"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate()
		if err != "":
			logger.warning("Ouch, SExtractor complains :")
			logger.warning(err)
		f = open(self.get_config_filepath(), 'w')
		f.write(out)
		f.close()
		logger.debug("Wrote %s" % (self.get_config_filepath()))


	def _write_default_conv(self):
		f = open(self.get_conv_filepath(), 'w')
		f.write("""CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1""")
		f.close()
		logger.debug("Wrote %s" % (self.get_conv_filepath()))

		

	def run(self, imgfilepath, catfilepath=None, log=True, assoccat=None, readcat=True):
		"""
		:param imgfilepath: Path to the FITS image I should run on
		:param catfilepath: Path of the catalog I should write. If None (default), I won't keep this catalog around.
		:param readcat: By default I read the SExtractor catalog and return it as an astropy table. If set to False, I skip this step.
		
		:returns: the astropy table of the output catalog (unless readcat is False).
		"""

		logging.info("Running SExtractor on %s..." % imgfilepath)
		starttime = datetime.now()
		
		self._write_params()
		self.config["PARAMETERS_NAME"] = self.get_params_filepath()
		
		self._write_default_config()
		
		self.config["CATALOG_NAME"] = self.get_cat_filepath()
		
		self.config["FILTER_NAME"] = self.get_conv_filepath()
		self._write_default_conv()
		
		

		popencmd = [self.sexpath, imgfilepath, "-c", self.get_config_filepath()]
		for (key, value) in self.config.items():
			popencmd.append("-"+str(key))
			popencmd.append(str(value))
			
		logging.debug("Running with command %s..." % (popencmd))
		p = subprocess.Popen(popencmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate()
		
		#print out
		"""
		if err != "":
			logger.warning("Ouch, SExtractor complains :")
			logger.warning(err)
		"""
		
		endtime = datetime.now()
		logger.info("Done, it took %.2f seconds." % ((endtime - starttime).total_seconds()))

		if readcat:
			
			cat = astropy.table.Table.read(self.get_cat_filepath(), format="ascii.sextractor")
			return cat
		
		
	def delete(self):
		"""
		Removes the temporary working dir (but not a problem if you do not use this).
		"""
		if self.tmp:
			if os.path.isdir(self.workdir):
				os.remove(self.workdir)


