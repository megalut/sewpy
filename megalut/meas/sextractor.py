"""
Yet another attempt to build a reusable and transparent SExtractor wrapper, this time:
- using only astropy
- using logging
- with support for ASSOC
- not only with MegaLUT in mind

"""



import os
import astropy
import subprocess
import tempfile

import logging
logger = logging.getLogger(__name__)


defaultparams = ["XWIN_IMAGE", "YWIN_IMAGE", "AWIN_IMAGE", "BWIN_IMAGE", "THETAWIN_IMAGE", "BACKGROUND", "FLUX_AUTO"]


class SExtractor():
	"""
	Holds together all the configuration and settings.
	"""
	
	def __init__(self, workdir = None, sexpath="sex"):
		"""
		:param workdir: where I'll write my files. If None, I use a temporary directory.
		"""
	
		if workdir is not None:
			self.workdir = workdir
			self.tmp = False
		else:
			self.workdir = tempfile.mkdtemp(prefix='sextractor_dot_py_workdir_')
			self.tmp = True
		
		self.sexpath = sexpath
		
		self.params = []
		self.conf = {}
		
		
		self.set_params(defaultparams)
	
	def __str__(self):
		return "SExtractor in %s" % (self.workdir)
	
	
	def set_params(self, params):
		"""
		
		"""
		self.params = params
	
		
		

	def run(self, imgfilepath, log=True):
		"""
		
		"""


		p = subprocess.Popen([self.sexpath, imgfilepath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		
	def delete(self):
		"""
		Removes the temporary working dir (but not a problem if you do not use this).
		"""
		if self.tmp:
			if os.path.isdir(self.workdir):
				os.remove(self.workdir)


"""
# We write the output into a log file...
logfile = open(exp + ".log", "w")
logfile.write("======= This is the output of createPSFcube =========\n")
out, err = p.communicate()
logfile.write(out)
logfile.write("======= And this is the stderr =========\n")
logfile.write(err)
logfile.close()
# And show errors as log.warnings :
if err != "":
logger.warning(err)

"""
