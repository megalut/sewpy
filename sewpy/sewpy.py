"""
sewpy: Source Extractor Wrapper for Python


Recent improvements (latest on top):

- new loglevel option to adjust sewpy's overall "verbosity" on instantiation.
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


To do:

- move "config" to run ?
- check that all masked columns of ASSOC do indeed share the same mask.
- implement _check_config()
- better detection of SExtractor failures
- implement raising Exceptions when SExtractor fails
- implement CHECK IMAGE "helper" ?
- give access to several conv and nnw settings (if needed)


"""


import os
import astropy
import astropy.table
import subprocess
import tempfile
import re
import copy
from datetime import datetime
import numpy as np
from astropy.io import fits

import logging
logger = logging.getLogger(__name__)


defaultparams = ["XWIN_IMAGE", "YWIN_IMAGE", "AWIN_IMAGE", "BWIN_IMAGE", "THETAWIN_IMAGE", "BACKGROUND",
                 "FLUX_AUTO"]
defaultconfig = {}


class SEW():
	"""
	Holds together all the settings to run SExtractor executable on one or several images.
	"""
	
	def __init__(self, workdir=None, sexpath="sex", params=None, config=None, configfilepath=None, nice=None, loglevel=None):
		"""
		All arguments have default values and are optional.
		
		:param workdir: where I'll write my files. Specify this (e.g., "test") if you care about the
                        output files.
			If None, I create a unique temporary directory myself, usually in /tmp.
			
		:param sexpath: path to the sextractor executable (e.g., "sex" or "sextractor", if in your PATH)
		:param params: the parameters you want SExtractor to measure (i.e., what you would write in the
                        "default.param" file)
		:type params: list of strings
		:param config: config settings that will supersede the default config (e.g., what you would
                        change in the "default.sex" file)
		:type config: dict
		:param configfilepath: specify this if you want me to use an existing SExtractor config file as
                        "default" (instead of the sextractor -d one)
		:param nice: niceness with which I should run SExtractor. Use e.g. ``19`` for set lowest priority.
		:type nice: int
		
		:param loglevel: verbosity, e.g. the python-level logging threshold for the sewpy module logger.
			For example, set this to "WARNING" and sewpy will no longer log simple INFOs.
			Choices are "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL". 
			To disable logging, set ``loglevel="CRITICAL"``
		:type loglevel: string or int or logging.level...
		
		
		To use an existing SExtractor param-, conv-, or nnw-file, simply specify these in the config
                dict, using the appropriate SExtractor keys (PARAMETERS_NAME, FILTER_NAME, ...)
		
		.. warning:: When using *vector*-type params resulting in multiple columns (such as "FLUX_RADIUS(3)"
			in the example above), do not put these in the last position of the params list, otherwise astropy
			fails reading the catalog! This is probably due to the fact that the SExtractor header doesn't give
			a hint that multiple columns are expected when a vector-type param comes last. A workaround would be
			way too complicated.

		
		
		"""

		# We start by setting the log "verbosity":
		if loglevel != None:
			logger.setLevel(loglevel)
		
		
		# We set up the trivial things:
		
		self.sexpath = sexpath
		self.configfilepath = configfilepath
		self.nice = nice
		
		logger.info("SExtractor version is %s" % (self.get_version()))
		
		# ... and the workdir
		
		if workdir is not None:
			self.workdir = workdir
			self.tmp = False
			if os.path.isdir(workdir):
				#logger.warning("SExtractor workdir '%s' exists, be careful! I will (maybe silently) delete or overwrite stuff." % (workdir))	
				pass
			else:
				logger.info("Making new SExtractor workdir '%s'..." % (workdir))
				os.makedirs(workdir)
		else:
			self.workdir = tempfile.mkdtemp(prefix='sewpy_workdir_')
			self.tmp = True
		

		#self._clean_workdir()
		# No, don't clean it ! This is an obvious race conditions when several processes use the same workdir !
		# Commenting this is just a quick fix, we need to clean this up.
		
		# ... and the params:
		
		if params == None:
			self.params = defaultparams
		else:
			self.params = params
		self._check_params()
			
		
		# ... and the config:
		
		if config == None:
			self.config = defaultconfig
		else:
			self.config = config
		
		self._set_instance_config() # Adds some fixed stuff to self.config
		self._check_config()
	
	
	def get_version(self):
		"""
		To find the SExtractor version, we call it without arguments and parse the stdout.
		
		:returns: a string (e.g. '2.4.4')
		"""
		try:
			p = subprocess.Popen([self.sexpath], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		except:
			raise RuntimeError("Could not run SExtractor. Is the path '%s' correct ? If not, specify sexpath='/path/to/sextractor'" % self.sexpath)
		out, err = p.communicate()
		version_match = re.search("[Vv]ersion ([0-9\.])+", err.decode(encoding='UTF-8'))
		if version_match is False:
			raise RuntimeError("Could not determine SExctractor version, check the output of running '%s'" % (self.sexpath))
		version = str(version_match.group()[8:])
		assert len(version) != 0
		return version
	
	
	def __str__(self):
		"""
		A string summary representing the instance
		"""
		return "'SEW object with workdir %s'" % (self.workdir)
	
	
	def _check_params(self):
		"""
		Compares the params to a list of known params, and spits out a useful warning if
		something seems fishy.
		"""
		strange_param_helper = False
		for param in self.params:
		
			# It could be that the param encapsulates several values (e.g., "FLUX_RADIUS(10)")
			# So we have to dissect this
			match = re.compile("(\w*)\(\d*\)").match(param)
			if match:
				cleanparam = match.group(1)
			else:
				cleanparam = param
				
			if cleanparam not in self.fullparamlist:
				logger.warning("Parameter '%s' seems strange and might be unknown to SExtractor" \
                                                   % (param))
				strange_param_helper = True
				
		if strange_param_helper:
			logger.warning("Known parameters are: %s" % (self.fullparamtxt))
			
	
	def _check_config(self):
		"""
		Not yet implemented
		"""
		pass

		
	def _set_instance_config(self):
		"""
		Sets config parameters that remain fixed for this instance.
		Called by __init__(). If needed, you could still mess with this config after __init__() has run.
		"""
		
		if "PARAMETERS_NAME" in self.config.keys():
			logger.info("You specified your own PARAMETERS_NAME, I will use it.")
		else:
			self.config["PARAMETERS_NAME"] = self._get_params_filepath()
		
		if "FILTER_NAME" in self.config.keys():
			logger.info("You specified your own FILTER_NAME, I will use it.")
		else:
			self.config["FILTER_NAME"] = self._get_conv_filepath()
		
		
		if "CATALOG_NAME" in self.config.keys():
			logger.warning("You specified your own CATALOG_NAME, but I will *NOT* use it !")
			del self.config["CATALOG_NAME"]

		if "PSF_NAME" in self.config.keys():
			logger.info("You specified your own PSF_NAME, I will use it.")
		else:
			self.config["PSF_NAME"] = self._get_psf_filepath()		

	def _get_params_filepath(self):
		"""
		Stays the same for a given instance.
		"""
		return os.path.join(self.workdir, "params.txt")
		
	def _get_config_filepath(self):
		"""
		Idem, stays the same for a given instance.
		Might return the non-default configfilepath, if set.
		"""
		if self.configfilepath is None:
			return os.path.join(self.workdir, "config.txt")
		else:
			return self.configfilepath

	def _get_conv_filepath(self):
		"""
		Stays the same for a given instance.
		"""
		return os.path.join(self.workdir, "conv.txt")
		
	def _get_psf_filepath(self):
		"""
		Stays the same for a given instance.
		"""
		return os.path.join(self.workdir, "default.psf")
	
	def _get_cat_filepath(self, imgname):
		"""
		This changes from image to image
		"""
		return os.path.join(self.workdir, imgname + ".cat.txt")
	
	def _get_assoc_filepath(self, imgname):
		"""
		Changes from image to image
		"""
		return os.path.join(self.workdir, imgname + ".assoc.txt")
	
	def _get_log_filepath(self, imgname):
		"""
		Changes from image to image
		"""
		return os.path.join(self.workdir, imgname + ".log.txt")
	
	
	def _write_params(self, force=False):
		"""
		Writes the parameters to the file, if needed.
		
		:param force: if True, I overwrite any existing file.
		"""
		if force or not os.path.exists(self._get_params_filepath()):
			f = open(self._get_params_filepath(), 'w')
			f.write("\n".join(self.params))
			f.write("\n")
			f.close()
			logger.debug("Wrote %s" % (self._get_params_filepath()))
		else:
			logger.debug("The params file already exists, I don't overwrite it.")


	def _write_default_config(self, force=False):
		
		"""
		Writes the *default* config file, if needed.
		I don't write this file if a specific config file is set.
		
		:param force: if True, I overwrite any existing file.
		"""
		
		if self.configfilepath is not None:
			logger.debug("You use the existing config file %s, I don't have to write one." % \
                                         (self._get_config_filepath()))
			return
		
		if force or not os.path.exists(self._get_config_filepath()):	
			p = subprocess.Popen([self.sexpath, "-dd"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			if err != "":
				logger.warning("Ouch, SExtractor complains :")
				logger.warning(err)
			f = open(self._get_config_filepath(), 'w')
			f.write(out.decode(encoding='UTF-8'))
			f.close()
			logger.debug("Wrote %s" % (self._get_config_filepath()))
		else:
			logger.debug("Default config file already exists, I don't overwrite it.")



	def _write_default_conv(self):
		"""
		Writes the default convolution matrix, if needed.
		"""
		
		if not os.path.exists(self._get_conv_filepath()):	
			f = open(self._get_conv_filepath(), 'w')
			f.write("""CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1""")
			f.close()
			logger.debug("Wrote %s" % (self._get_conv_filepath()))
		else:
			logger.debug("Default conv file already exists, I don't overwrite it.")

	def _write_default_psf(self):
		"""Writes the default psf file, if needed.
		"""
		if not os.path.exists(self._get_psf_filepath()):
			arr= np.array([[[  2.40381883e-06,   7.60069597e-05,   8.54981554e-05,
           7.54465946e-05,   8.13425067e-05,   8.50538563e-05,
           7.86255987e-05,   3.32223935e-05,  -2.53305316e-05,
          -2.43271165e-06,   1.11580441e-04,   1.85113196e-04,
           1.11781279e-04,  -5.22066912e-05,  -1.48545674e-04,
          -6.13458324e-05,   7.43288838e-05,   2.76304218e-05,
          -1.35628201e-04,  -1.57351358e-04,  -5.41317095e-05,
          -3.80859383e-05,  -8.01880524e-05,   1.93700293e-06,
           1.28125030e-04,   9.64317587e-05,  -1.00769712e-05,
          -3.07009032e-05,   2.51779438e-05,   7.89448313e-05,
           1.31651745e-04],
        [  1.07151063e-04,   7.80015980e-05,   5.20944668e-05,
           1.96801047e-05,   2.58915697e-05,   4.22252306e-05,
           2.61278237e-05,   1.09607536e-05,   2.22668241e-05,
           1.83793582e-05,   1.69202158e-05,   7.20969401e-05,
           1.13053175e-04,   2.51594101e-05,  -6.95991839e-05,
           4.80339113e-05,   2.57156673e-04,   2.34465711e-04,
           3.90260866e-05,   3.77108845e-05,   1.98102425e-04,
           1.63195917e-04,  -5.02058392e-05,  -9.12948817e-05,
           8.09195080e-06,  -3.04424921e-05,  -1.54503243e-04,
          -1.23093574e-04,   3.82040671e-05,   1.51542059e-04,
           1.55646892e-04],
        [  1.17171738e-04,   9.47605004e-05,   8.96408383e-05,
           1.84595974e-05,  -6.76616692e-05,  -6.04558227e-05,
           3.12718876e-05,   1.18101605e-04,   1.40572258e-04,
           7.31113250e-05,  -2.21763585e-05,  -1.00581055e-05,
           9.72477792e-05,   1.23102203e-04,   7.73502907e-05,
           1.37092138e-04,   2.56634929e-04,   2.35879590e-04,
           1.30117784e-04,   1.90375547e-04,   3.48268659e-04,
           2.85394461e-04,   1.23960363e-05,  -1.02018013e-04,
          -1.04167375e-05,  -6.52167273e-06,  -1.33219000e-04,
          -1.12960566e-04,   8.09253615e-05,   2.08761281e-04,
           1.66711587e-04],
        [  1.32234691e-05,   1.02787519e-04,   1.83809068e-04,
           9.58435703e-05,  -9.55243813e-05,  -9.70988549e-05,
           1.35439012e-04,   2.94848898e-04,   2.40708367e-04,
           1.22713420e-04,   5.07345976e-05,   6.42455634e-05,
           1.56476934e-04,   1.97447749e-04,   1.18205535e-04,
           4.94819542e-05,   6.14344171e-05,   9.58741002e-05,
           1.28766158e-04,   1.94652850e-04,   2.41393209e-04,
           1.60364230e-04,  -1.41999890e-05,  -6.13596058e-05,
           6.79572622e-05,   1.21097153e-04,  -9.76944037e-09,
          -5.49086435e-05,   5.58508291e-05,   1.36878050e-04,
           1.21550351e-04],
        [  5.26016038e-07,   1.24271712e-04,   2.35970729e-04,
           1.74356275e-04,  -1.32342866e-05,  -5.57084786e-05,
           1.29754219e-04,   2.72921199e-04,   2.43686882e-04,
           1.98139867e-04,   1.90322971e-04,   1.92822670e-04,
           2.19771260e-04,   1.89034996e-04,   2.74552949e-05,
          -1.07212341e-04,  -5.95865968e-05,   1.05343766e-04,
           2.24970223e-04,   2.19185677e-04,   1.25688792e-04,
           4.09510940e-05,   5.42922544e-06,   3.56244018e-05,
           1.26195810e-04,   1.44066478e-04,   3.65264459e-05,
          -3.75594900e-05,   6.59612351e-06,   4.50785592e-05,
           7.68707177e-05],
        [  9.01028179e-05,   1.57428833e-04,   2.07729943e-04,
           1.83902201e-04,   8.59787106e-05,  -2.93249741e-05,
          -4.01443904e-05,   8.65905531e-05,   2.74296646e-04,
           3.97606462e-04,   3.76625831e-04,   3.10488453e-04,
           3.03925684e-04,   2.50657264e-04,   8.78350329e-05,
          -1.62432534e-05,   8.49524658e-05,   2.89930234e-04,
           4.05957981e-04,   3.60714155e-04,   2.31550643e-04,
           1.86596124e-04,   2.22587347e-04,   1.94238193e-04,
           9.77661039e-05,   1.30580702e-05,  -3.79030898e-05,
          -2.51111196e-05,   5.09598831e-05,   1.17862386e-04,
           1.58784678e-04],
        [  1.26586805e-04,   1.47764193e-04,   1.21227618e-04,
           1.07567241e-04,   9.73561691e-05,  -2.62960634e-06,
          -7.90764097e-05,   8.23809023e-05,   4.21559729e-04,
           5.81062224e-04,   4.64305049e-04,   3.95985029e-04,
           5.04385040e-04,   5.14040526e-04,   3.40078899e-04,
           2.08895421e-04,   2.48318480e-04,   3.65095621e-04,
           4.52472566e-04,   4.69565508e-04,   4.18782613e-04,
           4.12627007e-04,   4.26682149e-04,   3.09842784e-04,
           1.15190043e-04,   1.26203922e-05,   1.67655326e-05,
           9.60616308e-05,   2.24732503e-04,   3.24832916e-04,
           3.07584094e-04],
        [  8.70921285e-05,   8.99088554e-05,   4.57445130e-05,
           1.47284136e-05,   2.48302695e-05,   3.75261916e-05,
           1.04437058e-04,   3.04165005e-04,   5.74388076e-04,
           6.31482864e-04,   5.24375762e-04,   6.48113259e-04,
           1.03455526e-03,   1.23596319e-03,   1.10962230e-03,
           9.19762475e-04,   8.06013879e-04,   7.29060732e-04,
           6.80224970e-04,   6.54205156e-04,   5.85820526e-04,
           5.15840831e-04,   4.28209722e-04,   3.26885085e-04,
           2.66392512e-04,   2.65364826e-04,   2.54026818e-04,
           2.83679867e-04,   3.82139726e-04,   4.34142858e-04,
           3.40681698e-04],
        [  2.67705618e-05,   3.52060233e-05,   5.01444010e-05,
           2.07116755e-05,  -1.47125493e-05,   2.03489180e-05,
           1.74667264e-04,   3.91519163e-04,   6.15471101e-04,
           7.03603029e-04,   7.29316904e-04,   9.43472434e-04,
           1.35301880e-03,   1.60744193e-03,   1.57620630e-03,
           1.42798119e-03,   1.28054479e-03,   1.16847351e-03,
           1.07864384e-03,   9.50956077e-04,   7.26292899e-04,
           5.19602501e-04,   3.67075088e-04,   3.11377604e-04,
           3.44564585e-04,   4.15781105e-04,   4.12696565e-04,
           3.71932460e-04,   3.34744342e-04,   2.89331714e-04,
           2.17575871e-04],
        [ -6.19701359e-06,   2.79479227e-05,   1.08629363e-04,
           1.10446766e-04,   4.28918393e-05,  -1.59790234e-05,
           6.37421181e-05,   2.94988626e-04,   6.72983995e-04,
           9.75804869e-04,   1.23918452e-03,   1.58980582e-03,
           2.15687836e-03,   2.84912600e-03,   3.44874267e-03,
           3.64848366e-03,   3.35466117e-03,   2.84754485e-03,
           2.30689673e-03,   1.70514034e-03,   1.04716129e-03,
           6.41009770e-04,   4.67646430e-04,   4.08837222e-04,
           2.92016397e-04,   2.86833878e-04,   3.29357979e-04,
           2.67972850e-04,   1.00064397e-04,   3.63361505e-05,
           9.96371964e-05],
        [ -5.67882466e-07,   5.08060548e-05,   1.30551591e-04,
           1.57827846e-04,   1.27582811e-04,   4.91817409e-05,
           9.40121536e-05,   3.15737299e-04,   7.65093428e-04,
           1.24707678e-03,   1.83331280e-03,   2.57003610e-03,
           3.45921540e-03,   4.29388555e-03,   4.88529447e-03,
           5.04447240e-03,   4.70967358e-03,   4.11674846e-03,
           3.41683673e-03,   2.55684974e-03,   1.53652672e-03,
           8.98711965e-04,   6.89776614e-04,   6.15979312e-04,
           3.43337771e-04,   1.73978668e-04,   1.21991347e-04,
           3.66753629e-05,  -8.58800995e-05,  -3.52677125e-05,
           1.08530170e-04],
        [  3.37872116e-05,   6.98544100e-05,   1.03699349e-04,
           1.31374065e-04,   1.63710123e-04,   1.81806943e-04,
           3.33123491e-04,   5.07342280e-04,   9.24890337e-04,
           1.47657737e-03,   2.67444947e-03,   4.68453299e-03,
           7.47570582e-03,   1.05615249e-02,   1.58804134e-02,
           1.49398511e-02,   1.19190756e-02,   9.18036141e-03,
           6.40049716e-03,   4.11758851e-03,   2.35498301e-03,
           1.46076444e-03,   9.99144861e-04,   8.69806274e-04,
           5.58527070e-04,   2.71378696e-04,   1.69032792e-05,
          -1.02712445e-04,  -7.02594843e-05,   7.04579725e-05,
           1.29349166e-04],
        [  1.50863212e-04,   1.31067063e-04,   1.07944987e-04,
           1.11346009e-04,   1.84717428e-04,   2.83048459e-04,
           5.35565021e-04,   6.50145288e-04,   1.27295032e-03,
           1.84973318e-03,   4.02115844e-03,   8.64035450e-03,
           1.85660869e-02,   2.67588310e-02,   4.13986742e-02,
           3.75211351e-02,   3.99451032e-02,   2.19815839e-02,
           1.46964323e-02,   7.59503990e-03,   3.71277868e-03,
           2.43205787e-03,   1.36608526e-03,   1.10033224e-03,
           7.15283619e-04,   4.24796075e-04,   9.07893118e-05,
          -3.93107002e-05,   3.73856092e-05,   1.34933362e-04,
           8.41469518e-05],
        [  3.12657416e-04,   2.02556097e-04,   1.18483724e-04,
           1.16338248e-04,   2.58727494e-04,   3.82861763e-04,
           6.47869776e-04,   7.22532743e-04,   1.85770064e-03,
           2.38643191e-03,   5.77186933e-03,   1.30834254e-02,
           3.21811736e-02,   5.02903536e-02,   7.60164186e-02,
           8.23163912e-02,   6.55187890e-02,   5.04501201e-02,
           2.74458304e-02,   1.12441424e-02,   5.38016111e-03,
           3.53081268e-03,   1.61524036e-03,   1.24434254e-03,
           7.16076291e-04,   4.93183674e-04,   2.25529220e-04,
           8.95171906e-05,   7.44353310e-05,   9.23250846e-05,
           5.80607739e-05],
        [  3.23481683e-04,   1.60379801e-04,   6.94236369e-05,
           1.16357034e-04,   3.34398763e-04,   4.58175607e-04,
           7.41398602e-04,   8.16271116e-04,   2.42795702e-03,
           2.75838282e-03,   7.24995369e-03,   2.07453221e-02,
           4.65088971e-02,   7.89608508e-02,   1.19249240e-01,
           1.33808449e-01,   1.13654882e-01,   8.49788636e-02,
           4.27414551e-02,   2.24782787e-02,   6.53822487e-03,
           4.25153133e-03,   1.61689415e-03,   1.29075698e-03,
           6.18802325e-04,   4.10659268e-04,   2.27874494e-04,
           1.22942773e-04,   4.52497516e-05,   1.28073461e-05,
           1.60166855e-05],
        [  1.58397612e-04,   4.24731224e-05,   3.01601358e-05,
           1.21667144e-04,   3.14779114e-04,   4.15842806e-04,
           7.71404011e-04,   8.69169366e-04,   2.61580618e-03,
           2.77252262e-03,   7.77245732e-03,   2.10336000e-02,
           5.16618006e-02,   8.90334770e-02,   1.32652745e-01,
           1.57979146e-01,   1.25354961e-01,   1.05975635e-01,
           4.87471484e-02,   2.40709223e-02,   6.75325841e-03,
           4.41299053e-03,   1.51851110e-03,   1.27687806e-03,
           5.12708386e-04,   2.49664765e-04,   1.04957260e-04,
           6.98268996e-05,   1.52905995e-05,  -5.66302078e-05,
          -1.03157217e-04],
        [  2.37388813e-05,  -9.44329622e-06,   7.28973610e-05,
           1.45914135e-04,   2.14597429e-04,   2.61879875e-04,
           6.45024644e-04,   7.76131405e-04,   2.32133828e-03,
           2.64355657e-03,   7.23631820e-03,   1.91156585e-02,
           3.87151986e-02,   7.85027444e-02,   1.04771100e-01,
           1.32386833e-01,   1.13936760e-01,   8.50518569e-02,
           4.59707938e-02,   1.98367629e-02,   6.28389278e-03,
           4.06758115e-03,   1.48181198e-03,   1.24518829e-03,
           5.25051320e-04,   2.38819033e-04,   6.78987781e-05,
           5.05482931e-05,   4.51459491e-05,  -2.35327320e-06,
          -8.10507045e-05],
        [  5.06010019e-06,   3.72947142e-07,   1.22372061e-04,
           1.55802540e-04,   1.24519778e-04,   1.23926759e-04,
           4.53011831e-04,   6.43062172e-04,   1.81168271e-03,
           2.53919838e-03,   6.02304796e-03,   1.33035155e-02,
           3.03080510e-02,   4.51316275e-02,   6.80428073e-02,
           7.87105411e-02,   6.83350638e-02,   5.41879274e-02,
           3.03043667e-02,   1.30727030e-02,   5.33299660e-03,
           3.26684839e-03,   1.41233020e-03,   1.17955718e-03,
           6.65704778e-04,   4.21560457e-04,   1.78371061e-04,
           9.97985844e-05,   1.39496973e-04,   2.18161294e-04,
           2.08089550e-04],
        [  2.39601777e-05,   2.21635291e-05,   1.14513583e-04,
           1.19043048e-04,   7.76750167e-05,   1.02026192e-04,
           4.02551843e-04,   6.65154657e-04,   1.42666139e-03,
           2.27654330e-03,   4.46224120e-03,   8.39890260e-03,
           1.30807683e-02,   2.55265832e-02,   3.33221108e-02,
           3.61914933e-02,   3.77908386e-02,   2.31504384e-02,
           1.68734826e-02,   8.33798666e-03,   4.00891248e-03,
           2.29511573e-03,   1.17400289e-03,   9.92648304e-04,
           7.22233963e-04,   5.48228098e-04,   2.67189520e-04,
           1.37489595e-04,   2.13896274e-04,   3.94635805e-04,
           4.39600262e-04],
        [  6.23589949e-05,   8.92760654e-05,   9.79828910e-05,
           6.91207024e-05,   8.57825144e-05,   2.02807147e-04,
           5.06308512e-04,   7.66013341e-04,   1.16996281e-03,
           1.70234416e-03,   2.81825359e-03,   4.55184840e-03,
           6.95564412e-03,   9.65203252e-03,   1.53076285e-02,
           1.43200746e-02,   1.65908057e-02,   1.03640771e-02,
           7.39802886e-03,   4.66627069e-03,   2.64005945e-03,
           1.50676281e-03,   8.05904623e-04,   6.65010943e-04,
           5.85204456e-04,   4.94408188e-04,   2.34044710e-04,
           1.05874708e-04,   2.03739924e-04,   3.64033622e-04,
           3.66085063e-04],
        [  8.67988638e-05,   1.41086144e-04,   1.11501467e-04,
           1.00657169e-04,   1.85763391e-04,   3.11478914e-04,
           4.91909974e-04,   6.40244340e-04,   8.55031773e-04,
           1.11193419e-03,   1.64164591e-03,   2.48288224e-03,
           3.66801536e-03,   4.93864343e-03,   5.92702301e-03,
           6.39919844e-03,   6.21232018e-03,   5.31890662e-03,
           3.92624969e-03,   2.61738058e-03,   1.64109841e-03,
           1.00608368e-03,   5.06966840e-04,   3.94159695e-04,
           4.27827297e-04,   3.86294152e-04,   1.48756444e-04,
           4.65777812e-05,   1.58774055e-04,   2.52897677e-04,
           2.00576906e-04],
        [  3.37751662e-05,   7.41251060e-05,   1.17368654e-04,
           2.04129843e-04,   3.13290016e-04,   3.23200540e-04,
           2.81042245e-04,   2.99956737e-04,   4.86716832e-04,
           7.54022039e-04,   1.09415397e-03,   1.52351952e-03,
           2.05146684e-03,   2.48926901e-03,   2.69005261e-03,
           2.72322842e-03,   2.66885129e-03,   2.44992133e-03,
           1.99378422e-03,   1.47749775e-03,   1.01475196e-03,
           6.71428046e-04,   4.07902524e-04,   3.35486606e-04,
           3.32384778e-04,   2.58078508e-04,   8.59883512e-05,
           4.60952906e-05,   1.45739934e-04,   1.80587085e-04,
           1.13000431e-04],
        [ -4.24637801e-06,  -1.40777156e-05,   8.70996955e-05,
           2.22028204e-04,   2.97378720e-04,   2.60808825e-04,
           1.91989529e-04,   1.85491386e-04,   3.34152573e-04,
           5.63626410e-04,   8.59337917e-04,   1.17473991e-03,
           1.52407016e-03,   1.90831302e-03,   2.20574974e-03,
           2.22596643e-03,   1.98955438e-03,   1.70866319e-03,
           1.43378158e-03,   1.09570054e-03,   7.16906914e-04,
           5.33116632e-04,   5.02312207e-04,   4.49821993e-04,
           2.42741953e-04,   9.23747939e-05,   7.99639820e-05,
           1.36120565e-04,   1.63458535e-04,   1.22343510e-04,
           3.93863447e-05],
        [  5.69609983e-05,   1.73523185e-05,   6.40260187e-05,
           1.20879384e-04,   1.43646175e-04,   2.16061788e-04,
           3.54024611e-04,   4.07661224e-04,   3.58933845e-04,
           3.38174548e-04,   4.71632462e-04,   6.41264895e-04,
           7.74718646e-04,   9.90612898e-04,   1.26875984e-03,
           1.31991541e-03,   1.06533815e-03,   8.48297495e-04,
           8.23326525e-04,   7.31756329e-04,   4.78470349e-04,
           4.04468388e-04,   5.31740312e-04,   4.86872537e-04,
           1.58969007e-04,   6.51523487e-06,   1.51630229e-04,
           2.65381881e-04,   1.91450192e-04,   7.43656929e-05,
          -1.58241182e-05],
        [  1.17583440e-04,   7.81745766e-05,   6.28807975e-05,
           5.76429666e-05,   8.64068497e-05,   2.51753721e-04,
           4.96915658e-04,   5.61750669e-04,   3.79482983e-04,
           1.83867858e-04,   1.96382141e-04,   3.05051944e-04,
           3.97037104e-04,   5.88906172e-04,   8.83791305e-04,
           9.97585594e-04,   7.86893943e-04,   5.53074176e-04,
           5.21226437e-04,   4.94938286e-04,   3.33534525e-04,
           2.57302396e-04,   3.20619816e-04,   2.87285424e-04,
           8.89307703e-05,   4.46845261e-05,   2.19037989e-04,
           3.19396699e-04,   2.09293052e-04,   8.27520344e-05,
           3.74381343e-05],
        [  8.29820856e-05,   6.90407614e-05,   7.88107500e-05,
           1.12841059e-04,   1.79206720e-04,   2.88871757e-04,
           3.90221598e-04,   3.79044010e-04,   2.56998552e-04,
           1.32744783e-04,   9.68408785e-05,   1.29585766e-04,
           2.00170762e-04,   3.32332536e-04,   4.80854040e-04,
           5.01023198e-04,   3.60588630e-04,   2.58127344e-04,
           3.06765636e-04,   3.57608777e-04,   2.83811503e-04,
           1.61550328e-04,   6.34512617e-05,   4.94756819e-07,
          -1.60358031e-05,   4.90231869e-05,   1.74128450e-04,
           2.60341651e-04,   2.19396825e-04,   1.49694664e-04,
           1.41322962e-04],
        [  1.49950465e-05,   7.63047574e-05,   1.41467142e-04,
           1.94792388e-04,   2.02353724e-04,   1.83372205e-04,
           1.58346622e-04,   1.04015133e-04,   7.28349551e-05,
           7.21577380e-05,   5.38176573e-05,   7.33939814e-05,
           1.82003147e-04,   2.95805483e-04,   2.92530225e-04,
           1.55183196e-04,   2.42199985e-05,   6.74540934e-05,
           2.36359934e-04,   3.27022077e-04,   2.74324237e-04,
           1.71450374e-04,   4.91856554e-05,  -4.54968213e-05,
          -4.42001874e-05,   1.27510330e-05,   7.70893821e-05,
           1.61725053e-04,   2.11200968e-04,   1.84710239e-04,
           1.30427681e-04],
        [  3.82523067e-05,   1.62776094e-04,   2.29700861e-04,
           2.13740321e-04,   8.46937110e-05,   2.72549073e-06,
           3.94645940e-05,   1.42122672e-05,  -6.94230403e-05,
          -9.39216843e-05,  -6.60941514e-05,   3.26379231e-05,
           2.00269962e-04,   2.96857586e-04,   2.09510399e-04,
           1.20158993e-05,  -1.01647180e-04,  -1.86980396e-05,
           1.49961386e-04,   2.03447606e-04,   1.50521664e-04,
           1.42593373e-04,   1.54100853e-04,   1.37039591e-04,
           1.27139734e-04,   1.02055506e-04,   4.87680190e-05,
           6.11553332e-05,   1.30138273e-04,   9.36011056e-05,
          -3.46084780e-05],
        [  1.45188038e-04,   2.12477811e-04,   2.30793623e-04,
           1.70628540e-04,   3.61743005e-05,   2.19909125e-05,
           1.38719668e-04,   1.01310732e-04,  -9.27221336e-05,
          -1.97339192e-04,  -1.34774629e-04,   3.83440201e-05,
           2.23026072e-04,   2.67371885e-04,   1.41579862e-04,
           4.79614926e-07,  -3.39062790e-05,   1.89720467e-05,
           8.16576139e-05,   7.70369734e-05,   4.06700310e-05,
           7.07829677e-05,   1.45996106e-04,   2.27258963e-04,
           2.88130395e-04,   2.59734894e-04,   1.29991880e-04,
           3.11676340e-05,   1.03388302e-05,  -7.12372421e-05,
          -1.81315132e-04],
        [  1.74572386e-04,   1.21468227e-04,   1.03557053e-04,
           9.54290008e-05,   1.19859134e-04,   2.27726166e-04,
           3.10518721e-04,   2.02525058e-04,  -9.97376901e-06,
          -1.25466555e-04,  -8.44838651e-05,   5.78080217e-05,
           1.99554896e-04,   1.96891691e-04,   6.75892225e-05,
           6.29681153e-07,   4.65245939e-05,   1.00516285e-04,
           1.16317417e-04,   9.37591758e-05,   5.00349415e-05,
           3.41828090e-05,   6.40527724e-05,   1.28187123e-04,
           2.03160991e-04,   2.64978939e-04,   2.48178956e-04,
           1.35484006e-04,  -2.65322942e-05,  -1.60183437e-04,
          -1.74708708e-04],
        [  7.79956536e-05,   4.45946062e-05,   6.88394939e-05,
           1.20791679e-04,   2.06197074e-04,   2.92086275e-04,
           2.88313808e-04,   1.82209784e-04,   7.11210378e-05,
           9.09970095e-06,  -1.04216479e-05,   2.26099164e-05,
           9.15258279e-05,   9.49135938e-05,   2.28152076e-05,
          -9.77456921e-07,   6.45484615e-05,   1.56735507e-04,
           2.20143396e-04,   2.00825802e-04,   1.07963591e-04,
           4.76653659e-05,   5.46586489e-05,   7.16102004e-05,
           9.44839630e-05,   1.83187702e-04,   2.52839702e-04,
           1.77972455e-04,  -1.64299454e-05,  -1.33052861e-04,
          -7.75216977e-05]]], dtype=np.float32)
			tbhdu=fits.BinTableHDU.from_columns(fits.ColDefs([fits.Column(name='PSF_MASK' , format='961E' , dim='(31, 31, 1)' , array=arr)]))
			hdr=tbhdu.header
			hdr.set('EXTNAME' , 'PSF_DATA', 'TABLE NAME')
			hdr.set('LOADED' , 36 , 'Number of loaded sources')
			hdr.set('ACCEPTED' , 32 , 'Number of accepted sources')
			hdr.set('CHI2' , 1.4190 , 'Final Chi2')
			hdr.set('POLNAXIS' , 0 , 'Number of context parameters')
			hdr.set('POLNGRP' , 0 , 'Number of context groups')
			hdr.set('PSF_FWHM' , 2.5813 , 'PSF FWHM')
			hdr.set('PSF_SAMP' , 0.5000 , 'Sampling step of the PSF data')
			hdr.set('PSFNAXIS' , 3 , 'Dimensionality of the PSF data')
			hdr.set('PSFAXIS1' , 31 , 'Number of element along this axis')
			hdr.set('PSFAXIS2' , 31 , 'Number of element along this axis')
			hdr.set('PSFAXIS3' , 1 , 'Number of element along this axis')             
			
			thdulist=fits.HDUList()
			thdulist.append(tbhdu)
			thdulist.writeto(self._get_psf_filepath())
			logger.debug("Wrote %s" % (self._get_psf_filepath()))
		else:
			logger.debug("Default psf file already exists, I don't overwrite it.")
			
	def _clean_workdir(self):
		"""
		Removes the config/param files related to this instance, to allow for a fresh restart.
		Files related to specific images are not removed.
		"""
		toremove = [self._get_config_filepath(), self._get_params_filepath(), self._get_conv_filepath(), self._get_psf_filepath()]
		for filepath in toremove:
			if os.path.exists(filepath):	
				logger.debug("Removing existing file %s..." % (filepath))
				os.remove(filepath)


	def _write_assoc(self, cat, xname, yname, imgname):
		"""
		Writes a plain text file which can be used as sextractor input for the ASSOC identification.
		And "index" for each source is generated, it gets used to identify galaxies.
		"""
		
		#if assoc_xname not in assoc_cat.colnames or assoc_yname not in assoc_cat.colnames:
		#	raise RuntimeError("I don't have columns %s or %s" % (assoc_xname, assoc_yname))
		
		if os.path.exists(self._get_assoc_filepath(imgname)):	
			logger.warning("ASSOC file already exists, I will overwrite it")

		lines = []
		for (number, row) in enumerate(cat):
			# Seems safe(r) to not use row.index but our own number.
			lines.append("%.3f\t%.3f\t%i\n" % (row[xname], row[yname], number))

		lines = "".join(lines)
		f = open(self._get_assoc_filepath(imgname), "w")
		f.writelines(lines)
		f.close()
		logger.debug("Wrote ASSOC file %s..." % (self._get_assoc_filepath(imgname)))

	
	def _add_prefix(self, table, prefix):
		"""
		Modifies the column names of a table by prepending the prefix *in place*.
		Skips the VECTOR_ASSOC stuff !
		"""
		if prefix == "":
			return
		for colname in table.colnames:
			if colname not in ["VECTOR_ASSOC", "VECTOR_ASSOC_1", "VECTOR_ASSOC_2"]:
				table.rename_column(colname, prefix + colname)

		

	def __call__(self, imgfilepath, imgname=None, assoc_cat=None, assoc_xname="x", assoc_yname="y",
                returncat=True, prefix="", writelog=True):
		"""
		Runs SExtractor on a given image.
		
		:param imgfilepath: Path to the input FITS image I should run on
		:param assoc_cat:  optional input catalog (astropy table), if you want to use the ASSOC helper
		:param assoc_xname: x coordinate name I should use in the ASSOC helper
		:param assoc_yname: idem
		
		:param returncat: by default I read the SExtractor output catalog and return it as an astropy
                        table.
			If set to False, I do not attempt to read it.
		:param prefix: will be prepended to the column names of the astropy table that I return
		:type prefix: string
		:param writelog: if True I save the sextractor command line input and output into a dedicated
                        log file in the workdir.
		
		:returns: a dict containing the keys:
		          
			* **catfilepath**: the path to the sextractor output catalog file
			* **table**: the astropy table of the output catalog (if returncat was not set to False)
			* **workdir**: the path to the workdir (all my internal files are there)
			* **logfilepath**: the path to the SExtractor log file (in the workdir)
		
		Everything related to this particular image stays within this method, the SExtractor instance
		(in particular config) is not modified !
		"""

		starttime = datetime.now()
		
		# Let's first check if the image file exists.
		if not os.path.exists(imgfilepath):
			raise IOError("The image file %s does not exist." % imgfilepath)
		logger.info("Preparing to run SExtractor on %s..." % imgfilepath)

		if imgname == None:
			imgname = os.path.splitext(os.path.basename(imgfilepath))[0]
		logger.debug("Using imgname '%s'..." % (imgname))		
		
		# We make a deep copy of the config, that we can modify with settings related to this particular
		# image.
		imgconfig = copy.deepcopy(self.config)
		
		# We set the catalog name :
		imgconfig["CATALOG_NAME"] = self._get_cat_filepath(imgname)
		if os.path.exists(self._get_cat_filepath(imgname)):
			logger.warning("Output catalog %s already exists, I will overwrite it" % (self._get_cat_filepath(imgname)))
		
		
		# We prepare the ASSOC catalog file, if needed
		if assoc_cat is not None:
			
			logger.info("I will run in ASSOC mode, trying to find %i sources..." % (len(assoc_cat)))
			if "VECTOR_ASSOC(3)" not in self.params:
				raise RuntimeError("To use the ASSOC helper, you have to add 'VECTOR_ASSOC(3)' to the params")
			if assoc_xname not in assoc_cat.colnames or assoc_yname not in assoc_cat.colnames:
				raise RuntimeError("I don't have columns %s or %s" % (assoc_xname, assoc_yname))
			if "VECTOR_ASSOC_2" in assoc_cat.colnames:
				raise RuntimeError("Do not give me an assoc_cat that already contains a column VECTOR_ASSOC_2")
			for param in self.params + [prefix + "assoc_flag"]:
				# This is not 100% correct, as some params might be vectors.
				if prefix + param in assoc_cat.colnames:
					raise RuntimeError("Your assoc_cat already has a column named %s, fix this" % (prefix + param))
			
			self._write_assoc(cat=assoc_cat, xname=assoc_xname, yname=assoc_yname, imgname=imgname)
		
			imgconfig["ASSOC_DATA"] = "1, 2, 3"
			imgconfig["ASSOC_NAME"] = self._get_assoc_filepath(imgname)
			imgconfig["ASSOC_PARAMS"] = "1, 2"
			if "ASSOC_RADIUS" not in imgconfig:
				logger.warning("ASSOC_RADIUS not specified, using a default of 10.0")
				imgconfig["ASSOC_RADIUS"] = 10.0
			if "ASSOC_TYPE" not in imgconfig:
				logger.warning("ASSOC_TYPE not specified, using a default NEAREST")
				imgconfig["ASSOC_TYPE"] = "NEAREST"
			if "ASSOCSELEC_TYPE" in imgconfig:
				raise RuntimeError("Sorry, you cannot mess with ASSOCSELEC_TYPE yourself when using the helper. I'm using MATCHED.")
			imgconfig["ASSOCSELEC_TYPE"] = "MATCHED"

		
		# We write the input files (if needed)
		self._write_default_config()
		self._write_params()
		self._write_default_conv()
		self._write_default_psf()
		
		# We build the command line arguments
		popencmd = [self.sexpath, imgfilepath, "-c", self._get_config_filepath()]
		if self.nice != None: # We prepend the nice command
			popencmd[:0] = ["nice", "-n", str(self.nice)]
		
		# We add the current state of config
		for (key, value) in imgconfig.items():
			popencmd.append("-"+str(key))
			popencmd.append(str(value).replace(' ',''))
		
		# And we run
		logger.info("Starting SExtractor now, with niceness %s..." % (self.nice))
		logger.debug("Running with command %s..." % (popencmd))
		p = subprocess.Popen(popencmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate()
		
		if writelog:
			logfile = open(self._get_log_filepath(imgname), "w")
			logfile.write("SExtractor was called with :\n")
			logfile.write(" ".join(popencmd))
			logfile.write("\n\nA nicer view of the config:\n")
			logfile.write("\n".join(["%30s : %30s" % (str(key), str(value)) for (key, value) in imgconfig.items()]))
			logfile.write("\n\n####### stdout #######\n")
			logfile.write(out.decode(encoding='UTF-8'))
			logfile.write("\n####### stderr #######\n")
			logfile.write(err.decode(encoding='UTF-8'))
			logfile.write("\n")
			logfile.close()
			
		logger.info("SExtractor stderr:")
		logger.info(err)
		
		if not "All done" in err.decode(encoding='UTF-8'):
			logger.warning("Ouch, something seems wrong, check SExtractor log: %s" % self._get_log_filepath(imgname))
		
		endtime = datetime.now()
		logger.info("Running SExtractor done, it took %.2f seconds." % \
                                ((endtime - starttime).total_seconds()))

		# Let's check if this worked.
		if not os.path.isfile(self._get_cat_filepath(imgname)):
			raise RuntimeError("It seems that SExtractor did not write the file '%s'. Check SExtractor log: %s" % (self._get_cat_filepath(imgname), self._get_log_filepath(imgname)))

		# We return a dict. It always contains at least the path to the sextractor catalog:
		output = {"catfilepath":self._get_cat_filepath(imgname), "workdir":self.workdir}
		if writelog:
			output["logfilepath"] = self._get_log_filepath(imgname)
		
		# And we read the output, if asked for:
		if returncat:
			if assoc_cat is None:
				sextable = astropy.table.Table.read(self._get_cat_filepath(imgname),
                                                                    format="ascii.sextractor")
				logger.info("Read %i objects from the SExtractor output catalog" % (len(sextable)))
				self._add_prefix(sextable, prefix)
				output["table"] = sextable
				
			else: # We have to process the output catalog, merging it.
			
				# We add the "number" column to the assoc_cat, calling it VECTOR_ASSOC_2:
				intable = copy.deepcopy(assoc_cat)
				intable["VECTOR_ASSOC_2"] = range(len(assoc_cat))
								
				# We read in the SExtractor output:					
				sextable = astropy.table.Table.read(self._get_cat_filepath(imgname),
                                                                    format="ascii.sextractor")
				logger.info("Read %i objects from the SExtractor output catalog" % (len(sextable)))
				self._add_prefix(sextable, prefix)
				sextable.remove_columns(["VECTOR_ASSOC", "VECTOR_ASSOC_1"])
				
				# Due to what seems to be a bug in SExtractor (version 2.19.5 and earlier),
				# we need to kick out "duplicated" (same VECTOR_ASSOC_2) rows.
				# That's weird, as in principle we asked to keep the NEAREST !
				sortedassoc = np.sort(sextable["VECTOR_ASSOC_2"].data)
				duplassoc = list(np.unique(sortedassoc[sortedassoc[1:] == sortedassoc[:-1]]))
				# The unique is here as there might be more than 2 identical numbers...
				if len(duplassoc) > 0:
					logger.warning("%i sources from the SExtractor catalog are strange duplicates (bug ?), I discard them." % (len(duplassoc)))
					rowindices_to_remove = []
					for row in sextable:
						if row["VECTOR_ASSOC_2"] in duplassoc:
							rowindices_to_remove.append(row.index)
					sextable.remove_rows(rowindices_to_remove)
					
				if len(sextable) == 0:
					raise RuntimeError("SExtractor has returned no ASSOC match")
							
				# We merge the tables, keeping all entries of the "intable"
				joined = astropy.table.join(intable, sextable,
					join_type='left', keys='VECTOR_ASSOC_2',
					 # raises an error in case of metadata conflict.
					metadata_conflicts = "error",
					 # Will only be used in case of column name conflicts.
					table_names = ['ASSOC', 'SEx'],
					uniq_col_name = "{table_name}_{col_name}"
					)
				
				# This join does not mix the order, as the output is sorted according to our own
				# VECTOR_ASSOC_2
				
				# We remove the last ASSOC column:
				joined.remove_columns(["VECTOR_ASSOC_2"])
				#assert len(intable) == len(joined)
				# More explicit:
				if not len(intable) == len(joined):
					raise RuntimeError("Problem with joined tables: intable has %i rows, joined has %i. %s %s" % (len(intable), len(joined), intable.colnames, joined.colnames))
				
				# The join might return a **masked** table.
				# In any case, we add one simply-named column with a flag telling if the
				# identification has worked.
				
				if joined.masked:
					logger.info("ASSOC join done, my output is a masked table.")
					joined[prefix + "assoc_flag"] = joined[joined.colnames[-1]].mask == False
					nfound = sum(joined[prefix + "assoc_flag"])
					logger.info("I could find %i out of %i sources (%i are missing)" % \
                                                        (nfound, len(assoc_cat), len(assoc_cat)-nfound))
				
				else:
					logger.info("ASSOC join done, I could find all your sources, my output is not masked.")
					joined[prefix + "assoc_flag"] = [True] * len(joined)
					
				
				output["table"] = joined
		
		return output
		
		
#	def destroy(self):
#		"""
#		Removes the complete working dir, careful with this.
#		"""
#		# No, this is way to dangerous, workdir could be "."
#		#shutil.rmtree(self.workdir)


	# Some class attributes:
	# We give this fullparamtxt here as some earlier versions of sextractor are not able to spit it out.
	# It's only used to check your params for typos, anyway.
	
	fullparamtxt = """
#NUMBER                 Running object number                                     
#EXT_NUMBER             FITS extension number                                     
#FLUX_ISO               Isophotal flux                                             [count]
#FLUXERR_ISO            RMS error for isophotal flux                               [count]
#MAG_ISO                Isophotal magnitude                                        [mag]
#MAGERR_ISO             RMS error for isophotal magnitude                          [mag]
#FLUX_ISOCOR            Corrected isophotal flux                                   [count]
#FLUXERR_ISOCOR         RMS error for corrected isophotal flux                     [count]
#MAG_ISOCOR             Corrected isophotal magnitude                              [mag]
#MAGERR_ISOCOR          RMS error for corrected isophotal magnitude                [mag]
#FLUX_APER              Flux vector within fixed circular aperture(s)              [count]
#FLUXERR_APER           RMS error vector for aperture flux(es)                     [count]
#MAG_APER               Fixed aperture magnitude vector                            [mag]
#MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
#FLUX_AUTO              Flux within a Kron-like elliptical aperture                [count]
#FLUXERR_AUTO           RMS error for AUTO flux                                    [count]
#MAG_AUTO               Kron-like elliptical aperture magnitude                    [mag]
#MAGERR_AUTO            RMS error for AUTO magnitude                               [mag]
#FLUX_PETRO             Flux within a Petrosian-like elliptical aperture           [count]
#FLUXERR_PETRO          RMS error for PETROsian flux                               [count]
#MAG_PETRO              Petrosian-like elliptical aperture magnitude               [mag]
#MAGERR_PETRO           RMS error for PETROsian magnitude                          [mag]
#FLUX_BEST              Best of FLUX_AUTO and FLUX_ISOCOR                          [count]
#FLUXERR_BEST           RMS error for BEST flux                                    [count]
#MAG_BEST               Best of MAG_AUTO and MAG_ISOCOR                            [mag]
#MAGERR_BEST            RMS error for MAG_BEST                                     [mag]
#FLUX_WIN               Gaussian-weighted flux                                     [count]
#FLUXERR_WIN            RMS error for WIN flux                                     [count]
#MAG_WIN                Gaussian-weighted magnitude                                [mag]
#MAGERR_WIN             RMS error for MAG_WIN                                      [mag]
#FLUX_SOMFIT            Flux derived from SOM fit                                  [count]
#FLUXERR_SOMFIT         RMS error for SOMFIT flux                                  [count]
#MAG_SOMFIT             Magnitude derived from SOM fit                             [mag]
#MAGERR_SOMFIT          Magnitude error derived from SOM fit                       [mag]
#ERROR_SOMFIT           Reduced Chi-square error of the SOM fit                   
#VECTOR_SOMFIT          Position vector of the winning SOM node                   
#KRON_RADIUS            Kron apertures in units of A or B                         
#PETRO_RADIUS           Petrosian apertures in units of A or B                    
#BACKGROUND             Background at centroid position                            [count]
#THRESHOLD              Detection threshold above background                       [count]
#FLUX_MAX               Peak flux above background                                 [count]
#ISOAREA_IMAGE          Isophotal area above Analysis threshold                    [pixel**2]
#ISOAREAF_IMAGE         Isophotal area (filtered) above Detection threshold        [pixel**2]
#XMIN_IMAGE             Minimum x-coordinate among detected pixels                 [pixel]
#YMIN_IMAGE             Minimum y-coordinate among detected pixels                 [pixel]
#XMAX_IMAGE             Maximum x-coordinate among detected pixels                 [pixel]
#YMAX_IMAGE             Maximum y-coordinate among detected pixels                 [pixel]
#XPEAK_IMAGE            x-coordinate of the brightest pixel                        [pixel]
#YPEAK_IMAGE            y-coordinate of the brightest pixel                        [pixel]
#XPEAK_WORLD            World-x coordinate of the brightest pixel                  [deg]
#YPEAK_WORLD            World-y coordinate of the brightest pixel                  [deg]
#ALPHAPEAK_SKY          Right ascension of brightest pix (native)                  [deg]
#DELTAPEAK_SKY          Declination of brightest pix (native)                      [deg]
#ALPHAPEAK_J2000        Right ascension of brightest pix (J2000)                   [deg]
#DELTAPEAK_J2000        Declination of brightest pix (J2000)                       [deg]
#ALPHAPEAK_B1950        Right ascension of brightest pix (B1950)                   [deg]
#DELTAPEAK_B1950        Declination of brightest pix (B1950)                       [deg]
#X_IMAGE                Object position along x                                    [pixel]
#Y_IMAGE                Object position along y                                    [pixel]
#X_IMAGE_DBL            Object position along x (double precision)                 [pixel]
#Y_IMAGE_DBL            Object position along y (double precision)                 [pixel]
#X_WORLD                Barycenter position along world x axis                     [deg]
#Y_WORLD                Barycenter position along world y axis                     [deg]
#X_MAMA                 Barycenter position along MAMA x axis                      [m**(-6)]
#Y_MAMA                 Barycenter position along MAMA y axis                      [m**(-6)]
#ALPHA_SKY              Right ascension of barycenter (native)                     [deg]
#DELTA_SKY              Declination of barycenter (native)                         [deg]
#ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
#DELTA_J2000            Declination of barycenter (J2000)                          [deg]
#ALPHA_B1950            Right ascension of barycenter (B1950)                      [deg]
#DELTA_B1950            Declination of barycenter (B1950)                          [deg]
#X2_IMAGE               Variance along x                                           [pixel**2]
#Y2_IMAGE               Variance along y                                           [pixel**2]
#XY_IMAGE               Covariance between x and y                                 [pixel**2]
#X2_WORLD               Variance along X-WORLD (alpha)                             [deg**2]
#Y2_WORLD               Variance along Y-WORLD (delta)                             [deg**2]
#XY_WORLD               Covariance between X-WORLD and Y-WORLD                     [deg**2]
#CXX_IMAGE              Cxx object ellipse parameter                               [pixel**(-2)]
#CYY_IMAGE              Cyy object ellipse parameter                               [pixel**(-2)]
#CXY_IMAGE              Cxy object ellipse parameter                               [pixel**(-2)]
#CXX_WORLD              Cxx object ellipse parameter (WORLD units)                 [deg**(-2)]
#CYY_WORLD              Cyy object ellipse parameter (WORLD units)                 [deg**(-2)]
#CXY_WORLD              Cxy object ellipse parameter (WORLD units)                 [deg**(-2)]
#A_IMAGE                Profile RMS along major axis                               [pixel]
#B_IMAGE                Profile RMS along minor axis                               [pixel]
#THETA_IMAGE            Position angle (CCW/x)                                     [deg]
#A_WORLD                Profile RMS along major axis (world units)                 [deg]
#B_WORLD                Profile RMS along minor axis (world units)                 [deg]
#THETA_WORLD            Position angle (CCW/world-x)                               [deg]
#THETA_SKY              Position angle (east of north) (native)                    [deg]
#THETA_J2000            Position angle (east of north) (J2000)                     [deg]
#THETA_B1950            Position angle (east of north) (B1950)                     [deg]
#ERRX2_IMAGE            Variance of position along x                               [pixel**2]
#ERRY2_IMAGE            Variance of position along y                               [pixel**2]
#ERRXY_IMAGE            Covariance of position between x and y                     [pixel**2]
#ERRX2_WORLD            Variance of position along X-WORLD (alpha)                 [deg**2]
#ERRY2_WORLD            Variance of position along Y-WORLD (delta)                 [deg**2]
#ERRXY_WORLD            Covariance of position X-WORLD/Y-WORLD                     [deg**2]
#ERRCXX_IMAGE           Cxx error ellipse parameter                                [pixel**(-2)]
#ERRCYY_IMAGE           Cyy error ellipse parameter                                [pixel**(-2)]
#ERRCXY_IMAGE           Cxy error ellipse parameter                                [pixel**(-2)]
#ERRCXX_WORLD           Cxx error ellipse parameter (WORLD units)                  [deg**(-2)]
#ERRCYY_WORLD           Cyy error ellipse parameter (WORLD units)                  [deg**(-2)]
#ERRCXY_WORLD           Cxy error ellipse parameter (WORLD units)                  [deg**(-2)]
#ERRA_IMAGE             RMS position error along major axis                        [pixel]
#ERRB_IMAGE             RMS position error along minor axis                        [pixel]
#ERRTHETA_IMAGE         Error ellipse position angle (CCW/x)                       [deg]
#ERRA_WORLD             World RMS position error along major axis                  [deg]
#ERRB_WORLD             World RMS position error along minor axis                  [deg]
#ERRTHETA_WORLD         Error ellipse pos. angle (CCW/world-x)                     [deg]
#ERRTHETA_SKY           Native error ellipse pos. angle (east of north)            [deg]
#ERRTHETA_J2000         J2000 error ellipse pos. angle (east of north)             [deg]
#ERRTHETA_B1950         B1950 error ellipse pos. angle (east of north)             [deg]
#XWIN_IMAGE             Windowed position estimate along x                         [pixel]
#YWIN_IMAGE             Windowed position estimate along y                         [pixel]
#XWIN_WORLD             Windowed position along world x axis                       [deg]
#YWIN_WORLD             Windowed position along world y axis                       [deg]
#ALPHAWIN_SKY           Windowed right ascension  (native)                         [deg]
#DELTAWIN_SKY           Windowed declination (native)                              [deg]
#ALPHAWIN_J2000         Windowed right ascension (J2000)                           [deg]
#DELTAWIN_J2000         windowed declination (J2000)                               [deg]
#ALPHAWIN_B1950         Windowed right ascension (B1950)                           [deg]
#DELTAWIN_B1950         Windowed declination (B1950)                               [deg]
#X2WIN_IMAGE            Windowed variance along x                                  [pixel**2]
#Y2WIN_IMAGE            Windowed variance along y                                  [pixel**2]
#XYWIN_IMAGE            Windowed covariance between x and y                        [pixel**2]
#X2WIN_WORLD            Windowed variance along X-WORLD (alpha)                    [deg**2]
#Y2WIN_WORLD            Windowed variance along Y-WORLD (delta)                    [deg**2]
#XYWIN_WORLD            Windowed covariance between X-WORLD and Y-WORLD            [deg**2]
#CXXWIN_IMAGE           Windowed Cxx object ellipse parameter                      [pixel**(-2)]
#CYYWIN_IMAGE           Windowed Cyy object ellipse parameter                      [pixel**(-2)]
#CXYWIN_IMAGE           Windowed Cxy object ellipse parameter                      [pixel**(-2)]
#CXXWIN_WORLD           Windowed Cxx object ellipse parameter (WORLD units)        [deg**(-2)]
#CYYWIN_WORLD           Windowed Cyy object ellipse parameter (WORLD units)        [deg**(-2)]
#CXYWIN_WORLD           Windowed Cxy object ellipse parameter (WORLD units)        [deg**(-2)]
#AWIN_IMAGE             Windowed profile RMS along major axis                      [pixel]
#BWIN_IMAGE             Windowed profile RMS along minor axis                      [pixel]
#THETAWIN_IMAGE         Windowed position angle (CCW/x)                            [deg]
#AWIN_WORLD             Windowed profile RMS along major axis (world units)        [deg]
#BWIN_WORLD             Windowed profile RMS along minor axis (world units)        [deg]
#THETAWIN_WORLD         Windowed position angle (CCW/world-x)                      [deg]
#THETAWIN_SKY           Windowed position angle (east of north) (native)           [deg]
#THETAWIN_J2000         Windowed position angle (east of north) (J2000)            [deg]
#THETAWIN_B1950         Windowed position angle (east of north) (B1950)            [deg]
#ERRX2WIN_IMAGE         Variance of windowed pos along x                           [pixel**2]
#ERRY2WIN_IMAGE         Variance of windowed pos along y                           [pixel**2]
#ERRXYWIN_IMAGE         Covariance of windowed pos between x and y                 [pixel**2]
#ERRX2WIN_WORLD         Variance of windowed pos along X-WORLD (alpha)             [deg**2]
#ERRY2WIN_WORLD         Variance of windowed pos along Y-WORLD (delta)             [deg**2]
#ERRXYWIN_WORLD         Covariance of windowed pos X-WORLD/Y-WORLD                 [deg**2]
#ERRCXXWIN_IMAGE        Cxx windowed error ellipse parameter                       [pixel**(-2)]
#ERRCYYWIN_IMAGE        Cyy windowed error ellipse parameter                       [pixel**(-2)]
#ERRCXYWIN_IMAGE        Cxy windowed error ellipse parameter                       [pixel**(-2)]
#ERRCXXWIN_WORLD        Cxx windowed error ellipse parameter (WORLD units)         [deg**(-2)]
#ERRCYYWIN_WORLD        Cyy windowed error ellipse parameter (WORLD units)         [deg**(-2)]
#ERRCXYWIN_WORLD        Cxy windowed error ellipse parameter (WORLD units)         [deg**(-2)]
#ERRAWIN_IMAGE          RMS windowed pos error along major axis                    [pixel]
#ERRBWIN_IMAGE          RMS windowed pos error along minor axis                    [pixel]
#ERRTHETAWIN_IMAGE      Windowed error ellipse pos angle (CCW/x)                   [deg]
#ERRAWIN_WORLD          World RMS windowed pos error along major axis              [deg]
#ERRBWIN_WORLD          World RMS windowed pos error along minor axis              [deg]
#ERRTHETAWIN_WORLD      Windowed error ellipse pos. angle (CCW/world-x)            [deg]
#ERRTHETAWIN_SKY        Native windowed error ellipse pos. angle (east of north)   [deg]
#ERRTHETAWIN_J2000      J2000 windowed error ellipse pos. angle (east of north)    [deg]
#ERRTHETAWIN_B1950      B1950 windowed error ellipse pos. angle (east of north)    [deg]
#NITER_WIN              Number of iterations for WIN centering                    
#MU_THRESHOLD           Detection threshold above background                       [mag * arcsec**(-2)]
#MU_MAX                 Peak surface brightness above background                   [mag * arcsec**(-2)]
#ISOAREA_WORLD          Isophotal area above Analysis threshold                    [deg**2]
#ISOAREAF_WORLD         Isophotal area (filtered) above Detection threshold        [deg**2]
#ISO0                   Isophotal area at level 0                                  [pixel**2]
#ISO1                   Isophotal area at level 1                                  [pixel**2]
#ISO2                   Isophotal area at level 2                                  [pixel**2]
#ISO3                   Isophotal area at level 3                                  [pixel**2]
#ISO4                   Isophotal area at level 4                                  [pixel**2]
#ISO5                   Isophotal area at level 5                                  [pixel**2]
#ISO6                   Isophotal area at level 6                                  [pixel**2]
#ISO7                   Isophotal area at level 7                                  [pixel**2]
#FLAGS                  Extraction flags                                          
#FLAGS_WEIGHT           Weighted extraction flags                                 
#FLAGS_WIN              Flags for WINdowed parameters                             
#IMAFLAGS_ISO           FLAG-image flags OR'ed over the iso. profile              
#NIMAFLAGS_ISO          Number of flagged pixels entering IMAFLAGS_ISO            
#FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
#FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
#ELONGATION             A_IMAGE/B_IMAGE                                           
#ELLIPTICITY            1 - B_IMAGE/A_IMAGE                                       
#POLAR_IMAGE            (A_IMAGE^2 - B_IMAGE^2)/(A_IMAGE^2 + B_IMAGE^2)           
#POLAR_WORLD            (A_WORLD^2 - B_WORLD^2)/(A_WORLD^2 + B_WORLD^2)           
#POLARWIN_IMAGE         (AWIN^2 - BWIN^2)/(AWIN^2 + BWIN^2)                       
#POLARWIN_WORLD         (AWIN^2 - BWIN^2)/(AWIN^2 + BWIN^2)                       
#CLASS_STAR             S/G classifier output                                     
#VIGNET                 Pixel data around detection                                [count]
#VIGNET_SHIFT           Pixel data around detection, corrected for shift           [count]
#VECTOR_ASSOC           ASSOCiated parameter vector                               
#NUMBER_ASSOC           Number of ASSOCiated IDs                                  
#THRESHOLDMAX           Maximum threshold possible for detection                   [count]
#FLUX_GROWTH            Cumulated growth-curve                                     [count]
#FLUX_GROWTHSTEP        Step for growth-curves                                     [pixel]
#MAG_GROWTH             Cumulated magnitude growth-curve                           [mag]
#MAG_GROWTHSTEP         Step for growth-curves                                     [pixel]
#FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
#XPSF_IMAGE             X coordinate from PSF-fitting                              [pixel]
#YPSF_IMAGE             Y coordinate from PSF-fitting                              [pixel]
#XPSF_WORLD             PSF position along world x axis                            [deg]
#YPSF_WORLD             PSF position along world y axis                            [deg]
#ALPHAPSF_SKY           Right ascension of the fitted PSF (native)                 [deg]
#DELTAPSF_SKY           Declination of the fitted PSF (native)                     [deg]
#ALPHAPSF_J2000         Right ascension of the fitted PSF (J2000)                  [deg]
#DELTAPSF_J2000         Declination of the fitted PSF (J2000)                      [deg]
#ALPHAPSF_B1950         Right ascension of the fitted PSF (B1950)                  [deg]
#DELTAPSF_B1950         Declination of the fitted PSF (B1950)                      [deg]
#FLUX_PSF               Flux from PSF-fitting                                      [count]
#FLUXERR_PSF            RMS flux error for PSF-fitting                             [count]
#MAG_PSF                Magnitude from PSF-fitting                                 [mag]
#MAGERR_PSF             RMS magnitude error from PSF-fitting                       [mag]
#NITER_PSF              Number of iterations for PSF-fitting                      
#CHI2_PSF               Reduced chi2 from PSF-fitting                             
#ERRX2PSF_IMAGE         Variance of PSF position along x                           [pixel**2]
#ERRY2PSF_IMAGE         Variance of PSF position along y                           [pixel**2]
#ERRXYPSF_IMAGE         Covariance of PSF position between x and y                 [pixel**2]
#ERRX2PSF_WORLD         Variance of PSF position along X-WORLD (alpha)             [deg**2]
#ERRY2PSF_WORLD         Variance of PSF position along Y-WORLD (delta)             [deg**2]
#ERRXYPSF_WORLD         Covariance of PSF position X-WORLD/Y-WORLD                 [deg**2]
#ERRCXXPSF_IMAGE        Cxx PSF error ellipse parameter                            [pixel**(-2)]
#ERRCYYPSF_IMAGE        Cyy PSF error ellipse parameter                            [pixel**(-2)]
#ERRCXYPSF_IMAGE        Cxy PSF error ellipse parameter                            [pixel**(-2)]
#ERRCXXPSF_WORLD        Cxx PSF error ellipse parameter (WORLD units)              [deg**(-2)]
#ERRCYYPSF_WORLD        Cyy PSF error ellipse parameter (WORLD units)              [deg**(-2)]
#ERRCXYPSF_WORLD        Cxy PSF error ellipse parameter (WORLD units)              [deg**(-2)]
#ERRAPSF_IMAGE          PSF RMS position error along major axis                    [pixel]
#ERRBPSF_IMAGE          PSF RMS position error along minor axis                    [pixel]
#ERRTHTPSF_IMAGE        PSF error ellipse position angle (CCW/x)                   [deg]
#ERRAPSF_WORLD          World PSF RMS position error along major axis              [pixel]
#ERRBPSF_WORLD          World PSF RMS position error along minor axis              [pixel]
#ERRTHTPSF_WORLD        PSF error ellipse pos. angle (CCW/world-x)                 [deg]
#ERRTHTPSF_SKY          Native PSF error ellipse pos. angle (east of north)        [deg]
#ERRTHTPSF_J2000        J2000 PSF error ellipse pos. angle (east of north)         [deg]
#ERRTHTPSF_B1950        B1950 PSF error ellipse pos. angle (east of north)         [deg]
#VECTOR_MODEL           Model-fitting coefficients                                
#VECTOR_MODELERR        Model-fitting coefficient uncertainties                   
#CHI2_MODEL             Reduced Chi2 of the fit                                   
#FLAGS_MODEL            Model-fitting flags                                       
#NITER_MODEL            Number of iterations for model-fitting                    
#FLUX_MODEL             Flux from model-fitting                                    [count]
#FLUXERR_MODEL          RMS error on model-fitting flux                            [count]
#MAG_MODEL              Magnitude from model-fitting                               [mag]
#MAGERR_MODEL           RMS error on model-fitting magnitude                       [mag]
#XMODEL_IMAGE           X coordinate from model-fitting                            [pixel]
#YMODEL_IMAGE           Y coordinate from model-fitting                            [pixel]
#XMODEL_WORLD           Fitted position along world x axis                         [deg]
#YMODEL_WORLD           Fitted position along world y axis                         [deg]
#ALPHAMODEL_SKY         Fitted position along right ascension  (native)            [deg]
#DELTAMODEL_SKY         Fitted position along declination (native)                 [deg]
#ALPHAMODEL_J2000       Fitted position along right ascension (J2000)              [deg]
#DELTAMODEL_J2000       Fitted position along declination (J2000)                  [deg]
#ALPHAMODEL_B1950       Fitted position along right ascension (B1950)              [deg]
#DELTAMODEL_B1950       Fitted position along declination (B1950)                  [deg]
#ERRX2MODEL_IMAGE       Variance of fitted position along x                        [pixel**2]
#ERRY2MODEL_IMAGE       Variance of fitted position along y                        [pixel**2]
#ERRXYMODEL_IMAGE       Covariance of fitted position between x and y              [pixel**2]
#ERRX2MODEL_WORLD       Variance of fitted position along X-WORLD (alpha)          [deg**2]
#ERRY2MODEL_WORLD       Variance of fitted position along Y-WORLD (delta)          [deg**2]
#ERRXYMODEL_WORLD       Covariance of fitted position X-WORLD/Y-WORLD              [deg**2]
#ERRCXXMODEL_IMAGE      Cxx error ellipse parameter of fitted position             [pixel**(-2)]
#ERRCYYMODEL_IMAGE      Cyy error ellipse parameter of fitted position             [pixel**(-2)]
#ERRCXYMODEL_IMAGE      Cxy error ellipse parameter of fitted position             [pixel**(-2)]
#ERRCXXMODEL_WORLD      Cxx fitted error ellipse parameter (WORLD units)           [deg**(-2)]
#ERRCYYMODEL_WORLD      Cyy fitted error ellipse parameter (WORLD units)           [deg**(-2)]
#ERRCXYMODEL_WORLD      Cxy fitted error ellipse parameter (WORLD units)           [deg**(-2)]
#ERRAMODEL_IMAGE        RMS error of fitted position along major axis              [pixel]
#ERRBMODEL_IMAGE        RMS error of fitted position along minor axis              [pixel]
#ERRTHETAMODEL_IMAGE    Error ellipse pos.angle of fitted position (CCW/x)         [deg]
#ERRAMODEL_WORLD        World RMS error of fitted position along major axis        [deg]
#ERRBMODEL_WORLD        World RMS error of fitted position along minor axis        [deg]
#ERRTHETAMODEL_WORLD    Error ellipse pos.angle of fitted position (CCW/world-x)   [deg]
#ERRTHETAMODEL_SKY      Native fitted error ellipse pos. angle (east of north)     [deg]
#ERRTHETAMODEL_J2000    J2000 fitted error ellipse pos. angle (east of north)      [deg]
#ERRTHETAMODEL_B1950    B1950 fitted error ellipse pos. angle (east of north)      [deg]
#X2MODEL_IMAGE          Variance along x from model-fitting                        [pixel**2]
#Y2MODEL_IMAGE          Variance along y from model-fitting                        [pixel**2]
#XYMODEL_IMAGE          Covariance between x and y from model-fitting              [pixel**2]
#E1MODEL_IMAGE          Ellipticity component from model-fitting                  
#E2MODEL_IMAGE          Ellipticity component from model-fitting                  
#EPS1MODEL_IMAGE        Ellipticity component (quadratic) from model-fitting      
#EPS2MODEL_IMAGE        Ellipticity component (quadratic) from model-fitting      
#CONCENTRATION_MODEL    Concentration parameter from model-fitting                
#CLASS_STAR_MODEL       S/G classifier from model-fitting                         
#FLUX_BACKOFFSET        Background offset from fitting                             [count]
#FLUXERR_BACKOFFSET     RMS error on fitted background offset                      [count]
#FLUX_SPHEROID          Spheroid total flux from fitting                           [count]
#FLUXERR_SPHEROID       RMS error on fitted spheroid total flux                    [count]
#MAG_SPHEROID           Spheroid total magnitude from fitting                      [mag]
#MAGERR_SPHEROID        RMS error on fitted spheroid total magnitude               [mag]
#SPHEROID_REFF_IMAGE    Spheroid effective radius from fitting                     [pixel]
#SPHEROID_REFFERR_IMAGE RMS error on fitted spheroid effective radius              [pixel]
#SPHEROID_REFF_WORLD    Spheroid effective radius from fitting                     [deg]
#SPHEROID_REFFERR_WORLD RMS error on fitted spheroid effective radius              [deg]
#SPHEROID_ASPECT_IMAGE  Spheroid aspect ratio from fitting                        
#SPHEROID_ASPECTERR_IMA RMS error on fitted spheroid aspect ratio                 
#SPHEROID_ASPECT_WORLD  Spheroid aspect ratio from fitting                        
#SPHEROID_ASPECTERR_WOR RMS error on fitted spheroid aspect ratio                 
#SPHEROID_THETA_IMAGE   Spheroid position angle (CCW/x) from fitting               [deg]
#SPHEROID_THETAERR_IMAG RMS error on spheroid position angle                       [deg]
#SPHEROID_THETA_WORLD   Spheroid position angle (CCW/world-x)                      [deg]
#SPHEROID_THETAERR_WORL RMS error on spheroid position angle                       [deg]
#SPHEROID_THETA_SKY     Spheroid position angle (east of north, native)            [deg]
#SPHEROID_THETA_J2000   Spheroid position angle (east of north, J2000)             [deg]
#SPHEROID_THETA_B1950   Spheroid position angle (east of north, B1950)             [deg]
#SPHEROID_SERSICN       Spheroid Sersic index from fitting                        
#SPHEROID_SERSICNERR    RMS error on fitted spheroid Sersic index                 
#FLUX_DISK              Disk total flux from fitting                               [count]
#FLUXERR_DISK           RMS error on fitted disk total flux                        [count]
#MAG_DISK               Disk total magnitude from fitting                          [mag]
#MAGERR_DISK            RMS error on fitted disk total magnitude                   [mag]
#DISK_SCALE_IMAGE       Disk scalelength from fitting                              [pixel]
#DISK_SCALEERR_IMAGE    RMS error on fitted disk scalelength                       [pixel]
#DISK_SCALE_WORLD       Disk scalelength from fitting (world coords)               [deg]
#DISK_SCALEERR_WORLD    RMS error on fitted disk scalelength (world coords)        [deg]
#DISK_ASPECT_IMAGE      Disk aspect ratio from fitting                            
#DISK_ASPECTERR_IMAGE   RMS error on fitted disk aspect ratio                     
#DISK_ASPECT_WORLD      Disk aspect ratio from fitting                            
#DISK_ASPECTERR_WORLD   RMS error on disk aspect ratio                            
#DISK_INCLINATION       Disk inclination from fitting                              [deg]
#DISK_INCLINATIONERR    RMS error on disk inclination from fitting                 [deg]
#DISK_THETA_IMAGE       Disk position angle (CCW/x) from fitting                   [deg]
#DISK_THETAERR_IMAGE    RMS error on fitted disk position angle                    [deg]
#DISK_THETA_WORLD       Disk position angle (CCW/world-x)                          [deg]
#DISK_THETAERR_WORLD    RMS error on disk position angle                           [deg]
#DISK_THETA_SKY         Disk position angle (east of north, native)                [deg]
#DISK_THETA_J2000       Disk position angle (east of north, J2000)                 [deg]
#DISK_THETA_B1950       Disk position angle (east of north, B1950)                 [deg]
#DISK_PATTERN_VECTOR    Disk pattern fitted coefficients                          
#DISK_PATTERNMOD_VECTOR Disk pattern fitted moduli                                
#DISK_PATTERNARG_VECTOR Disk pattern fitted arguments                              [deg]
#DISK_PATTERN_SPIRAL    Disk pattern spiral index  
"""

	# We turn this text block into a list of the parameter names:
	fullparamlist = list(map(lambda s: s[1:-1], re.compile("#\w*\s").findall(fullparamtxt)))

