Installing sewpy
----------------

The requirements are `astropy <http://www.astropy.org>`_ (version 0.4.2 or later), and obviously `SExtractor <http://www.astromatic.net/software/sextractor>`_ (``sewpy`` does not require a specific version).


To get ``sewpy``, clone the `project on github <http://github.com/megalut/sewpy>`_, or click `here <https://github.com/megalut/sewpy/tarball/master>`_ to get a tarball of the latest version. You now have (at least) two options to install the module.

Method 1, if you just want to use the module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Extract the tarball if required, ``cd`` into it, and

.. code-block:: bash

	python setup.py install --user


For a system-wide install, remove the ``--user``. More information about how to install python modules with distutil can be found `in the official Distutils documentation <https://docs.python.org/2/install/index.html#install-index>`_.

Method 2, if you intend to develop the module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simple solution is to add your local clone of ``sewpy`` to your python path, for instance by modifying your ``.bash_profile`` or equivalent:

.. code-block:: bash

	export PYTHONPATH=${PYTHONPATH}:/path/to/megalut-sewpy-1.0

