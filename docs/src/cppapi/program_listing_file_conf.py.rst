
.. _program_listing_file_conf.py:

Program Listing for File conf.py
================================

|exhale_lsh| :ref:`Return to documentation for file <file_conf.py>` (``conf.py``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: py

   # Configuration file for the Sphinx documentation builder.
   #
   # This file only contains a selection of the most common options. For a full
   # list see the documentation:
   # http://www.sphinx-doc.org/en/master/config
   
   # -- Path setup --------------------------------------------------------------
   
   # If extensions (or modules to document with autodoc) are in another directory,
   # add these directories to sys.path here. If the directory is relative to the
   # documentation root, use os.path.abspath to make it absolute, like shown here.
   #
   import os
   import textwrap
   import sys
   sys.path.insert(0, os.path.abspath('..'))
   
   print(sys.path)
   
   # -- Project information -----------------------------------------------------
   
   project = 'angle_calibration'
   copyright = '2024, CPS Detector Group'
   author = 'CPS Detector Group'
   version = '0.1.0'
   
   # -- General configuration ---------------------------------------------------
   
   # Add any Sphinx extension module names here, as strings. They can be
   # extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
   # ones.
   extensions = ['breathe',
                 'exhale',
                 'sphinx.ext.autodoc',
                 'sphinx.ext.napoleon',
                 'sphinx.ext.autosummary'
                ]
   
   breathe_default_project = "angle_calibration"
   napoleon_use_ivar = True
   
   master_doc = 'index'
   
   # Add any paths that contain templates here, relative to this directory.
   #templates_path = ['_templates']
   
   # List of patterns, relative to source directory, that match files and
   # directories to ignore when looking for source files.
   # This pattern also affects html_static_path and html_extra_path.
   #exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
   
   
   # -- Options for HTML output -------------------------------------------------
   
   # The theme to use for HTML and HTML Help pages.  See the documentation for
   # a list of builtin themes.
   #
   html_theme = "furo"
   
   # Add any paths that contain custom static files (such as style sheets) here,
   # relative to this directory. They are copied after the builtin static files,
   # so a file named "default.css" will overwrite the builtin "default.css".
   #html_static_path = ['static']
   
   
   #def setup(app):
       #app.add_css_file('css/extra.css')  # may also be an URL
   
   
   # Exhale config
   
   exhale_args = {
       "containmentFolder": "/home/mazzol_a/Documents/angle_calibration/docs/src/cppapi",
       "rootFileName": "generated_cpp_api.rst",
       "rootFileTitle": "C++ API",
       "doxygenStripFromPath": "../",  # remove common prefix from paths
   
       # Tells Exhale to run Doxygen automatically
       #"createDoxygenXml": False,  # set to True if you want Exhale to run Doxygen for you
   
       # Breathe config
       #"breathe_projects_source": {
           #"angle_calibration": "../include"  # or src/, as appropriate
       #},
       #"breathe_default_project": "angle_calibration",
   
       # Reorganize output for nicer layout
       "exhaleExecutesDoxygen": False,
       "exhaleUseDoxyfile": True,
   }
