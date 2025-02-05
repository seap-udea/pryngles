# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../src/pryngles'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Pryngles'
copyright = '2025, Jorge I. Zuluaga, Allard Veenstra, Sebastian Numpaque, Mario Sucerquia and Jaime A. Alvarado-Montes'
author = 'Zuluaga et al.'
release = '0.9.10'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [#'nbsphinx', 
              'sphinx_mdinclude', 'sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.viewcode']

# nbsphinx_execute = 'never'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygment_style = 'sphinx'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
