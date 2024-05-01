import os
import sys
sys.path.insert(0, os.path.abspath('../src/seidart'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SeidarT'
copyright = '2024, Steven Bernsen'
author = 'Steven Bernsen'
version = '0.1.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


# extensions = ['sphinx.ext.autodoc', ]
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "recommonmark",
]

templates_path = ['_templates']
exclude_patterns = []

from recommonmark.parser import CommonMarkParser

source_parser = {
    '.md': CommonMarkParser
}

source_suffix = ['.rst', '.md']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'haiku'
# html_static_path = ['_static']

html_theme = "furo"
html_title = f"{project} documentation v{version}"

# Disable the generation of the various indexes
html_use_modindex = False
html_use_index = False