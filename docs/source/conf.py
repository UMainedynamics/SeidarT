from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SeidarT'
copyright = '2024, Steven Bernsen'
author = 'Steven Bernsen'
version = '3.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


# extensions = ['sphinx.ext.autodoc', ]
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "sphinx.ext.graphviz",
    "sphinx_copybutton",
    "myst_parser",
]

templates_path = ['_templates']
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

autodoc_mock_imports = [
    "dill",
    "disba",
    "glob2",
    "matplotlib",
    "mpl_toolkits",
    "mplstereonet",
    "PIL",
    "pyevtk",
    "scipy",
    "seaborn",
    "seidart.routines.prjrun",
    "trimesh",
]

suppress_warnings = ["docutils", "myst.header"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'haiku'
# html_static_path = ['_static']

html_theme = "furo"
html_title = f"{project} documentation v{version}"

# Disable the generation of the various indexes
html_use_modindex = False
html_use_index = False
