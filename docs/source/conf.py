# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import pathlib
import sys
subfolder = pathlib.Path(__file__).parents[2].resolve()/'src/gfapy'
sys.path.insert(0, subfolder.as_posix())

project = 'iGFA'
copyright = '2024, Shriramprasad V, Rudiyanto Gunawan'
author = 'Shriramprasad V, Rudiyanto Gunawan'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration', 
    'sphinx.ext.doctest', 
    'sphinx.ext.autodoc', 
    'sphinx.ext.autosummary'
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_theme_options = {
    'page_width': 'auto',
}
html_static_path = ['_static']
html_css_files = [
    'custom.css'  # Include your custom CSS file
]
