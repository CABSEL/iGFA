# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import pathlib
import sys
code_folder = pathlib.Path(__file__).parents[2].resolve()/'code/src'
sys.path.insert(0, code_folder.as_posix())

# subfolder = pathlib.Path(__file__).parents[2].resolve()/'code/src/gfapy'
# sys.path.insert(0, subfolder.as_posix())

tutorials_folder = pathlib.Path(__file__).parents[2].resolve()/'tutorials'
sys.path.insert(0, tutorials_folder.as_posix())

project = 'iGFA'
copyright = 'MIT License'
author = 'Shriramprasad V, Rudiyanto Gunawan'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration', 
    'sphinx.ext.doctest', 
    'sphinx.ext.autodoc', 
    'sphinx.ext.autosummary',
    'sphinx.ext.githubpages',
    'sphinx.ext.todo',
    'myst_nb'
]
autosummary_generate = True
autosummary_ignore_module_all = False
templates_path = ['_templates']
exclude_patterns = []
todo_include_todos = True
nb_execution_mode = "off"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_theme_options = {
    'page_width': 'auto',
}
html_static_path = ['_static']
html_logo = "_static/igfa_logo.png"
html_favicon = '_static/igfa_logo.png'
html_css_files = [
    'custom.css'  # Include your custom CSS file
]


# Example configuration for intersphinx: refer to the Python standard library.
# intersphinx_mapping = {'tutorial_nb': ('https://github.com/CABSEL/glycoTS', None)}
# intersphinx_disabled_reftypes = ["*"]