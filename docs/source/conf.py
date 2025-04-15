import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))  # Adjust path to your source code

# -- Project information -----------------------------------------------------

project = 'disccofan'
author = 'Simon Gazagnes, Michael H.F. Wilkilson'
release = '2.0.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = "images/disccofan.png"
html_theme_options = {
    'logo_only': True,
    'display_version': False,
}

from docutils import nodes
from docutils.parsers.rst import roles

def smallcaps_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    node = nodes.inline(text.upper(), text.upper())
    node['classes'].append('smallcaps')
    return [node], []

roles.register_canonical_role('sc', smallcaps_role)
def setup(app):
    app.add_css_file('custom.css')