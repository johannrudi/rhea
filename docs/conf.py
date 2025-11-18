# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Folder containing this file
DOCS_FOLDER = os.path.basename(os.path.dirname(os.path.realpath(__file__)))

# Add utilities directory to path
sys.path.insert(0, os.path.abspath("./_utils"))

from doxygen_utils import doxygenfile_section

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Rhea"
copyright = "2025, Johann Rudi"
author = "Johann Rudi, Max Heldman, Leonid Pereiaslov, Jiaqi Fang, Jiashun Hu"
release = "2.0.0"

# Custom project information
github_namespace = "johannrudi"
github_repo = "rhea"
github_version = "main"
project_url = f"https://github.com/{github_namespace}/{github_repo}/"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # Sphinx's own extensions
    "sphinx.ext.autosectionlabel",  # Auto-generate section labels
    "sphinx.ext.todo",  # Enables `admonition-todo`
    # External
    "breathe",  # Integrates Doxygen with Sphinx
    "myst_parser",  # Markdown support
    "sphinx_copybutton",  # Copy button for code blocks
    "sphinx_design",  # Responsive web components
    "sphinx_inline_tabs",  # Inline tabbed content
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "_doxygen", "develop", "Thumbs.db", ".DS_Store"]

# Support both .rst and .md files
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Suppress the cache warning for myst_substitutions
suppress_warnings = ["config.cache"]

# -- MyST configuration ------------------------------------------------------
# https://myst-parser.readthedocs.io/en/latest/syntax/optional.html

myst_enable_extensions = [
    "amsmath",  # LaTeX amsmath package
    "attrs_inline",  # Add dditional information to elements
    "colon_fence",  # Colon syntax `:::{some_keyword}`
    "deflist",  # Definition lists
    "dollarmath",  # Insert math with `$` and `$$`
    "smartquotes",  # Convert to opening/closing quotes
    "substitution",  # Jinja templates
    "tasklist",  # Task lists with `[ ]` and `[x]`
]

# Define substitutions (variables and functions)
myst_substitutions = {
    "project": project,
    "project_url": project_url,
    "github_namespace": github_namespace,
    "github_repo": github_repo,
    "github_version": github_version,
    "docs_dirname": DOCS_FOLDER,
    "doxygenfile_section": doxygenfile_section,
}

# -- Breathe configuration ---------------------------------------------------
# https://breathe.readthedocs.io/en/latest/quickstart.html
# https://breathe.readthedocs.io/en/latest/directives.html#config-values

breathe_projects = {project: "./_doxygen/xml"}
breathe_default_project = project

breathe_show_define_initializer = True
breathe_show_enumvalue_initializer = True

# Specify domains for particular files according to their extension
breathe_domain_by_extension = {
    "h": "c",
    "c": "c",
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_title = project
html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
language = "en"

html_context = {
    # GitHub context
    "display_github": True,
    "github_user": github_namespace,
    "github_repo": github_repo,
    "github_version": github_version,
    "conf_py_path": f"/{DOCS_FOLDER}/",
}

# -- Furo configuration ------------------------------------------------------
# https://pradyunsg.me/furo/

html_theme_options = {
    "navigation_with_keys": True,
    "top_of_page_buttons": ["view", "edit"],
    "source_repository": project_url,
    "source_branch": github_version,
    "source_directory": f"{DOCS_FOLDER}/",
    "footer_icons": [
        {
            "name": f"{project} source code",
            "url": project_url,
            "html": """
<svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
</svg>
                """,
            "class": "",
        },
    ],
}
