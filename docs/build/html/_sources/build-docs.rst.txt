Building documentation
#############################

This documentation is built with Sphinx.

There has been a separate documentation environment to limit the dependencies in *SeidarT*. The YML file to build the documentation environment is in the root directory of the GitHub repo, *documentation-environment.yml*. For installing dependencies outside of conda, refer to this document to install the required packages. To build the documentation, ::
    
    conda env create -f documentation-environment.yml

This installs required dependencies to build the documentation from `Sphinx docstrings <https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html>`_ or for `Numpy docstrings <https://numpydoc.readthedocs.io/en/latest/format.html>`_. After building the environment, activate it::
     
    conda activate documentation

Using the *sphinx-apidoc* command you can generate the documentation. The source files are in *docs/source* and the code that we want to generate the autodocs fore is in *src/* so use the correct relative paths when running the command. For example, from the *docs* folder, run the command::
    
    sphinx-apidoc -o ./source ../src -e

Then, while still in the *docs* folder, generate the html and PDF files::
    
    make html 

This will build documentation and move these
items to *docs/build/* folder where GitHub will look to render them (``docs/``). You can preview any of the HTML files in the *docs/build/html/* folder in a web viewer like *Firefox*.

.. note:: 
    When contributing, use robust comments, type hints, docstrings, and include examples. Self-documenting code is very helpful, but to avoid any ambiguity, add comments. As the source code evolves, some comments will need to be updated. 

    
