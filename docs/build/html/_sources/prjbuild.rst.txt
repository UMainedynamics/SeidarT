prjbuild
#######################

*generates project text file template*

**Usage**

.. code-block:: bash

    prjbuild -h -i IMAGEFILE -p PROJECTFILE

**Inputs**

* ``-h,--help``

* ``-i,--imagefile``: .png file of cross sectional view of study site

    * No antialiasing
    * Use pixels dimensions proportional to unit dimensions
    * Start with small dimensions (150x150, 200x500), as larger
      dimensions quickly become more time intensive

* ``-o,--output_json``: filename, including .json extension


**Outputs**

* Project file to fill in based on site characteristics

