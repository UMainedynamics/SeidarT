im2anim
#######################

*generates project text file template*

**Usage**

.. code-block:: bash
    
    im2anim -p PROJECTFILE  -c CHANNEL -n NUMSTEPS -d DELAY -a ALPHA -z -D -P PLOTTYPE

**Inputs**


* ``-p``: filename, including .prj extension

* ``-i``: .png file of cross sectional view of study site

    * No antialiasing
    * Use pixels dimensions proportional to unit dimensions
    * Start with small dimensions (150x150, 200x500), as larger
      dimensions quickly become more time intensive


**Outputs**

* Project file to fill in based on site characteristics

