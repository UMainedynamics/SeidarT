Creating a project file
~~~~~~~~~~~~~~~~~~~~~~~

A project file can be generating using *prjbuild*.

.. code-block::

    from seidart.routines.prjbuild import prjbuild
    prjbuild(*image_file_path*, *name_of_project*)


This command will create a blank project file and populate some of the domain fields. Alternatively, ``prjbuild`` can be called from the command line

.. code-block::

    prjbild -i IMAGE_FILE_PATH -o OUTPUT_JSON 

