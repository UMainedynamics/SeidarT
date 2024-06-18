Creating a project file
~~~~~~~~~~~~~~~~~~~~~~~

A project file can be generating using *prjbuild*.

.. :code-block::

    from seidart.routines.prjbuild import prjbuild
    
    prjbuild(*image_file_path*, *name_of_project*)

This command will create a blank project file and populate some of the domain fields. 
