``dendrocat``
=============

Introduction
------------

This is the documentation for ``dendrocat``, a package for detecting and processing sources in radio images. ``dendrocat`` provides classes and methods to do the following core functions:

- Pipeline radio images to source catalogs
- Match sources between multiple catalogs
- Aperture photometry on source catalogs

Note that this package relies on other astronomy-related Python packages to function, which themselves may be in development.

- Code: `Github repository <https://github.com/cmcclellan1010/dendrocat>`__
- Docs: ``dendrocat`` documentation
- Contributors: https://github.com/cmcclellan1010/dendrocat/graphs/contributors

Getting Started
---------------

The procedure for creating radio source objects, generating dendrograms, and creating source catalogs is demonstrated below. This example uses default settings--more in-depth documentation is located in :ref:`the_radiosource_class` and :ref:`the_mastercatalog_class` sections.

``RadioSource`` is the starting place for any analysis to be done using dendrocat. It takes a radio image, extracts information from the FITS header, and sets attributes that are necessary for further processing.

.. code-block:: python

    >>> import dendrocat
    >>> from astropy.io import fits
    >>> source_object = dendrocat.RadioSource(fits.open('/path/to/file.fits'))
    >>> source_object
    <dendrocat.radiosource.RadioSource at 0x7fda2f851080>

.. Note:: FITS header extraction is undergoing development. Currently, only specific headers from EVLA and ALMA are supported, but this will change soon. See :ref:`the_radiosource_class` for details.

A catalog of sources can then be created using dendrograms.

.. code-block:: python

    >>> source_object.to_catalog()
    Computing catalog for 113 structures
    [========================================>] 100%
    <Table masked=True length=113>
    _idx _index _name  ... rejected 226.1GHz_detected
     ...    ...   ...  ...      ...              ...
     
    >>> source_object.catalog
    <Table masked=True length=113>
    _idx _index _name  ... rejected 226.1GHz_detected
     ...    ...   ...  ...      ...              ...
     
``to_catalog`` returns a source catalog, which is also saved as instance attribute. If ``to_catalog`` is called and the ``RadioSource`` object has no existing dendrogram, one will be generated using default parameters.

Using ``dendrocat``
-------------------

More in-depth documentation is linked here.

The RadioSource Class
~~~~~~~~~~~~~~~~~~~~~
.. toctree::
  :maxdepth: 2
  
  the_radiosource_class

The MasterCatalog Class
~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   the_mastercatalog_class

Reference/API
~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
    
   api
   
   
