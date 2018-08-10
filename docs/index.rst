.. _dendrocat_documentation:

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
- Docs: :ref:`dendrocat_documentation`
- Contributors: https://github.com/cmcclellan1010/dendrocat/graphs/contributors

Getting Started
---------------

The procedure for creating radio source objects, generating dendrograms, and creating source catalogs is demonstrated below. This example uses default settings. In-depth documentation is located in :ref:`the_radiosource_class` and :ref:`the_mastercatalog_class` sections.

``RadioSource`` is the starting place for any analysis to be done using dendrocat. It takes a radio image, extracts information from the FITS header, and sets attributes that are necessary for further processing.

.. code-block:: python

    >>> from dendrocat import RadioSource
    >>> from astropy.io import fits
    >>> source_object = RadioSource(fits.open('/path/to/file.fits'))
    >>> source_object
    <dendrocat.radiosource.RadioSource at 0x7fda2f851080>

.. Note:: FITS header extraction is undergoing development. Currently, only specific headers from EVLA and ALMA are supported, but this will change soon. See :ref:`the_radiosource_class` for details.

A catalog of sources can then be created using dendrograms.

.. code-block:: python

    >>> source_object.to_catalog()
    Generating dendrogram using 738,094 of 49,787,136 pixels (1.4824994151099593% of data)
    [========================================>] 100%
    Computing catalog for 113 structures
    [========================================>] 100%
    <Table masked=True length=113>
    _idx _index _name  ... rejected 226.1GHz_detected
     ...    ...   ...  ...      ...              ...

Most false detections (i.e., noise) can be filtered out using the `~dendrocat.RadioSource.autoreject` method, with the signal-to-noise requirement specified as a keyword argument ``threshold``. To view all the image cutouts around each source in a grid, call the `~dendrocat.RadioSource.plot_grid` method. Rejected sources are greyed out, and the source apertures used to calculate the signal to noise are overplotted on top of the image.

.. code-block:: python

    >>> source_object.autoreject(threshold=6.)
    >>> source_object.plot_grid(skip_rejects=False)
    
.. image:: ./_figures/plot_grid_skip_rejects_false.png
    :width: 400
    :alt: A grid of extracted sources from a radio image, showing overlaid elliptical and annular apertures. Some of the squares in the grid are greyed out to represent which sources are rejected.
    
Sources can then be manually accepted or rejected using `~dendrocat.RadioSource.accept` or `~dendrocat.RadioSource.reject` (and reset entirely using `~dendrocat.RadioSource.reset`). 

.. code-block:: python

    >>> source_object.reject([226000, 226008, 226016, 226023, 226028, 226085, 226124, 226135, 226137])
    >>> source_object.accept([226024, 226043, 226123])
    >>> source_object.plot_grid(skip_rejects=False)

.. image:: ./_figures/plot_grid_skip_rejects_false_manual.png
    :width: 400
    :alt: A grid of extracted sources from a radio image, showing overlaid elliptical and annular apertures. Some of the squares in the grid are greyed out to represent which sources are rejected.
    
The identifiers used to accept and reject sources are those stored under the ``_name`` column in the `~dendrocat.RadioSource` object's catalog. More information about accessing Astropy tables can be found in the `Accessing a table <http://docs.astropy.org/en/stable/table/access_table.html>`__ section of the Astropy documentation.

`dendrocat.utils.saveregions` can be used to save a DS9 region file of all the apertures in a catalog. With an image open in DS9, the region file can be loaded in to check each aperture and find the proper identifier (saved as ``_name`` in the source catalog) to use for manual acceptance or rejection.

.. code-block:: python

    >>> dendrocat.utils.saveregions(source_object.catalog, '/path/to/outputfile.reg')

The `~dendrocat.RadioSource` object is (more or less) complete after source rejection is finished. It can be combined with other `~dendrocat.RadioSource` objects to form a `~dendrocat.MasterCatalog`. The catalogs of all the `~dendrocat.RadioSource` objects are combined using `dendrocat.utils.match`.

.. code-block:: python

    from astropy.io import fits
    from dendrocat import RadioSource, MasterCatalog
    from dendrocat.utils import match
    
    source_object1 = RadioSource(fits.open('/path/to/file1.fits'))
    source_object2 = RadioSource(fits.open('/path/to/file2.fits'))

    
    source_object1.autoreject()
    source_object2.autoreject()
    
    combined_catalog = match(source_object1.catalog, source_object2.catalog)

    mastercat_object = MasterCatalog(source_object1, source_object2, catalog=combined_catalog)
    




Using ``dendrocat``
-------------------

In-depth documentation is linked below.

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

Apertures
~~~~~~~~~

.. toctree::
   :maxdepth: 2

   apertures

Reference/API
~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
    
   api
   
   
