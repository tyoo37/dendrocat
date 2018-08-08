Getting Started
===============

Introduction
------------

``dendrocat`` provides classes and methods to do the following:

- Radio image to source catalog pipeline
- Matching sources between multiple catalogs
- Aperture photometry

Note that this package relies on other astronomy-related Python packages to function, which themselves may be in development.

The RadioSource Class
---------------------

The ``RadioSource`` class is the starting place for any analysis to be done using ``dendrocat``. It takes a radio image, extracts information from the fits header, and sets attributes that are necessary for further processing. 

.. code-block:: python
    
    >>> import dendrocat
    >>> from astropy.io import fits
    >>> source_object = dendrocat.RadioSource(fits.open('/path/to/file.fits'))

Initializing a ``RadioSource`` object will attempt to extract data from the FITS file header of the image. Since FITS headers vary between telescopes and observations, support for new header formats will be added as needed.

Currently, ``dendrocat`` supports FITS header information extraction from the following telescopes:

- EVLA
- ALMA

*Developer note: In a future update, this functionality will be replaced by a ``fitsconfig`` file that will tell the program where to find the information it needs. Users can add to their ``fitsconfig`` to handle unsupported formats.*



