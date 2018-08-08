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

If the FITS header is not recognized, the missing attributes will be printed to the console so they can be set manually. Typically, this will include

- ``nu`` : The central frequency of the observation. Specified as an ``astropy.units.quantity.Quantity``.
- ``freq_id`` : A string used to identify different images at a glance. For a 45 GHz observation, for example, the `freq_id` could be ``45.0GHz``.
- ``metadata`` : Image metadata required for generating a dendrogram.

.. code-block:: python

    >>> import astropy.units as u
    >>> source_object.nu = 45.0*u.GHz
    >>> source_object.freq_id = '45.0GHz'
    >>> source_object.set_metadata()
    
If the wcs, beam, and pixel scale are successfully extracted from the FITS header, after manually setting ``nu`` the method ``RadioSource.set_metadata`` can be called to fill in the rest of the metadata.

