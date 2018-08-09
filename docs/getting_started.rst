Getting Started
===============



Using RadioSource
-----------------

``RadioSource`` is the starting place for any analysis to be done using dendrocat. It takes a radio image, extracts information from the fits header, and sets attributes that are necessary for further processing.

.. code-block:: python

    >>> import dendrocat
    >>> from astrop
    y.io import fits
    >>> source_object = dendrocat.RadioSource(fits.open('/path/to/file.fits'))
    >>> source_object
    <dendrocat.radiosource.RadioSource at 0x7fda2f851080>

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
