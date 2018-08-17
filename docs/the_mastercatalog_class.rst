.. _the_mastercatalog_class:

The MasterCatalog Class
=====================

A `~dendrocat.MasterCatalog` object combines multiple `~dendrocat.RadioSource` objects, retaining each object's attributes while creating a new, combined catalog from the sources contained in their individual catalogs.

The `~dendrocat.MasterCatalog` enables analysis not possible with a single `~dendrocat.RadioSource` object, such as photometry and catalog cross-matching. 

Adding Objects
--------------

Additional `~dendrocat.RadioSource` or `~dendrocat.MasterCatalog` objects can be added to an existing `~dendrocat.MasterCatalog` using `~dendrocat.MasterCatalog.add_objects`.

.. code-block:: python
    from dendrocat import RadioSource, MasterCatalog
    from dendrocat.utils import match
    from astropy.io import fits

    source_object1 = RadioSource(fits.open('file1.fits'), name='so1')
    source_object2 = RadioSource(fits.open('file2.fits'), name='so2')
    source_object3 = RadioSource(fits.open('file3.fits'), name='so3')

If ``source_object3`` is a much lower-resolution, noisier image, you may want to forgo generating a dendrogram for it. The source regions from the other two higher-resolution images may be used instead.

.. code-block:: python

    source_object1.autoreject()
    source_object2.autoreject()

    mastercatalog = match(source_object1, source_object2)

Now, adding other `~dendrocat.RadioSource` objects or `~dendrocat.MasterCatalog` objects will preserve the existing `~dendrocat.MasterCatalog`'s source catalog, while adding all the associated `~dendrocat.RadioSource` objects to its dictionary of attributes.

.. code-block:: python

    >>> mc.add_objects(source_object3)
    >>> mc.__dict__.keys()
    dict_keys(['catalog', 'accepted', 'so1', 'so2', 'so3'])

At this point, performing photometry yields photometry data for all three images, though only two images were used to detect the sources in the first place.

.. Note:: The `~dendrocat.MasterCatalog` that calls `~dendrocat.MasterCatalog.add_objects` will always have its catalog preserved. If `~dendrocat.MasterCatalog.add_objects` adds another `~dendrocat.MasterCatalog`, all of the other `~dendrocat.MasterCatalog`'s associated `~dendrocat.RadioSource` objects will be appropriated by the `~dendrocat.MasterCatalog` calling `~dendrocat.MasterCatalog.add_objects`. The other's source catalog will be abandoned.

Renaming Sources
----------------

Any source can have its ``_name`` set for distinguishability. Suppose we have a `~dendrocat.MasterCatalog` object named ``mc``.

.. code-block:: python

    #Index of the row in the source catalog containing the old name
    index_of_old_name = [mc.catalog['_name'] == '226007']
    
    # Assign new name to the row
    mc.catalog['_name'][index_of_old_name] = 'w51d2'

This can also be done in one line.

.. code-block:: python
    
    >>> mc.catalog['_name'][mc.catalog['_name'] == '226007'] = 'w51d2'

Here, we first access the ``_name`` column of the catalog. Then, we index the column with the location in the catalog the ``_name`` is equal to the old name, ``'226007'``. Then the entry in the catalog is set to be equal to ``'w51d2'``
