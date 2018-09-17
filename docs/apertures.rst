.. apertures:

Apertures
=========

Using Preset Apertures
----------------------

The presets for an `~dendrocat.aperture.Aperture` are `~dendrocat.aperture.Ellipse`, `~dendrocat.aperture.Circle`, and `~dendrocat.aperture.Annulus`. Each of these subclasses inherit properties from the base class, `~dendrocat.aperture.Aperture`.

When using `~dendrocat.MasterCatalog.photometer`, an aperture must be specified. To use variable-width apertures that change with the FWHM of the sources in the catalog, use any of the `~dendrocat.aperture.Aperture` subclasses as the ``aperture`` argument in `~dendrocat.MasterCatalog.photometer`.

For this example, suppose we have a `~dendrocat.MasterCatalog` object called ``mc`` (For directions on how to create this object, see the :ref:`getting_started` section) and we want to do some photometry.

.. code-block :: python

    from dendrocat.aperture import Ellipse, Circle, Annulus

    mc.photometer(Ellipse, Circle, Annulus)

This will take the `~dendrocat.MasterCatalog`'s catalog source entries and use their major and minor FWHM as aperture radii, where applicable. 
 - For a `~dendrocat.aperture.Circle`, the radius is the source's major FWHM.
 - For a `~dendrocat.aperture.Ellipse`, the major and minor FWHM of the source, as well as its position angle, are used directly as the parameters of the elliptical aperture.
 - For a `~dendrocat.aperture.Annulus`, the inner and outer radii of the aperture are determined by the major FWHM of the source, as well as the ``annulus_padding`` and ``annulus_width`` attributes of the `~dendrocat.RadioSource` object storing the image data. 

Annulus Inner Radius = Source Major FWHM + Annulus Padding

Annulus Outer Radius = Source Major FWHM + Annulus Padding + Annulus Width

These two parameters ensure that annular apertures don't overlap with source apertures, and can be tuned within each `~dendrocat.RadioSource` object.

Defining Custom Apertures
-------------------------

To use fixed apertures instead, create an instance of any of the `~dendrocat.aperture.Aperture` subclasses with the desired parameters. Parameters can be specified in pixel or degree coordinates.

.. code-block:: python

    import astropy.units as u

    # Define a fixed-radius elliptical aperture in pixels
    fixed_ellipse_pix = Ellipse([0,0], 15*u.pix, 10*u.pix, 30*u.deg, name=ellipsepix)

    # Define a fixed-radius elliptical aperture in degrees
    fixed_ellipse_deg = Ellipse([0,0], 15*u.arcsec, 10*u.arcsec, 30*u.deg, name=ellipsedeg)

Photometry can then be performed exactly as if these were new aperture presets.

.. code-block:: python

    mc. photometer(fixed_ellipse_pix, fixed_ellipse_deg)

.. note::

    The first argument of any `~dendrocat.aperture.Aperture` subclass is always ``center``. When creating an instance of the `~dendrocat.aperture.Aperture` subclasses, this argument can be filled with any two coordinates---they will be overwritten with the source objects' center coordinates when photometry is performed.
