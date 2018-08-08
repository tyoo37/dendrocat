dendrocat
--------------------------------------------------------------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

A package for detecting and processing sources in radio images. The core of$dendrocat contains infrastructure for source detection using $astrodendro, noisy source rejection, cataloging, and photometry. Auxiliary classes and methods provide more specialized analysis on already-photometered data.

Features
--------

- Source detection using $astrodendro
- A suite of customizable apertures, each with masking functionality
- Signal-to-noise calculation
- SNR-based source rejection
- Aperture photometry
- Automatically-generated source catalogs
- Flux-Flux plotting
- Spectral energy distribution plotting

Upcoming Features
-----------------

- Matching with external catalogs
- Aperture centering on peak flux
- Beam aperture and Beam Background aperture

Installation
------------
dendrocat requires the following packages to install:
* `numpy <http://www.numpy.org>`__ 1.14.3 or later
* `Astropy <http://www.astropy.org>`__ 1.2 or later

dendrocat requires the following packages to function:
* `regions <https://github.com/astropy/regions>`__ 0.3.dev662 or later
* `astrodendro <https://github.com/dendrograms/astrodendro>`__

$dendrocat can be cloned from its `repository <http://github.com/cmcclellan1010/dendrocat/>`__ on GitHub.

.. code-block:: bash
    
    git clone https://github.com/cmcclellan1010/dendrocat.git
    cd dendrocat
    python setup.py install

Contribute
----------

- Issue Tracker: github.com/cmcclellan1010/dendrocat/issues
- Source Code: github.com/cmcclellan1010/dendrocat

License
-------

This project is Copyright (c) B. Connor McClellan and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause license. See the licenses folder for
more information.
