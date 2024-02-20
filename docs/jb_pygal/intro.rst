Introduction
============

py-galaxies is a semi-analytic galaxy formation model within the L-Galaxies family: https://lgalaxiespublicrelease.github.io.

It features 3 major developments:

* A python front-end to make parameter selection and I/O much easier.
* The use of merger trees rather than a merger graph, to allow splitting of halos and subhalos.
* The formation of galaxies within subhalos, rather than halos, and the replacement of "Type 0" galaxies with galaxies within **central** subhalos, which are not guaranteed to exist, for example in merging systems.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage
   developers
   api

