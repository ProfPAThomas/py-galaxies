The galaxy formation model
==========================

Overview
--------

Halos
^^^^^

Halos are containers of diffuse, hot (virial temperature) gas **hot gas** (this is the intracluster gas, ICG)  and of subhalos and/or galaxies.  The gas is presumed predominantly ionised and therefore too hot for galaxy formation.

The may also contain **stars** (intracluster stars) that have been stripped from galaxies and/or subhalos.

Accretion
:::::::::

On each timestep, the baryonic content of each halo is checked and, if less than the global mean, it is topped up to that value.

Cooling
:::::::
  
* If there is a central subhalo, then the hotgas will cool onto it.
* If there are one or more subhalos but no central subhalo then the cluster is deemed to be in a dynamically turbulent state and there will be no cooling.
* If there are no subhalos then, for now, no cooling takes place: this may be relaxed in later versions of the code to allow galaxy formation to begin in lower mass halos.

Subhalos
^^^^^^^^

Subhalos are the dark mater halos within which galaxies form and reside.

Their virial temperature is typically lower than that of the host halo.  They contain **hot gas** (the galactic corona).

Galaxy merging
::::::::::::::

At the beginning of each timestep, check for merging of satellite galaxies onto the central galaxy.  This will trigger a rearrangement of both stars and cold gas within the central galaxy, dependent upon the change in angular momentum.  It may or may not require imposition of an explicit starburst phase, or that may arise naturally from the compression of cold gas -- to be determined.  The black hole of the satellite will merge onto that of the central galaxy.

Stripping
:::::::::

Satellite subhalos may have hot gas removed and transferred to the host halo either by tidal or ram pressure stripping.  Similarly, stars may be stripped from the central galaxy and/or whole satellite galaxies disrupted and transferred to the stellar component of the host halo.

Accretion
:::::::::

The subhalo, if it is a central subhalo, may have accreted hot gas from the host halo.

Cooling
:::::::

Hot gas will cool onto the central galaxy.

Galaxies
^^^^^^^^

Galaxies are comprised of **cold gas** (the interstellar medium, ISM), **disc** stars, **bulge** stars and a **black hole**.

Accretion
:::::::::

The cold gas may have accreted cooling gas from the containing subhalo.

Star Formation
::::::::::::::

The cold gas disc is broken up into annular rings.  Within each ring, star formation may occur (dependent upon the particular physical model for star formation, but typically requiring the gas to exceed a critical surface density threshold).  The stars will form within the stellar disc.

Star formation will trigger feedback of mass and metals from the cold gas into the hot gas of the enclosing subhalo.  In extreme cases, feedback may push material out of the subhalo and into the host halo.

Black hole accretion
::::::::::::::::::::

The black hole of the central galaxy will accrete some diffuse gas from the hot gas of the subhalo.  This will result in radio mode feedback of energy that pushes hot gas from the subhalo to the host halo.

Halos
-----

Properties
^^^^^^^^^^

In MEGA, galaxies do not come with a virial speed, needed to determine the virial temperature.  However, they have 2 measures that can be used for this purpose:

* 3D_velocity_dispersion, :math:`\sigma_3`
* mass / half_mass_radius

Consider the simple isothermal sphere (SIS) for which :math:`m=2\sigma^2r/G`, where :math:`m(r)` is the mass within radius :math:`r`, :math:`\sigma` is the (constant) 1-D velocity dispersion, and :math:`G` is the gravitational constant.
Then :math:`\sigma^2=Gm/2r`, where we can evaluate at any radius.

We then expect :math:`\sigma_3^2=3\sigma^2`.

Here is a plot showing the results of evaluating :math:`\sigma` these two different ways:

.. image:: figs/vdisp.png
   :width: 600
   :alt: halo 1-D velocity dispersions versus halo mass

The turn-down at small masses is almost certainly due to the finite softening affecting halos with small particle number.

The offset between the two methods is not yet clear: it could be due to the fact that the halos are not isothermal spheres, or it could be due to a bug / misunderstanding of the input data.  The degree of offset (a factor of approximately :math:`(3/2)^{1/2}`) corresponds to a density profile of :math:`\rho\propto r^{-1}` within the half mass radius, which seems rather shallow: this could be investigated by taking NFW profiles of different concentrations and solving the Jeans equations.

For now, we park this uncertainty and use the measured velocity dispersion of the halo as an indicator of its virial (ie the hot gas) temperature: :math:`k_\mathrm{B}T/\mu m_\mathrm{H}=\sigma^2`, where :math:`k_\mathrm{B}` is the Boltzmann constant and :math:`\mu m_\mathrm{H}\approx 10^{-27}` kg is the mass per particle in a fully ionised gas of cosmic metallicity.  That then gives:

.. image:: figs/temp.png
   :width: 600
   :alt: virial temperature versus halo mass

