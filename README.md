# DoubletFinder.jl
Given a high resolution input spectrum with columns of wavelength, flux, error in flux, the program returns possible doublets of ions of choice.

Astronomers study astrophysical sources by taking their spectra. From these spectra, a plethora of information about the astrophysical source and its environment can be gleaned. This module is meant to identify possible doublets of choice ions in a given spectra. In addition, it also gives the locations of these doublets, the corresponding velocities. The module also has functions within it which can be called for different purposes, see the revelant documenation for more information.

Shown here is an example of a high resolution spectrum ![high resolution spectrum](https://github.com/sameeresque/DoubletFinder.jl/blob/master/norm_spec.pdf)
