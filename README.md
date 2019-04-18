# DoubletFinder.jl
Given a moderate/high resolution input spectrum with columns of wavelength, flux, error in flux, the program returns possible doublets of ions of choice.

Astronomers study astrophysical sources by taking their spectra. From these spectra, a plethora of information about the astrophysical source and its environment can be gleaned. This module is meant to identify possible doublets of choice ions in a given spectra. In addition, it also gives the locations of these doublets, the corresponding velocities and, the module also has functions within it which can be called for different purposes, see the revelant documenation for more information.

Shown here is an example of a ![high resolution spectrum](https://github.com/sameeresque/DoubletFinder.jl/blob/master/examples/norm_spec.pdf). We would like to find the ![needles](https://github.com/sameeresque/DoubletFinder.jl/blob/master/examples/Merged_doublets.pdf) in the ![haystack](https://github.com/sameeresque/DoubletFinder.jl/blob/master/examples/norm_spec.pdf). This is precisely the purpose of the module: to find potential doublets present in a given spectrum.

The module has been constructed to utilize parallel computing capabilities if available at your disposal. The following graph demonstrates the performance comparison on an `Intel(R) Xeon(R) CPU E5-2680 v2 @ 2.80GHz` with different number of cores and for different number of species.

![performance](https://github.com/sameeresque/DoubletFinder.jl/blob/master/examples/comparison.png)

The plot on the left shows the computation times for different number of species for different number of processors. On the right, we see the expected times assuming an ideal case, where the expected time with p processors is T/p, T being time for computation on a single processor. The difference between the actual and expected times, the overhead, is because of the communication between the master and the slaves, load unbalance, and also partly due to the fact that some portion of the code is inherently sequential. 

# How to setup 

From the command prompt issue 
```
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/sameeresque/DoubletFinder.jl.git"))'
```
If you have access to multiple processors, n, you can start your Julia session like so

```
julia -p n
```

# How to use the package 
- Demonstrated using an example: `raw_data.txt` is the spectrum file (csv format) containing columns of wavelength, flux, error in flux. emission_redshift is the emission redshift of the quasar, and the third argument is a dictionary of species and the search windows in km/s.
```
using DoubletFinder
@Find("raw_data.txt",emission_redshift=2.28,Dict("SiIV"=>(-70000,10000),"CIV"=>(-70000,10000)))
```

# Other functions available within this module

There are other functions that are available within this module and they can be invoked by using for e.g., `DoubletFinder.getblocks1`. To know more about what each function does access the documentation using for e.g.,
`?getblocks1`

![other functions](https://github.com/sameeresque/DoubletFinder.jl/blob/master/examples/Screenshot%20from%202019-04-12%2016-01-10.png)


