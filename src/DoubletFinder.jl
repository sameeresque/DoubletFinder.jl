module DoubletFinder

import Pkg
import PyCall

# export functions

export findV
export findZAbs
export getblocks1
export getblocks2
export sep_blocks
export three_sigmadet
export atomicdata
export rb_limits
export df_interest
export datacleanandtag
export iscontained
export getINFO
export possible_doublets_parallel
export possible_doublets
export Find


using Distributed
@everywhere using Distributed
using Pkg
@everywhere using Pkg
Pkg.activate(".")
@everywhere Pkg.activate(".")
#@everywhere Pkg.instantiate()
using CSV
@everywhere using CSV
using PyCall
@everywhere using PyCall
using Test
@everywhere using Test
using Statistics
@everywhere using Statistics
using DataFrames, Query, DataFramesMeta
@everywhere using DataFrames, Query, DataFramesMeta
using Images
@everywhere using Images
using LinearAlgebra
@everywhere using LinearAlgebra
using DSP
@everywhere using DSP
using Printf
@everywhere using Printf
ss = pyimport("scipy.signal")
@everywhere ss = pyimport("scipy.signal")
using Distributions
@everywhere using Distributions
@everywhere include("src/linefinder.jl")

@everywhere const c_light=299792.458; # speed of light in km/s
@everywhere const LYA = 1215.6; # Lyman alpha rest-wavelength in angstrom


@everywhere begin
    """
    Find(file,zEm,search_ions)

    Returns an array of tuples with the possible doublets and the other information of indice location, maximum and minimum velocities of the troughs, and component redshifts. 
    """
    function Find(file::String,zEm::Float64,search_ions::Dict{String,Tuple{Int64,Int64}})
        df_spectral_data = CSV.read(file,header=true,delim='\t',allowmissing=:none);
        # dictionary of species
        species_dict=atomicdata();
        # red and blue limits 
        rb_limit=rb_limits(df_spectral_data[:WAVELENGTH],zEm,species_dict);
        # DataFrame of interest
        df_int=df_interest(df_spectral_data,zEm,rb_limit[1],rb_limit[2],species_dict);
        # Cleaned DataFrame of interest
        df_clean=datacleanandtag(df_int);
        # Locations of possible absorption troughs
        blocks=getblocks1(df_clean[:FLUX]);
        blocks_3s=sep_blocks(blocks)[pmap((x->threesigmadet(df_clean,x)),sep_blocks(blocks),batch_size=Int(round(length(sep_blocks(blocks))/nworkers())))];
        return possible_doublets_parallel(df_clean,zEm,blocks_3s,species_dict,search_ions)
    end
end

end #module