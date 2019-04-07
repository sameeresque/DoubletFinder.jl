module DoubletFinder

import Pkg
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.instantiate()
@everywhere using Distributions
@everywhere include("./linefinder.jl")

@everywhere using Distributed
@everywhere using CSV
@everywhere using PyCall
@everywhere using Test
@everywhere using Statistics
@everywhere using DataFrames, Query, DataFramesMeta
@everywhere using Images
@everywhere using LinearAlgebra
@everywhere using DSP
@everywhere using Printf
@everywhere using Random
@everywhere @pyimport scipy.signal as ss

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
export DoubletFinder

                        
end #module