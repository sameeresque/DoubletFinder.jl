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
Pkg.instantiate()
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


@everywhere const c_light=299792.458; # speed of light in km/s
@everywhere const LYA = 1215.6; # Lyman alpha rest-wavelength in angstrom


@everywhere begin
    """
        findV(zAbs, zEm)

    "Computes the velocity at a redshift, zAbs, relative to the emission redshift, zEm."

    # Examples
    ```julia-repl
    julia> findV(3.12,3.12)
    0.0
    ```
    """
    function findV(zAbs::Float64, zEm::Float64)
        @assert length(zAbs) == length(zEm)
        v = ((((1.0+zEm)^2) - ((1.0+zAbs)^2))/(((1.0+zEm)^2) + ((1.0+zAbs)^2))) * c_light
        return v
    end
end

@everywhere begin
    """
        findZAbs(v, zEm)

    "Computes the redshift at a given velocity, v, relative to the emission redshift, zEm."

    # Examples
    ```julia-repl
    julia> findZAbs(0, 3.12)
    3.12
    ```
    """
    function findZAbs(v::Float64, zEm::Float64)
        @assert length(v) == length(zEm)
        temp = ((1.0+zEm)^2)*((c_light-v)/(c_light+v))
        zAbs = temp^0.5 - 1.0
        return zAbs
    end
end

@everywhere begin
    """
        getblocks1(vs)

    "Returns an array of tuple_arrays that are indices of the input vector, vs, which have values less than 1"

    """
    function getblocks1(vs::Array{Float64,1})
        @assert length(vs)>0
        blocks = Tuple{Int, Int}[]
        inblock, start = false, 0, 0
        for (i, v) in enumerate(vs)
            if inblock
                if v >= 1.0
                    push!(blocks, (start, i-1))
                    inblock = false
                end
            else
                if v < 1.0
                    start = i
                    inblock = true
                end
            end
        end
        inblock && push!(blocks, (start, length(vs)))
        blocks
    end
end

@everywhere begin
    """
        getblocks2(vs)

    "Returns an array of tuple_arrays that are indices of the input vector, vs, which have values less than 1. A more Pythonic implementation."

    """
    function getblocks2(vs::Array{Float64,1})
        @assert length(vs)>0
        t = [false; vs .< 1.0; false]
        dt = diff(t)
        f = findall(==(1), dt)
        l = findall(==(-1), dt) .- 1
        collect(zip(f, l))
    end
end

@everywhere begin
    """
        getblocks2(vs)

    "Returns an array of tuples, where each tuple is separated by more than 1 element"

    """
    function sep_blocks(blocks::Array{Tuple{Int64,Int64},1})
        n_blocks=Array{Tuple{Int64,Int64},1}(undef,length(blocks))
        n_n_blocks=1
        for i in 1:length(blocks)
            if blocks[i][2]-blocks[i][1]>1
                n_blocks[n_n_blocks]=blocks[i]
                n_n_blocks=n_n_blocks+1
            end
        end
        resize!(n_blocks,n_n_blocks-1)
        return n_blocks
    end
end


@everywhere begin
    """
        threesigmadet(df,tuple_array)

    Returns true if the trough (absorption line) is a 3 sigma detection. Takes input as df, the data frame containing the spectral data, and the `tuple_array` which identifies the indices corresponding to beginning and ending of the trough.

    """
    function threesigmadet(df::DataFrame,tuple_array::Tuple{Int64,Int64})
        @assert size(df,2)>0
        @assert length(tuple_array)>0

        val=minimum(df[:FLUX][tuple_array[1]:tuple_array[2]])
        loc_min=findall(x->x==val,df[:FLUX])[1]
        if (1-val)/df[:ERROR][loc_min][1] >= 3
            return true
        else
            return false
        end
    end
end


@everywhere begin
    """
        atomicdata(all_ions)

    Returns a tuple consisting of dictionary of species in `all_ions`. The file `atoms_UVES.csv` is read 
    from the current path. It is a library containing information on commonly observed absorption lines in a quasar spectrum. It contains information such as rest-wavelength, oscillator strength, atomic mass. The file is ingested into DataFrame format for further manipulation. `all_ions` denotes the most commonly observed ions; it is an array string, e.g. ["NV","CIV","SiIV","OVI"]

    """
    function atomicdata(all_ions::Array{String,1}=["NV","CIV","SiIV","OVI","MgII","MgI","FeII","AlII","AlIII","MnII","CrII","ZnII","Lya","HI","NeVIII"])
        #species = Dict()
        # The data file containing atomic information of interest
        
        df_atomic_data=CSV.read("src/atoms_UVES.csv",header=["Ion","Rest_Wavelength","f","Gamma","Mass","Column6"],delim='\t')

        #delete unused column
        deletecols!(df_atomic_data, :Column6)

        #create a new column with only the Ion Identifier name
        df_atomic_data[:Identifier_Ion]=[split(i,r"[0-9]+")[1] for i in df_atomic_data[:Ion]]

        #create a new column with only the Wavelength Identifier name
        df_atomic_data[:Identifier_Wavelength]=[split(i,r"[a-zA-Z]+")[2] for i in df_atomic_data[:Ion]]

        #get the properties of all_ions from the atomic_data table
        species = Dict{String,Dict{String,Tuple{Float64,Float64}}}()
        for search_trans in all_ions
            #species[search_trans] = Dict()
            species[search_trans] = Dict{String,Tuple{Float64,Float64}}()
            for i in 1:1:length(df_atomic_data[1])
                atom = split(df_atomic_data[:Ion][i],r"[0-9]+")[1]
                tranwav = split(df_atomic_data[:Ion][i],r"[a-zA-Z]+")[2]
                current_oscillator = df_atomic_data[:f][i]
                if ( atom == search_trans )
                    species[atom][tranwav] = df_atomic_data[:Rest_Wavelength][i] , current_oscillator
                end
            end
        end

        for i in species
            species[i[1]]=sort(species[i[1]])
        end
        return species
    end
end

@everywhere begin
    """
        rb_limits(wavelength,zEm,species)

    Returns a tuple of blue limit and red limit of the spectrum. Takes as input a wavelength vector, emission redshift, and the dictionary of species.

    """
    function rb_limits(wavelength::Array{Float64,1},zEm::Float64,species::Dict{String,Dict{String,Tuple{Float64,Float64}}})
        # determining the blue limit
        @assert length(wavelength)>0
        @assert zEm>0
        if LYA * (1. + zEm) > wavelength[1]
           blue_limit = LYA * (1. + zEm)
        else
           blue_limit = wavelength[1]
        end
        #Emission redshift of quasar - 3000 km/s
        emission_limit = 1. + zEm - 0.01
        red_limit =  Dict{String,Float64}()
        for specie in species
             # Choose either the reddest_transition or the end of the spectrum
            reddest = maximum(species[specie[1]])[1]
            if ((species[specie[1]][reddest][1] * emission_limit) < wavelength[end])
                red_limit[specie[1]] = species[specie[1]][reddest][1] * emission_limit
            else
                red_limit[specie[1]] = wavelength[end]
            end
        end
        max_red = maximum(values(red_limit)); 
        return blue_limit, max_red
    end
end


@everywhere begin
    """
        df_interest(df,zEm,blue_limit,red_limit)

    Returns a DataFrame within the limits of interest. Takes as input the full DataFrame, emission redshift, and the limits.

    """
    function df_interest(df::DataFrame,zEm::Float64,blue_limit::Float64,max_red::Float64,
            species::Dict{String,Dict{String,Tuple{Float64,Float64}}})
        @assert max_red>blue_limit
        @assert zEm>0
        wavelengths=df[:WAVELENGTH]::Array{Float64,1}

        #dataframe of interest between the blue and red limits
        df_int=df[blue_limit.<= wavelengths .<=max_red,:]
        #adding the column for rest wavelength
        df_int[:Rest_Wavelength]=df_int[:WAVELENGTH]/(1.0+zEm);
        #create velocity columns appended to the dataframe for each of the species in the library
        for specie in species
            for transition in keys(specie[2])
                df_int[Symbol(join([specie[1],transition,"_Vel"]))]=c_light*(df_int[:Rest_Wavelength].^2 .- species[specie[1]][transition][1]^2)./(df_int[:Rest_Wavelength].^2 .+ species[specie[1]][transition][1]^2)
            end
        end
        return df_int
    end
end

@everywhere begin
    """
        datacleanandtag(df,telluriclines)

    Returns a DataFrame with telluric regions set to a value of 1 and data is cleaned for negative errors and fluxes.

    """
    function datacleanandtag(df::DataFrame,telluriclines::Array{Tuple{Float64,Float64},1}=[(9300.0,9630.0),(7594.0,7700.0),(6868.0,6932.0),(6277.0,6318.0)])
        # clean cosmic rays/bad high flux values
        df[:ERROR][df[:FLUX].>1.5].=0
        df[:FLUX][df[:FLUX].>1.5].=1
        for i in telluriclines
            # Avoid Telluric features, in order from strongest to weakest by setting those regions to a value of 1.
            df[:FLUX][(df[:WAVELENGTH].>i[1]) .& (df[:WAVELENGTH].<i[2])].=1;
        end
        # Remove negative errors
        df[:FLUX][df[:ERROR].<0] .= 1;
        df[:ERROR][df[:ERROR].<0] .= 0;
        # Remove negative flux
        df[:ERROR][df[:FLUX].<0] .= 0;
        df[:FLUX][df[:FLUX].<0] .= 0;
        # tag regions with flux less than 1
        df[:TAG]=df[:FLUX] .< 1.0;
        df[:MASK]=zeros(size(df[:TAG])[1]);

        df[:MASK][df[:TAG].==true] .=1;
        df[:MASK][df[:TAG].==false] .=0;    
        return df
    end
end

@everywhere begin
    """
        iscontained(array,value)

    Returns true if the inputted array contains the value of interest within it.

    """
    function iscontained(array::Array{Float64,1},value::Float64)
        if minimum(array) <= value <= maximum(array)
            return true
        else
            return false
        end
    end                                    
end


@everywhere begin
    """
        getINFO(df,zEm,blocks,indice,ion,transition,species_dictionary)

    Returns the indice, minimum velocity, maximum velocity, median redshift of the trough, an array of redshifts corresponding to components if a doublet is detected. Takes as input a dataframe, emission redshift, the tuple of troughs, index for the tuple describing a trough, ion of interest, transition of interest, and the `species_dictionary`
    """
    function getINFO(df::DataFrame,zEm::Float64,blocks::Array{Tuple{Int64,Int64},1},indice::Int64,ion::AbstractString,transition::String,species::Dict{String,Dict{String,Tuple{Float64,Float64}}}=atomicdata(df,zEm))

        d=df[blocks[indice][1]:blocks[indice][2],:]

    #     if length(d[:FLUX])>=5
    #         x=savitzkyGolay(Float64[ a for a in d[:FLUX] if !ismissing(a) ],5,1) #smoothing the array
    #     elseif length(d[:FLUX])>=3 & length(d[:FLUX])<5
    #         x=savitzkyGolay(Float64[ a for a in d[:FLUX] if !ismissing(a) ],3,1) #smoothing the array
    #     else
    #         x=d[:FLUX]
    #     end

        length(d[:FLUX])>=5 ? x=ss.savgol_filter(d[:FLUX],5,1) : x=ss.savgol_filter(d[:FLUX],3,1)

        xs=findlocalminima(x)
        yvals = [idx[1] for idx in xs]
        xvals=d[Symbol(join([ion,transition,"_Vel"]))][yvals]


        ele_list=Array{String}(undef,length(keys(species[ion])))

        n_ele_list=1

        for i in keys(species[ion])
            ele_list[n_ele_list]= i
            n_ele_list=n_ele_list+1
        end
        resize!(ele_list, n_ele_list-1)

        filter!(e->e != minimum(ele_list),ele_list)

        zlocs=Array{Float64}(undef, length(xvals)*length(ele_list))

        n_zlocs=1
        for j in xvals
            for ele in ele_list
                if length(findall(x->x==true,iscontained.([df[Symbol(join([ion,ele,"_Vel"]))][blocks[i][1]:blocks[i][2]] for i in 1:length(blocks)],j))) > 0
                    zlocs[n_zlocs]= findZAbs(-j, zEm) 
                    n_zlocs=n_zlocs+1
                end
            end
        end
        resize!(zlocs, n_zlocs-1)

        if transition==minimum(species[ion])[1]
            return (indice,blocks[indice],minimum(d[Symbol(join([ion,transition,"_Vel"]))]),maximum(d[Symbol(join([ion,transition,"_Vel"]))]),
                Statistics.median(d[:WAVELENGTH])/species[ion][transition][1]-1,zlocs)
        else
            return (indice,blocks[indice],minimum(d[Symbol(join([ion,transition,"_Vel"]))]),maximum(d[Symbol(join([ion,transition,"_Vel"]))]),
                Statistics.median(d[:WAVELENGTH])/species[ion][transition][1]-1,NaN)
        end
    end
end

@everywhere begin
    """
        possible_doublets_parallel(df,zEm,blocks,species_dictionary,search_ions)

    Returns an array of tuples with the possible doublets and the other information of indice location, maximum and minimum velocities of the troughs, and component redshifts. Parallel version.

    """
    function possible_doublets_parallel(df::DataFrame,zEm::Float64,blocks::Array{Tuple{Int64,Int64},1},dict_atomicdata::Dict{String,Dict{String,Tuple{Float64,Float64}}},search_ions::Dict{String,Tuple{Int64,Int64}}=Dict("NV" => (-70000,10000), "CIV" => (-70000,10000), "SiIV" => (-70000,10000),"OVI" => (-70000,10000)))
        doublet_array=Dict{String,Array{Tuple{Float64,Tuple{Int64,Int64},Float64,Float64,Float64,Array{Float64,1}},1}}()
        for specie in keys(search_ions)
            doublet_array[specie]=pmap(y->getINFO(df,zEm,blocks,y,specie,minimum(keys(dict_atomicdata[specie])),dict_atomicdata),1:length(blocks),batch_size=Int(round(length(blocks)/nworkers())))  
            doublet_array[specie] = [x for x in doublet_array[specie] if (x[3]>=search_ions[specie][1]) & (x[4]<=search_ions[specie][2]) & (length(x[6])!=0)]
        end
        return doublet_array
    end
end           

@everywhere begin
    """
        possible_doublets(df,zEm,blocks,species_dictionary,search_ions)

    Returns an array of tuples with the possible doublets and the other information of indice location, maximum and minimum velocities of the troughs, and component redshifts. Serial version.

    """
    function possible_doublets(df::DataFrame,zEm::Float64,blocks::Array{Tuple{Int64,Int64},1},dict_atomicdata::Dict{String,Dict{String,Tuple{Float64,Float64}}},search_ions::Dict{String,Tuple{Int64,Int64}}=Dict("NV" => (-70000,10000), "CIV" => (-70000,10000), "SiIV" => (-70000,10000),"OVI" => (-70000,10000)))
        doublet_array=Dict{String,Array{Tuple{Float64,Tuple{Int64,Int64},Float64,Float64,Float64,Array{Float64,1}},1}}()
        for specie in keys(search_ions)
            doublet_array[specie]=map(y->getINFO(df,zEm,blocks,y,specie,minimum(keys(dict_atomicdata[specie])),dict_atomicdata),1:length(blocks))  
            doublet_array[specie] = [x for x in doublet_array[specie] if (x[3]>=search_ions[specie][1]) & (x[4]<=search_ions[specie][2]) & (length(x[6])!=0)]
        end
        return doublet_array
    end        
end

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
