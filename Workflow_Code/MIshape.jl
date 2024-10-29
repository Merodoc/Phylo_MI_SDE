#Check the shape of the Multiple Imputed data to see how reasonable it is

using DataFrames
using Phylo
using CSV
using StatsKit
using Plots
using Statistics

dir = "C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/Data/"

extensions = ["1lPMMnoPVR/", "Cart/", "impnopvrpmm/", "impnopvrpan/", "Midas/", "PMMPVR/"]

function Phybridge_MIDict(dir)
    try 
        readdir(dir)
    catch
        println("Invalid directory")
    end
    Files = readdir(dir)

    traitsource = CSV.read(string(dir, Files[1]), DataFrame)
    traits = names(select(traitsource, Not([:Species, :Individual])))

    df_dict = Dict{String, DataFrame}()
# Collates all the simulated data for each trait, can do stuff to the data frames in the dictionary
for i in traits 
    #file_idx = findall(x -> i == split(x, "_")[1], Files)
    trait_data = DataFrame()
    iter = 1
    for file in Files
        data = CSV.read(string(dir, file), DataFrame)
        if iter == 1
            trait_data[!, :Species] = data[!, :Species]
        end
        iter = iter + 1
        trait_data[!, file] = data[!, i]
    end
    df_dict[i] = trait_data
end
return df_dict
end 

newdir = string(dir, "MI_", extensions[6])

PMMnoPVRone = Phybridge_MIDict(newdir)

femur = PMMnoPVRone["femur"]

gdf = groupby(femur, :Species)

p = plot(size = (1600, 1000))
for group in gdf
    nospec = select(group, Not([:Species]))
    iter = 1
    vec = Vector{Float64}()
    for row in eachrow(nospec)
        current = collect(row)
        if iter == 1
            vec = current
        else
            vec = [vec; current]
        end
        iter = iter + 1
    end
    U = kde(vec)
    if maximum(U.density) > 1
        println(group.Species[1])
        continue
    end
    plot!(U.x, U.density, label = group.Species[1])
end
display(p)


savefig("PMMPVR")

