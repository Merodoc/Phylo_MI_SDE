using DataFrames
using Phylo
using CSV
using Plots
using StatsKit
using PlotlyJS
using Statistics


function Phybridge_Dict(dir)
    try 
        readdir(dir)
    catch
        println("Invalid directory")
    end
    Files = readdir(dir)
    traits = Vector{String}()

    for file in Files
        trait = split(file, "_")[1]
        if trait ∉ traits
            push!(traits, trait)
        end
    end

    df_dict = Dict{String, DataFrame}()
# Collates all the simulated data for each trait, can do stuff to the data frames in the dictionary
for i in traits 
    file_idx = findall(x -> i == split(x, "_")[1], Files)
    trait_data = DataFrame()
    iter = 1
    for j in file_idx
        data = CSV.read(string(dir, Files[j]), DataFrame)
        titles = names(data)
        for k in 1:length(eachcol(data))
            if iter == 1
                trait_data[!, titles[k]] = eachcol(data)[k]
            else
                if k != 1
                    trait_data[!, string(titles[k], j)] = eachcol(data)[k]
                end
            end
            iter = iter + 1 
        end

    end
    df_dict[i] = trait_data
end
return df_dict
end 



function MICompare(df, heightdf, label)
com_df = innerjoin(heightdf, df, on =:Species)

gdf = groupby(com_df, :Depth)

root = select(gdf[1], Not([:Species, :Depth]))
rootvals = Matrix(root)
rootvals = vec(rootvals)
rootkde = kde(rootvals)
Plots.plot!(rootkde.x, rootkde.density, label = label)
end

function RootCompare(df, heightdf, label)
    com_df = innerjoin(heightdf, df, on =:Species)
    
    gdf = groupby(com_df, :Depth)
    idx = length(gdf)
    root = select(gdf[idx], Not([:Species, :Depth]))
    rootvals = Matrix(root)
    rootvals = vec(rootvals)
    rootkde = kde(rootvals)
    Plots.plot!(rootkde.x, rootkde.density, label = label)
    end
mars_tree = open(parse(RootedTree), Phylo.path("C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/newtree.nwk"))
#dir = "C:/PhD/Phylo_MI_SDE/Workflow_Code/Sampled2l_241024/"
dir = "C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/results_180125/"

Folder = readdir(dir)

string(dir, readdir(dir)[1], "/")

heights = nodeheights(mars_tree)

heightdf = DataFrame(Species = heights.axes[1][:], Depth = collect(heights))

p = Plots.plot(title = "Leaf Density comparison between MI strategies")
for file in Folder
    newdir = string(dir, file)
    Phy_data = Phybridge_Dict(dir)
    femur = Phy_data["femur"]
    MICompare(femur, heightdf, file)
end

display(p)
Plots.savefig("MILeafComparison")

p = Plots.plot(title = "Root Density comparison between MI strategies")
for file in Folder
    newdir = string(dir, file)
    Phy_data = Phybridge_Dict(dir)
    femur = Phy_data["femur"]
    RootCompare(femur, heightdf, file)
end

display(p)
Plots.savefig("MIRootComparisons")
    
