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
mars_tree = open(parse(RootedTree), Phylo.path("C:/Users/uqrelso1/Documents/GitHub/Phylo_MI_SDE/Workflow_Code/newtree.nwk"))
#dir = "C:/PhD/Phylo_MI_SDE/Workflow_Code/Sampled2l_241024/"
dir = "C:/Users/uqrelso1/Documents/GitHub/Phylo_MI_SDE/Workflow_Code/results_200225combined/"


function TreeSTDPlot(tree, df, Title)
    std_dict = Dict{String, Float64}()
    mean_dict = Dict{String, Float64}()
    species_list = Vector{String}()
    std_list = Vector{Float64}()
    means = Vector{Float64}()
    p = Plots.plot(title = Title)
    for species in reverse(getnodenames(tree))
        println(species)
        row = filter(:Species => ==(species), df)
        vals = select(row, Not([:Species]))
        vals = collect(eachrow(vals)[1])
        x = std(vals)
        xhat = mean(vals)
        xhat = round(xhat, digits = 3)
        std_dict[species] = x/xhat
        mean_dict[species] = xhat
        push!(species_list, species)
        push!(means, xhat)
        push!(std_list, round(x/xhat, digits = 2))
    end
    Plots.plot!(mars_tree, size = (1400, 800), linewidth = 5, marker_z = std_dict, markersize = 10)
    display(p)
end


function TreeMeanPlot(tree, df, Title)
    mean_dict = Dict{String, Float64}()
    species_list = Vector{String}()
    means = Vector{Float64}()
    p = Plots.plot(title = Title)
    for species in reverse(getnodenames(tree))
        println(species)
        row = filter(:Species => ==(species), df)
        vals = select(row, Not([:Species]))
        vals = collect(eachrow(vals)[1])
        xhat = mean(vals)
        xhat = round(xhat, digits = 3)
        mean_dict[species] = xhat
        push!(species_list, species)
        push!(means, xhat)
    end
    Plots.plot!(mars_tree, size = (1400, 800), linewidth = 5, marker_z = mean_dict, markersize = 10)
    display(p)
end

function PhyGetMeans(df)
    newdf = select(df, Not([:Species]))
    meanlist = Vector{Float64}()
    for row in eachrow(newdf)
        xhat = mean(row)
        push!(meanlist, xhat)
    end
    return meanlist
end

function PhyGetSTD(df)
    newdf = select(df, Not([:Species]))
    stdlist = Vector{Float64}()
    for row in eachrow(newdf)
        σ = std(row)
        push!(stdlist, σ)
    end
    return stdlist
end






Folder = readdir(dir)

string(dir, readdir(dir)[1], "/")

heights = nodeheights(mars_tree)

heightdf = DataFrame(Species = heights.axes[1][:], Depth = collect(heights))

p = Plots.plot(title = "Leaf Density comparison between MI strategies")
for file in Folder
    newdir = string(dir, file, "/")
    Phy_data = Phybridge_Dict(newdir)
    femur = Phy_data["femur"]
    MICompare(femur, heightdf, file)
end

display(p)
Plots.savefig("MILeafComparison")

p = Plots.plot(title = "Root Density comparison between MI strategies")
for file in Folder
    newdir = string(dir, file, "/")
    Phy_data = Phybridge_Dict(newdir)
    femur = Phy_data["femur"]
    RootCompare(femur, heightdf, file)
end

display(p)
Plots.savefig("MIRootComparisons")

for file in Folder
    newdir = string(dir, file, "/")
    Phy_data = Phybridge_Dict(newdir)
    femur = Phy_data["femur"]

    com_df = innerjoin(heightdf, femur, on =:Species)
    TreeSTDPlot(mars_tree, femur, string("Standard Deviation for tree: ", file)) 
    TreeMeanPlot(mars_tree, femur, string("Means for tree: ", file))
    means = PhyGetMeans(femur)
    stds = PhyGetSTD(femur)
    x = com_df[!, :Depth]
    p = Plots.scatter(x, stds, title = string("Standard Deviation vs Depth for: ", file))
    display(p)
    p2 = Plots.scatter(x, means, title = string("Mean vs Depth for: ", file), yerror = stds)   
    display(p2)
end

testfile = string(dir, Folder[1], "/")

Phy_data = Phybridge_Dict(testfile)
femur = Phy_data["femur"]
com_df = innerjoin(heightdf, femur, on =:Species)

femur

x = com_df[!, :Depth]
newdata = select(com_df, Not([:Species, :Depth]))
newdata

stdlist = Vector{Float64}()
meanlist = Vector{Float64}()



depthdf = select(com_df, Not([:Species]))
gdf = groupby(depthdf, :Depth)

depthmeans = Vector{Float64}()
depthstds = Vector{Float64}()
depths = Vector{Float64}()

for group in gdf
    push!(depths, mean(group[!, :Depth]))
    mat = Matrix(select(group, Not([:Depth])))
    xhat = mean(mat)
    σ = std(mat)
    push!(depthmeans, xhat)
    push!(depthstds, σ)
end

LMtest = DataFrame(Depth = x, Means = meanlist, STDs = stdlist)

using GLM

ols = lm(@formula(STDs ~ Depth), LMtest)

Plots.scatter(depths, depthmeans, smooth = true, ribbon = depthstds)

df = dataset(DataFrame, "tips")

PlotlyJS.plot(com_df, x=:Depth, y=:"11326", xbingroup="x", ybingroup="y", kind="histogram2d")