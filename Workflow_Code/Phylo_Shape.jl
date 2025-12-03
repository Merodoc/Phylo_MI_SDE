using DataFrames
using Phylo
using CSV
using Plots

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


dir = "Workflow_Code/results_200225combined/"
#dir = "C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/results_180125/"

Phy_data = Phybridge_Dict(dir)

using StatsKit

femur = Phy_data["femur"]

gdf = groupby(femur, :Species)

aep_femur = select(gdf[1], Not([:Species]))

aep_femur = collect(eachrow(aep_femur)[1])

aep_femur_kde = kde(aep_femur)


using StatsKit
using PlotlyJS
using Plots
x = aep_femur_kde.x
y = aep_femur_kde.density

Plots.plot(x, y)


femur = Phy_data["femur"]
das_m = filter(:Species => ==("Dasyurus_maculatus"), femur)
das_m = select(das_m, Not([:Species]))
das_m = collect(eachrow(das_m)[1])
dasm_kde = kde(das_m)
dasm_x = dasm_kde.x 
dasm_y = dasm_kde.density

Plots.plot(dasm_x, dasm_y)

das_v = filter(:Species => ==("Dasyurus_viverrinus"), femur)
das_v = select(das_v, Not([:Species]))
das_v = collect(eachrow(das_v)[1])
dasv_kde = kde(das_v)
dasv_x = dasv_kde.x 
dasv_y = dasv_kde.density

Plots.plot!(dasv_x, dasv_y)

thy_c = filter(:Species => ==("Thylacinus_cynocephalus"), femur)
thy_c = select(thy_c, Not([:Species]))
thy_c = collect(eachrow(thy_c)[1])
thyc_kde = kde(thy_c)
thyc_x = thyc_kde.x 
thyc_y = thyc_kde.density

Plots.plot(dasm_x, dasm_y, title = "Density Approximation of Femur Sagittal Head Length", label = "Dasyurus maculatus")
Plots.plot!(dasv_x, dasv_y, label = "Dasyurus viverrinus")
Plots.plot!(thyc_x, thyc_y, label = "Thylacinus cynocephalus")
Plots.savefig("Kernel_Density_femur_1LPMMnoPVR")


U = kde(thy_c, bandwidth = 1.0)
Plots.plot(U.x, U.density)
Plots.plot!(thyc_x, thyc_y)

using Phylo
try
    mars_tree = open(parse(RootedTree), Phylo.path("C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/newtree.nwk"))
catch
    mars_tree = open(parse(RootedTree), Phylo.path("C:/Users/uqrelso1/Documents/GitHub/Phylo_MI_SDE/Workflow_Code/newtree.nwk"))
end

import Random
Random.seed!(123)

Plots.plot(mars_tree)

das_parent = "'14'"
node14 = filter(:Species => ==("'14'"), femur)
node14 = select(node14, Not([:Species]))
node14 = collect(eachrow(node14)[1])
node14_kde = kde(node14)

#Extant Dasyurus - Parent
Plots.plot(node14_kde.x, node14_kde.density, label = "Dasyurus Ancestor", title = "Density Comparison between Child nodes and Parent")
Plots.plot!(dasm_x, dasm_y, label = "Dasyurus maculatus")
Plots.plot!(dasv_x, dasv_y, label = "Dasyurus viverrinus")
Plots.savefig("dasyurusancestordensity_femur1LPMMnoPVR")

thy_parent = "'13'"
node13 = filter(:Species => ==(thy_parent), femur)
node13 = select(node13, Not([:Species]))
node13 = collect(eachrow(node13)[1])
node13_kde = kde(node13)

Plots.plot(node13_kde.x, node13_kde.density, label = "Thylacine ancestor")
Plots.plot!(node14_kde.x, node14_kde.density, label = "Dasyurus ancestor")
Plots.plot!(thyc_x, thyc_y, label = "Thylacinus cynocephalus")

Plots.savefig("thylacancestorfemur1LPMMnoPVR")

root = "Node 69"

root_val = filter(:Species => ==(root), femur)
root_val = select(root_val, Not([:Species]))
root_val = collect(eachrow(root_val)[1])
root_kde = kde(root_val)

Plots.plot(root_kde.x, root_kde.density)

using Statistics

std_dict = Dict{String, Float64}()
mean_dict = Dict{String, Float64}()
species_list = Vector{String}()
std_list = Vector{Float64}()
means = Vector{Float64}()
p = Plots.plot()
for species in reverse(getnodenames(mars_tree))
    println(species)
    std_dict2 = Dict{String, Float64}()
    row = filter(:Species => ==(species), femur)
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

display(p)


std_dict["Node 69"]
mean_dict["Node 69"]

Plots.savefig("Tree_Femur_STDev")

Plots.plot(mars_tree, title = "Mean sampled femur trochantericfossa length", size = (1400, 800), linewidth = 2, marker_z = means, markersize = 10, linecolor = :purple)
Plots.savefig("1LPMMnoPVR_TreeMeansbyvar")
Plots.plot(mars_tree, title = "Standard Deviation femur trochentericfossa length", size = (1400, 800), linewidth = 5, marker_z = std_dict, markersize = 10)

Plots.savefig("Tree_femur_sdev_1LPMMnoPVR")


Plots.plot(mars_tree, title = "Standard Deviation femur trochentericfossa length", size = (1400, 800), linewidth = 2, marker_z = means, markersize = means.*std_list)

species_list

std_dict

using Bridge
function Phylo_Bridge_Plot(start, final, dt, t, n)
    bridge_data = DataFrame()
    for i in 1:n
        startval = sample(start)
        finalval = sample(final)
        B = sample(0:dt:t, WienerBridge(t, finalval), startval)
        bridge_data[!, string(i)] = Vector(B.yy)
    end
    return(bridge_data)
end

test = Phylo_Bridge_Plot(das_v, das_m, 0.01, 12.52, 50)

np = length(eachrow(test)[1])
p = Plots.plot(np, title = "Bridge Density between D.maculatus and D.viverrinus", xlabel = "Tree Depth", ylabel = "Femur Trochantericfossa Length", legend = false)
iter = 1
maxy = 0
miny = 50
for i in 1:np
    y = collect(eachcol(test)[i])
    #println(y)
    #println(y)
    x = 0:0.01:12.52
    Plots.plot!(x, y, marker_z = (6.26, y[626]), marker_size = 10)
end
x = [6.26, 6.26]
y = [5, 22]
Plots.plot!(x, y, label = "Predicted Ancestor Time", linewidth = 5, linecolor = :red)
display(p)

Plots.savefig("BridgeFigure1LPMMnoPVR")




y = collect(eachcol(test)[1])
x = 0:0.01:12.52
Plots.plot(x, y)
y = collect(eachcol(test)[2])


using PlotlyJS

df = dataset(DataFrame, "tips")

PlotlyJS.plot(df, y=:total_bill, kind = "box")

femur



femur = select(permutedims(femur, 1), Not([:Species]))
p1 = PlotlyJS.plot(femur, y=:Thylacinus_cynocephalus, kind = "box")
p2 = PlotlyJS.plot(femur, y=:Dasyurus_maculatus, kind = "box")

p = [p1 p2]

p1 = PlotlyJS.plot(femur, y=:Thylacinus_cynocephalus, kind = "scatter")

x = das_m
y = das_v

PlotlyJS.plot(histogram2dcontour(x=x,y=y))

PlotlyJS.plot(femur, x=:Dasyurus_viverrinus, y=:Dasyurus_maculatus, xbingroup = "x", ybingroup = "y", kind ="histogram2d")

std(das_m)

getnodedata(mars_tree, "Node 69")

newdf = mapcols(col -> std(col), femur)
names(femur)

newdf2 = DataFrame(Species = names(femur))
newdf2[:, :Mean] = collect(mapcols(col -> mean(col), femur)[1,:])
newdf2
newdf2[:, :Stdev] = collect(mapcols(col -> std(col), femur)[1,:])

heights = nodeheights(mars_tree)

heightdf = DataFrame(Species = heights.axes[1][:], Depth = collect(heights))

test_df = innerjoin(heightdf, newdf2, on = :Species)

y = test_df[!,4]
y2 = test_df[!,4]
x = test_df[!,2]
Plots.plot(x,y./y2, seriestype=:scatter)

plotlyjs()
gr()
p = Plots.plot()
scatter!(x, y)

PlotlyJS.plot(test_df, x=:Depth, y=:Mean, xbingroup="x", ybingroup="y", kind="histogram2d")
#Get Bounds for each node
#Split the depths then get KDEs for the depths of the tree? 


femur2 = Phy_data["femur"]

femur_df = innerjoin(heightdf, femur2, on =:Species)

gdf = groupby(femur_df, :Depth)

gdf[1]

root = select(gdf[1], Not([:Species, :Depth]))
root2 = select(gdf[1], Not([:Depth]))
rootvals = Matrix(root)
rootvals = vec(rootvals)
rootkde = kde(rootvals)


#Plot of the density at leaves vs root
x = mean.(eachrow(root))
y = zeros(Float64, 21, 1)
Plots.scatter(x, y, label = "Mean Leaf Value", title = "Leaf KDE compared to Root")
Plots.plot!(rootkde.x, rootkde.density, label = "Density at Leaves")

gdf[34]

node = select(gdf[34], Not([:Species, :Depth]))
nodevals = vec(Matrix(node))
nodekde = kde(nodevals)

Plots.plot!(nodekde.x, nodekde.density, label = "Density at Root")
Plots.savefig("RootLeafKDECombined")
Plots.scatter(x = mean.(eachrow(root)), y = zeros(Float64, 21, 1))
maximum(nodevals)
minimum(nodevals)

root2
root2 = select(gdf[1], Not([:Depth]))
root2 = select(permutedims(root2, 1), Not([:Species]))


p = Plots.plot(nodekde.x, nodekde.density, label = "Root", linewidth = 4, xlimits = (0, 50), size = (600,800), title = "Root density compared to Leaves")

for i in 1:length(eachcol(root2))
    x = root2[!, i]
    xkde = kde(x, bandwidth = 1)
    Plots.plot!(xkde.x, xkde.density, label = names(root2)[i])
end

display(p)

Plots.savefig("Root_Leaf_KDE")


heights = nodeheights(mars_tree)

heightdf = DataFrame(Species = heights.axes[1][:], Depth = collect(heights))

function MICompare(df, heightdf)
com_df = innerjoin(heightdf, df, on =:Species)

gdf = groupby(com_df, :Depth)

root = select(gdf[1], Not([:Species, :Depth]))
rootvals = Matrix(root)
rootvals = vec(rootvals)
rootkde = kde(rootvals)
Plots.plot(rootkde.x, rootkde.density, label = "Density at Leaves")
end

for i in values(Phy_data)
    MICompare(df)
end
