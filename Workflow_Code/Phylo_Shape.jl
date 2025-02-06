using DataFrames
using Phylo
using CSV

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


#dir = "C:/PhD/Phylo_MI_SDE/Workflow_Code/Sampled2l_241024/"
dir = "C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/MI_1lPMMnoPVR_121124/"

Phy_data = Phybridge_Dict(dir)

using StatsKit

femur = Phy_data["femur"]

gdf = groupby(femur, :Species)

aep_femur = select(gdf[1], Not([:Species]))

aep_femur = collect(eachrow(aep_femur)[1])

aep_femur_kde = kde(aep_femur)


test = Phybridge_Dict(dir)

Phy_data = test

using StatsKit

femur = Phy_data["femur"]

gdf = groupby(femur, :Species)

aep_femur = select(gdf[1], Not([:Species]))

aep_femur = collect(eachrow(aep_femur)[1])

aep_femur_kde = kde(aep_femur)

using Plots

x = aep_femur_kde.x
y = aep_femur_kde.density

plot(x, y)

aep_femur

femur = Phy_data["femur"]
das_m = filter(:Species => ==("Dasyurus_maculatus"), femur)
das_m = select(das_m, Not([:Species]))
das_m = collect(eachrow(das_m)[1])
dasm_kde = kde(das_m)
dasm_x = dasm_kde.x 
dasm_y = dasm_kde.density

plot(dasm_x, dasm_y)

das_v = filter(:Species => ==("Dasyurus_viverrinus"), femur)
das_v = select(das_v, Not([:Species]))
das_v = collect(eachrow(das_v)[1])
dasv_kde = kde(das_v)
dasv_x = dasv_kde.x 
dasv_y = dasv_kde.density

plot!(dasv_x, dasv_y)

thy_c = filter(:Species => ==("Thylacinus_cynocephalus"), femur)
thy_c = select(thy_c, Not([:Species]))
thy_c = collect(eachrow(thy_c)[1])
thyc_kde = kde(thy_c)
thyc_x = thyc_kde.x 
thyc_y = thyc_kde.density

plot(dasm_x, dasm_y, title = "Density Approximation of Femur Sagittal Head Length", label = "Dasyurus maculatus")
plot!(dasv_x, dasv_y, label = "Dasyurus viverrinus")
plot!(thyc_x, thyc_y, label = "Thylacinus cynocephalus")
savefig("Kernel_Density_femur_2lonlypmm")


U = kde(thy_c, bandwidth = 1.0)
plot(U.x, U.density)
plot!(thyc_x, thyc_y)

using Phylo

mars_tree = open(parse(RootedTree), Phylo.path("C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/newtree.nwk"))
import Random
Random.seed!(123)

plot(mars_tree)

das_parent = "'14'"
node14 = filter(:Species => ==("'14'"), femur)
node14 = select(node14, Not([:Species]))
node14 = collect(eachrow(node14)[1])
node14_kde = kde(node14)

#Extant Dasyurus - Parent
plot(node14_kde.x, node14_kde.density, label = "Dasyurus Ancestor", title = "Density Comparison between Child nodes and Parent")
plot!(dasm_x, dasm_y, label = "Dasyurus maculatus")
plot!(dasv_x, dasv_y, label = "Dasyurus viverrinus")
savefig("dasyurusancestordensity_femur")

thy_parent = "'13'"
node13 = filter(:Species => ==(thy_parent), femur)
node13 = select(node13, Not([:Species]))
node13 = collect(eachrow(node13)[1])
node13_kde = kde(node13)

plot(node13_kde.x, node13_kde.density, label = "Thylacine ancestor")
plot!(node14_kde.x, node14_kde.density, label = "Dasyurus ancestor")
plot!(thyc_x, thyc_y, label = "Thylacinus cynocephalus")

savefig("thylacancestorfemur")

root = "Node 69"

root_val = filter(:Species => ==(root), femur)
root_val = select(root_val, Not([:Species]))
root_val = collect(eachrow(root_val)[1])
root_kde = kde(root_val)

plot(root_kde.x, root_kde.density)

using Statistics

std_dict = Dict{String, Float64}()
mean_dict = Dict{String, Float64}()
species_list = Vector{String}()
std_list = Vector{Float64}()
means = Vector{Float64}()
p = plot()
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

savefig("Tree_Femur_STDev")

plot(mars_tree, title = "Mean sampled femur trochantericfossa length", size = (1400, 800), linewidth = 2, marker_z = means, markersize = 35*std_list, linecolor = :purple)
savefig("2lonlyFemur_TreeMeansbyvar")
plot(mars_tree, title = "Standard Deviation femur trochentericfossa length", size = (1400, 800), linewidth = 5, line_z = std_dict)

savefig("Tree_femur_sdev_2lonly")


plot(mars_tree, title = "Standard Deviation femur trochentericfossa length", size = (1400, 800), linewidth = 2, marker_z = means, markersize = means.*std_list)

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
p = plot(np, title = "Bridge Density between D.maculatus and D.viverrinus", xlabel = "Tree Depth", ylabel = "Femur Trochantericfossa Length", legend = false)
iter = 1
maxy = 0
miny = 50
for i in 1:np
    y = collect(eachcol(test)[i])
    #println(y)
    #println(y)
    x = 0:0.01:12.52
    plot!(x, y, marker_z = (6.26, y[626]), marker_size = 10)
end
x = [6.26, 6.26]
y = [5, 22]
plot!(x, y, label = "Predicted Ancestor Time", linewidth = 5, linecolor = :red)
display(p)

savefig("BridgeFigure")

for i in getleaves(mars_tree)                                                                                       
    name = i.name                                                                                                       
    name = split(name, "_")                                                                                             
    name = string(name[1][1], ".", name[2])                                                                             
    
    i.name = name
end


y = collect(eachcol(test)[1])
x = 0:0.01:12.52
plot(x, y)
y = collect(eachcol(test)[2])


MI_Dir = "C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/MI_Data2lonly/"
MIdict = Phybridge_Dict(MI_Dir)

MIdata = MIdict["MarsMI"]

MIdata.femur1


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
            trait_data[!, :Species] = data[!, 1]
        end
        iter = iter + 1
        trait_data[!, file] = data[!, i]
    end
    df_dict[i] = trait_data
end
return df_dict
end 

MI_dict = Phybridge_MIDict(MI_Dir)

MI_femur = MI_dict["femur"]

length(names(MI_femur))
p = plot(legend = false)
for i in 2:length(names(MI_femur))
    vals = collect(eachcol(MI_femur)[i])
    U = kde(vals)
    plot!(U.x, U.density)
end

display(p)

osp = filter(:Species => ==("Osphranter_robustus"), femur)
osp = select(osp, Not([:Species]))
osp = collect(eachrow(osp)[1])

U = kde(osp)

plot!(U.x, U.density, linewidth = 5, linecolor = :red)

osp = filter(:Species => ==("Onychogalea_unguifera"), femur)
osp = select(osp, Not([:Species]))
osp = collect(eachrow(osp)[1])

U = kde(osp)

plot(U.x, U.density, linewidth = 5, linecolor = :blue)
savefig("Onychogalea_femur")

std(osp)

