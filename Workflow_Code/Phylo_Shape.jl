
dir = "C:/PhD/Phylo_MI_SDE/Workflow_Code/Sampled_Results_211024/"

Files = readdir(dir)


traits = Vector{String}()
for file in Files
    trait = split(file, "_")[1]
    if trait ∉ traits
        push!(traits, trait)
    end
end



traits

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


dir = "C:/PhD/Phylo_MI_SDE/Workflow_Code/Sampled_231024/"

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

dentary = Phy_data["femur"]
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
savefig("Kernel_Density_femur")


U = kde(thy_c, bandwidth = 1.0)
plot(U.x, U.density)
plot!(thyc_x, thyc_y)

using Phylo

mars_tree = open(parse(RootedTree), Phylo.path("C:/PhD/Phylo_MI_SDE/Workflow_Code/newtree.nwk"))
import Random
Random.seed!(123)

plot(mars_tree)

das_parent = "'14'"
node14 = filter(:Species => ==("'14'"), dentary)
node14 = select(node14, Not([:Species]))
node14 = collect(eachrow(node14)[1])
node14_kde = kde(node14)

#Extant Dasyurus - Parent
plot(node14_kde.x, node14_kde.density)
plot!(dasm_x, dasm_y)
plot!(dasv_x, dasv_y)

thy_parent = "'13'"
node13 = filter(:Species => ==(thy_parent), dentary)
node13 = select(node13, Not([:Species]))
node13 = collect(eachrow(node13)[1])
node13_kde = kde(node13)

plot(node13_kde.x, node13_kde.density)
plot!(node14_kde.x, node14_kde.density)
plot!(thyc_x, thyc_y)



root = "Node 69"

root_val = filter(:Species => ==(root), dentary)
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
for species in getnodenames(mars_tree)
    println(species)
    std_dict2 = Dict{String, Float64}()
    row = filter(:Species => ==(species), dentary)
    vals = select(row, Not([:Species]))
    vals = collect(eachrow(vals)[1])
    x = std(vals)
    xhat = mean(vals)
    xhat = round(xhat, digits = 3)
    std_dict[species] = x
    mean_dict[species] = xhat
    push!(species_list, species)
    push!(means, xhat)
    if x < 1
        x = 1
    end
    push!(std_list, round(x/xhat, digits = 2))
end

display(p)

plot(mars_tree, showtips = false, size = (800, 600), linewidth = 5, line_z = std_dict, series_annotations = text.mean_dict)

std_dict["Node 69"]
mean_dict["Node 69"]

savefig("Tree_Femur_STDev")

plot(mars_tree, showtips = false, size = (800, 600), linewidth = 2, marker_z = reverse(means), markersize = reverse(20*std_list))

savefig("Tree_Dentary_Means_byvar2")

species_list

getnodenames(mars_tree)

test = Phybridge_Dict(dir)

test["MarsMI"]