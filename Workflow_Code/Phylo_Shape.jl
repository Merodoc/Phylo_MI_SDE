
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

Files = readdir(dir)


traits = Vector{String}()
for file in Files
    trait = split(file, "_")[1]
    if trait ∉ traits
        push!(traits, trait)
    end
end