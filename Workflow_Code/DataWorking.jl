include("Leaf_Prune.jl")


mars_tree = open(parse(RootedTree), Phylo.path("C:/PhD/Phylo_MI_SDE/Workflow_Code/newtree.nwk"))
import Random
Random.seed!(123)

function PhyMIR_Analyze_dir(start_dir, end_dir, tree)
    #start_dir = Directory containing Phylogenetic MI data
    #end_dir = Location for final data to end up
    #tree = a Phylo.jl tree object containing a phylogenetic tree of species in data
    Files = readdir(start_dir)
    try
        mkdir(end_dir)
    catch
        return println("Invalid Return Directory")
    end
    for file in Files
        filename = string(start_dir, file)
        data = CSV.read(filename, DataFrame)
        means = select!(data, Not([:Species]))
        iter = 1
        varstring = split(file, ".")[1]
        try
            Pooledmeans = DataFrame()
            PooledVar = DataFrame()
            for i in eachcol(means)
            start_vals = i
            prune = Leaf_Prune(tree, start_vals)
            newfile = string(end_dir, varstring, "_", iter, ".csv")
            title = string("TraitVal ", iter)
            title2 = string("TraitVar ", iter)
            #if iter == 1
            #    Pooledmeans[!, :Species] = prune.Species
            #end
            Pooledmeans[!, string("TraitVal ", iter)] = prune.Vals
            #PooledVar[!, string("TraitVar ", iter)] = prune.BridgeVar
            CSV.write(newfile, prune)
            iter = iter + 1
            end
            pooledfilename = string(end_dir, varstring, "_pooledmean.csv")
            #pooledfilename2 = string(end_dir, varsting, "_pooledvar.csv")
            CSV.write(pooledfilename, Pooledmeans)
            #CSV.write(pooledfilename2, PooledVar)
        catch
        println(string("File: ", file , " invalid"))
        end

    end
    
end

dir = "C:/PhD/Phylo_MI_SDE/Workflow_Code/MI_Data/"
enddir = "C:/PhD/Phylo_MI_SDE/Workflow_Code/Sampled_Results_1/"
Random.seed!(123)
PhyMIR_Analyze_dir(dir, enddir, mars_tree)
mkdir(enddir)

Files = readdir(dir)

data = CSV.read(string(dir, Files[1]), DataFrame)
data = select!(data, Not([:Individual]))
gdf = groupby(data, :Species)

sampled_df = DataFrame()
titles = names(gdf[1])
for i in titles 
    sampled_df[!, i] = []
end
for species in gdf
    row = Vector()
    for col in eachcol(species)
        sample = rand(1:length(col))
        push!(row, col[sample])
    end
    push!(sampled_df, row)
end



# Leaf_Prune2 is an advancement on Leaf_Prune allowing us to import a DataFrame containing multiple variables and run the prune bridge method on them

# Read the files in from directory
data = CSV.read(string(dir, Files[1]), DataFrame)
# Remove the "Individual" factor as this is functionally useless in this sampling regime
data = select!(data, Not([:Individual]))

#Remove Bad Traits now

data = select!(data, Not([:Calcaneus],[:"Calcaneus.4"]))
# Group each species by their resulting traits 
gdf = groupby(data, :Species)
sampled_df = DataFrame()

# Iterate over the grouped data frame to sample a random measurement for each variable from each Species
titles = names(gdf[1])
Random.seed!(123)
for i in titles 
    sampled_df[!, i] = []
end

for species in gdf
    iter = 1
    row = Vector()
    for col in eachcol(species)
        if iter == 1
            push!(row, col[1])
        else 
        val = sample(col)
        push!(row, val)
        end
    iter = iter + 1
    end
    push!(sampled_df, row)
end

dentary_mi1_1 = Leaf_Prune2(mars_tree, sampled_df, "dentary")

plot_dict = Dict()

for i in eachrow(dentary_mi1_1)
    push!(plot_dict, i[1] => i[2])
end

plot(mars_tree, showtips = false, marker_z = plot_dict, linewidth = 5, markersize = 15)


# Now need to make the bridges sample from the point
# Bridge Should sample and use those samples
# Cannot pull bridge path out yet

data_test = CSV.read(string(dir, Files[1]), DataFrame)
data_test = select!(data_test, Not([:Individual]))

sampled_df = DataFrame()
titles = names(gdf[1])
for i in titles
    if "NA" ∈ data_test[!, i]
        println(string("Trait: ", i, " removed due to missing values"))
        data_test = select!(data_test, Not([i]))
    else
    sampled_df[!, i] = []
    end
end

gdf = groupby(data_test, :Species)


for species in gdf
    iter = 1
    row = Vector()
    for col in eachcol(species)
        if iter == 1
            push!(row, col[1])
        else 
        val = sample(col)
        push!(row, val)
        end
    iter = iter + 1
    end
    push!(sampled_df, row)
end

function Phy_Bridge_Sim(start_dir, end_dir, tree, max_iter, init_samples)
    #Read MI files from start_dir
    Files = readdir(start_dir)
    try
        #mkdir(end_dir)
    catch
        return println("Invalid Return Directory")
    end

    for file in Files
        println(file)
        data = CSV.read(string(start_dir, file), DataFrame)
        # Remove the "Individual" factor as this is functionally useless in this sampling regime
        data = select!(data, Not([:Individual]))
        iter_data = DataFrame()
        iter_data[!, :Species] = names(data)
        # Iterate over all traits
        # Create a new empty data frame that will contain the sampled trait values
        # Remove any traits that the Multiple Imputation couldn't handle
        titles = names(data)

   
        for i in 1:init_samples
            sampled_df = DataFrame()
            for title in titles
                if "NA" ∈ data[!, title]
                    #println(string("Trait: ", title, " removed due to missing values"))
                    data = select!(data, Not([title]))
                else
                    #println("In Else Loop")
                    sampled_df[!, title] = []
                end
            end
            titles = names(sampled_df)
        # Group Data by Species for sampling
            gdf = groupby(data, :Species)
        
        # Sample the Species data for each trait and add to the sampled dataframe
            for species in gdf
                iter = 1
                row = Vector()
                    for col in eachcol(species)
                        if iter == 1
                        push!(row, col[1])
                        else 
                            val = sample(col)
                            push!(row, val)
                        end
                        iter = iter + 1
                    end
                push!(sampled_df, row)
            end
        # We should now have a data frame that has a sample per species from every variable

            for trait in names(sampled_df)
                if trait == "Species"
                    continue
                else
                    trait_data = DataFrame()
                    for j in 1:max_iter
                                       
                        results = Leaf_Prune2(mars_tree, sampled_df, trait)
                        if j == 1
                            species = results[!, :Species]
                            trait_data[!, :Species] = species
                        end
                        title = string("Iter", j)
                        #filename = string(end_dir, title, "sample", i, ".csv")
                        #CSV.write(filename, results)
                        trait_vals = results[!, :Variable]
                        trait_data[!, title] = trait_vals
                    end
                end
                trait_split = split(trait, ".")
                if length(trait_split) == 2
                    trait_name = string(trait_split[1], trait_split[2])
                    filename = string(end_dir, trait_name, "_MISample", i, ".csv")
                    CSV.write(filename, trait_data) 
                else
                    filename = string(end_dir, trait, "_MISample", i, ".csv")
                    CSV.write(filename, trait_data)   

                end               
            end
        end
    end
end



dir = "C:/PhD/Phylo_MI_SDE/Workflow_Code/MI_Data/"
enddir = "C:/PhD/Phylo_MI_SDE/Workflow_Code/Sampled_Results_211024/"

Phy_Bridge_Sim(dir, enddir, mars_tree, 5, 5)


