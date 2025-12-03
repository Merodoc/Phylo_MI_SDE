include("Leaf_Prune.jl")


mars_tree = Phylo.open(parsenewick, Phylo.path("C:/Users/uqrelso1/Documents/GitHub/Phylo_MI_SDE/Workflow_Code/newtree.nwk"))
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



function Phy_Bridge_Sim(start_dir, end_dir, tree, max_iter, init_samples)
    #Read MI files from start_dir
    Files = readdir(start_dir)
    try
        mkdir(end_dir)
    catch
        return println("Invalid Return Directory")
    end
    file_number = 0
    for file in Files
        file_number = file_number + 1
        file_time = time()
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
            sample_time = time()
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
                            val = mean(col)
                            push!(row, val)
                        end
                        iter = iter + 1
                    end
                push!(sampled_df, row)
            end
        # We should now have a data frame that has a sample per species from every variable
            trait_data = DataFrame()

            for trait in names(sampled_df)
                trait_time = time()
                if trait == "Species"
                    continue
                else
                    for j in 1:max_iter              
                        results = Leaf_Prune2(mars_tree, sampled_df, trait)
                        if j == 1
                            trait_data[!, :Species] = results[!, :Species]
                        end
                        raw_data = select(results, Not([:Species]))
                        title = string(i,j)
                        trait_data[!, title] = results[!, :Variable]
                    end
                end
                trait_split = split(trait, ".")
                if length(trait_split) == 2
                    trait_name = string(trait_split[1], trait_split[2])
                    filename = string(end_dir, trait_name, "_MI", file_number, "_Sample", i, ".csv")
                    CSV.write(filename, trait_data) 
                else
                    filename = string(end_dir, trait, "_MI", file_number, "_Sample", i, ".csv")
                    CSV.write(filename, trait_data)   

                end
                elapsed_trait = time() - trait_time
                println("Time for trait ", trait, ": ", elapsed_trait, " seconds")
            end
            elapsed_sample = time()-sample_time
            println("Time for sample ", i, ": ", elapsed_sample, " seconds")
        end
        elapsed_file = time() - file_time
        println("Time for file - ", file, ": ", elapsed_file, " seconds")
    end
end


function Phy_Bridge_SimSDE(start_dir, end_dir, tree, max_iter, init_samples)
    #Read MI files from start_dir
    Files = readdir(start_dir)
    try
        mkdir(end_dir)
    catch
        return println("Invalid Return Directory")
    end
    file_number = 0
    for file in Files
        file_number = file_number + 1
        file_time = time()
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
            sample_time = time()
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
                        #println(col)
                        if iter == 1
                        push!(row, col[1])
                        else 
                            val = log(mean(col))
                            push!(row, val)
                        end
                        iter = iter + 1
                    end
                push!(sampled_df, row)
            end
        # We should now have a data frame that has a sample per species from every variable
            trait_data = DataFrame()
            #println(sampled_df[!, :Species])
            for trait in names(sampled_df)
                trait_time = time()
                if trait == "Species"
                    continue
                else
                    for j in 1:max_iter              
                        results = Leaf_PruneSDE(mars_tree, sampled_df, trait)
                        if j == 1
                            trait_data[!, :Species] = results[!, :Species]
                        end
                        raw_data = select(results, Not([:Species]))
                        title = string(i,j)
                        trait_data[!, title] = results[!, :Variable]
                    end
                end
                trait_split = split(trait, ".")
                if length(trait_split) == 2
                    trait_name = string(trait_split[1], trait_split[2])
                    filename = string(end_dir, trait_name, "_MI", file_number, "_Sample", i, ".csv")
                    CSV.write(filename, trait_data) 
                else
                    filename = string(end_dir, trait, "_MI", file_number, "_Sample", i, ".csv")
                    CSV.write(filename, trait_data)   

                end
                elapsed_trait = time() - trait_time
                println("Time for trait ", trait, ": ", elapsed_trait, " seconds")
            end
            elapsed_sample = time()-sample_time
            println("Time for sample ", i, ": ", elapsed_sample, " seconds")
        end
        elapsed_file = time() - file_time
        println("Time for file - ", file, ": ", elapsed_file, " seconds")
    end
end



#dir = "C:/Users/uqrelso1/Documents/GitHub/Phylo_MI_SDE/Workflow_Code/Data/MultipleImputes/MI_Midas/"
#enddir = "C:/Users/uqrelso1/Documents/GitHub/Phylo_MI_SDE/Workflow_Code/"
#Files = readdir(dir)


#df_dict = Dict{String, DataFrame}()
#data = CSV.read(string(dir, Files[1]), DataFrame)
#select(data, [:Species, :humerus])

 
dir2 = "C:/Users/uqrelso1/Documents/GitHub/Phylo_MI_SDE/Workflow_Code/MI_Data/"
#Files2 = readdir(dir2)
#data2 = CSV.read(string(dir2, Files2[7]), DataFrame)

#Phy_Bridge_Sim(dir, string("MI_Midas_Results_030325/"), mars_tree, 5, 1)

Phy_Bridge_SimSDE(dir2, string("MI_Midas_ResultsLog_031225/"), mars_tree, 5, 1)

