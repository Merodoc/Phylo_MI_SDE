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

function Leaf_Prune2(tree, df, variable, species = "Species", dt = 0.01, samples = 100)
    #Currently this iterates through all the nodes in the tree from leaves to root and returns the list
root = first(nodenamefilter(isroot, tree))
data = DataFrame()
data[!, :Species] = df[!, species]
data[!, :Variable] = df[!, variable]
v = []
children = ["NA", "NA"]
push!(v, children)
data[!, :Children] .= v
data[!,:Child1] .= 0.
data[!,:Child1Val] .= 0.
data[!, :Child2] .= 0.
data[!,:Child2Val] .= 0.
data[!,:BridgeVar] .= 0.
iter = 0
idx = 0
vals = start_vals
for leaf in leaves
    idx = idx + 1
    if isroot(mars_tree, leaf)
        #bridgelen = Vector{Float64}()
        #print(getnodename(mars_tree, parent))
        #optimizing this depends on what inputs the bridge needs
                #for branch in getoutbounds(tree, parent)
                #    len = getlength(tree, branch)
                #    push!(bridgelen, len)
                #end
                #time = sum(bridgelen)
                #push!(data, (leaf,0, bridgelen[1], bridgelen[2]))
        return data
        break
    end
    #should be a dictionary of leaves and their values
    parent = getparent(tree, leaf)
    if parent ∈ leaves 
        continue
        #println("PARENT IN LEAVES")
    end
    #print(iter)
    iter = iter +1
    #might be good time for a try/catch
        #test that all the children of the parent leaf are leaves
    ch = getchildren(tree, parent)
    children = [getnodename(tree, ch[1]), getnodename(tree,ch[2])]
    if getnodename(tree, ch[1]) ∈ leaves && getnodename(tree, ch[2]) ∈ leaves

        #Pull existing trait values from the data frame
        row1 = filter(row -> row.Species == getnodename(tree, ch[1]), data)
        row2 = filter(row -> row.Species == getnodename(tree, ch[2]), data)
        val1 = row1.Variable[1]
        val2 = row2.Variable[1]
        bridgelen = Vector{Float64}()
    #print(getnodename(mars_tree, parent))
    #optimizing this depends on what inputs the bridge needs

            for branch in getoutbounds(tree, parent)
                len = getlength(tree, branch)
                push!(bridgelen, len)
            end

            time = sum(bridgelen)
            # Need some Error work in here to guarantee that timescale is divisible by dt
            time = round(time, digits = 2)
            anc_time = round(bridgelen[1], digits = 2)
            bridgesim = Phylo_Bridge(val1, val2, time, anc_time, dt, samples)
            Xhat = mean(bridgesim)
            Var = var(bridgesim)

        #bridgelen works - we can now do the diffusion bridge on them
    push!(data, (parent,Xhat, children, bridgelen[1], val1, bridgelen[2], val2, Var))
    push!(leaves, parent)
    #Should be able to use bridgelen for bridge operations we wanna try

            
    #Currently breaks on root node
    #If Statement before the push for isroot: it breaks the loop there
end
    #if i get here, get branch lengths for the outbound branches of parent,
    #get the bridge here and take correct point as new value for parent
    #remove children from the list of names
    #add (parent, val) to the list of names

end
return data
end
leaves = getleafnames(mars_tree)
Base.load_InteractiveUtils
sort!(sampled_df, Species = leaves)

Leaf_Prune2(mars_tree, sampled_df, "dentary")