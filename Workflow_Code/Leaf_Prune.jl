using Phylo


mars_tree = open(parse(RootedTree), Phylo.path("C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/Mars_TimeTree.nwk"))

plot(mars_tree)
plot(hummers, marker_z = trait)

using CSV
using DataFrames

# Read in Data
# Will need to do similar data matching to that in R
mars = CSV.read("C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/mars.csv", DataFrame)
mars_avg = CSV.read("C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/Imp_Mars_PVR25.csv", DataFrame)
using Statistics



leaves = Vector{String}()


# for loop over all children from root, if child name in getleafnames(tree) then we see if other branch is a child
#Turning this into a function

function Leaf_Prune(tree, start_vals)
    #Currently this iterates through all the nodes in the tree from leaves to root and returns the list
leaves = getleafnames(tree)
root = first(nodenamefilter(isroot, tree))
iter = 0
for leaf in leaves
    if isroot(mars_tree, leaf)
        return leaves
    end
    #should be a dictionary of leaves and their values
    parent = getparent(tree, leaf)
    children = Vector{String}()
    #print(iter)
    iter = iter +1
    #might be good time for a try/catch
    try
        #test that all the children of the parent leaf are leaves
        for ch in getchildren(tree, parent)
            name = getnodename(tree, ch)
            if name in getleafnames(tree)
                push!(children, name)
                #print(name)
                continue
            else 
                break
            end
        end
    catch
        print("CAUGHT")
        continue
    end
    # not all children are leaves, so return to testing leaves until we find one that does
    #if all children of parent leaf are leaves then do stuff
    bridgelen = Vector{Float64}()
    #print(getnodename(mars_tree, parent))
    #optimizing this depends on what inputs the bridge needs
        for branch in getoutbounds(tree, parent)
            len = getlength(tree, branch)
            push!(bridgelen, len)
        end
        print(bridgelen)

        #bridgelen works - we can now do the diffusion bridge on them
    
    push!(leaves, getnodename(tree, parent))
    #Should be able to use bridgelen for bridge operations we wanna try

    #Currently breaks on root node
    #If Statement before the push for isroot: it breaks the loop there
end
    #if i get here, get branch lengths for the outbound branches of parent,
    #get the bridge here and take correct point as new value for parent
    #remove children from the list of names
    #add (parent, val) to the list of names


return leaves
end

start_vals = mars_avg[!,2]
start_vals
pruned = Leaf_Prune(mars_tree)

for i in 1:length(getleafnames(mars_tree))
    idx = findall(x -> x == getleafnames(mars_tree)[i], mars_avg[:,1])
    #idx = parse.(Int, idx)
    if length(idx) >= 1
        idx = idx[1]
    end
    print(idx)
    val = mars_avg[idx,2]
    #print(val)
end

mars_avg
findall(x -> x == getleafnames(mars_tree), mars_avg[!,1])

species = mars_avg[!,1]
leaves = getleafnames(mars_tree)

findall(x -> x == species, leaves)