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
leaves = getleafnames(mars_tree)
iter = 0
for leaf in leaves
    #should be a dictionary of leaves and their values
    parent = getparent(mars_tree, leaf)
    children = Vector{String}()
    #print(iter)
    iter = iter +1
    #might be good time for a try/catch
    try
        #test that all the children of the parent leaf are leaves
        for ch in getchildren(mars_tree, parent)
            name = getnodename(mars_tree, ch)
            if name in getleafnames(mars_tree)
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
        for branch in getoutbounds(mars_tree, parent)
            len = getlength(mars_tree, branch)
            push!(bridgelen, len)
        end
        #print(bridgelen)
        #print(children)
        #leaves <- filter!(e->e∉children, leaves)
    
    push!(leaves, getnodename(mars_tree, parent))
    #Should be able to use bridgelen for bridge operations we wanna try

    #Currently breaks on root node
    #If Statement before the push for isroot: it breaks the loop there
end
    #if i get here, get branch lengths for the outbound branches of parent,
    #get the bridge here and take correct point as new value for parent
    #remove children from the list of names
    #add (parent, val) to the list of names


leaves

#Turning this into a function

function = Leaf_Prune(tree)