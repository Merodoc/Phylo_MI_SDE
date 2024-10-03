using Phylo


mars_tree = open(parse(RootedTree), Phylo.path("C:/PhD/Phylo_MI_SDE/Workflow_Code/newtree.nwk"))
#"C:/PhD/Phylo_MI_SDE/Workflow_Code/Data/Mars_TimeTree.nwk"
#C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/Mars_TimeTree.nwk
plot(mars_tree)

using CSV
using DataFrames

# Read in Data
# Will need to do similar data matching to that in R
mars = CSV.read("C:/PhD/Phylo_MI_SDE/Workflow_Code/mars.csv", DataFrame)
mars_avg = CSV.read("C:/PhD/Phylo_MI_SDE/Workflow_Code/Imp_Mars_PVR25new5.csv", DataFrame, types = [String, Float64])
using Statistics
#C:/PhD/Phylo_MI_SDE/Workflow_Code
#C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code

leaves = Vector{String}()

function Phylo_Bridge(start, fin, fin_time, anc, dt, samples)
    N = 1:samples
    Xhat = Vector{Float64}()

    for n in N
        B = sample(0:dt:start, WienerBridge(fin_time,fin), start)
        idx = findall(x -> x == anc, B.tt)
        val = B.yy[idx][1]
        push!(Xhat, val)
    end

return Xhat
end


# for loop over all children from root, if child name in getleafnames(tree) then we see if other branch is a child
#Turning this into a function

function Leaf_Prune(tree, start_vals)
    #Currently this iterates through all the nodes in the tree from leaves to root and returns the list
leaves = getleafnames(tree)
root = first(nodenamefilter(isroot, tree))
data = DataFrame(Nodes = leaves, Vals = start_vals)
iter = 0
idx = 0
for leaf in leaves
    idx = idx + 1
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
        time = sum(bridgelen)
        #have to figure out a better way to assign values to child nodes, maybe make a val, node dict in initialization
        #bridgesim = Phylo_Bridge()
        #bridgelen works - we can now do the diffusion bridge on them
    
    push!(data, (parent,0))
    #Should be able to use bridgelen for bridge operations we wanna try

    #Currently breaks on root node
    #If Statement before the push for isroot: it breaks the loop there
end
    #if i get here, get branch lengths for the outbound branches of parent,
    #get the bridge here and take correct point as new value for parent
    #remove children from the list of names
    #add (parent, val) to the list of names


return data
end

start_vals = mars_avg[!,2]
start_vals
pruned = Leaf_Prune(mars_tree, start_vals)

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

test = findall(x -> x == species, leaves)


mars_avg2 = DataFrame(Species = leaves)
data = Vector{Float64}()
for i in leaves
    idx = findall(x -> x == i, species)
    println(idx)
    try
        push!(data, mars_avg[idx[1],2])
    catch
        println(i)
        continue
    end
end

data