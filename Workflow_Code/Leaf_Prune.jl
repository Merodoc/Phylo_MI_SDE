using Phylo
using CSV
using DataFrames

# Read in Data
# Will need to do similar data matching to that in R
mars = CSV.read("C:/PhD/Phylo_MI_SDE/Workflow_Code/mars.csv", DataFrame)
mars_avg = CSV.read("C:/PhD/Phylo_MI_SDE/Workflow_Code/Imp_Mars_PVR25new5.csv", DataFrame, types = [String, Float64])
using Statistics
using Bridge
#C:/PhD/Phylo_MI_SDE/Workflow_Code
#C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code



function Phylo_Bridge(start, fin, fin_time, anc, dt, samples)
    #start = start value
    #fin = final value
    #fin_time = total time
    #anc = time point of the ancestor node
    #samples = number of repeats
    N = 1:samples
    Xhat = Vector{Float64}()

    for n in N
        B = sample(0:dt:fin_time, WienerBridge(fin_time,fin), start)
        idx = findall(x -> x == anc, B.tt)
        val = B.yy[idx][1]
        push!(Xhat, val)
    end

return Xhat
end



# for loop over all children from root, if child name in getleafnames(tree) then we see if other branch is a child
#Turning this into a function

function Leaf_Prune(tree, start_vals, dt = 0.01, samples = 100)
    #Currently this iterates through all the nodes in the tree from leaves to root and returns the list
leaves = getleafnames(tree)
root = first(nodenamefilter(isroot, tree))
data = DataFrame(Nodes = leaves, Vals = start_vals)
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
        row1 = filter(row -> row.Nodes == getnodename(tree, ch[1]), data)
        row2 = filter(row -> row.Nodes == getnodename(tree, ch[2]), data)
        val1 = row1.Vals[1]
        val2 = row2.Vals[1]
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


pruned = Leaf_Prune(mars_tree, start_vals)
plot(mars_tree, showtips = false, marker_z = plot_dict, linewidth = 5, markersize = 15)


plot_dict = Dict()

for i in eachrow(pruned)
    push!(plot_dict, i[1] => i[2])
end


function Leaf_Prune2(tree, df, variable, species = "Species", dt = 0.01, samples = 100)
    #Currently this iterates through all the nodes in the tree from leaves to root and returns the list
root = first(nodenamefilter(isroot, tree))
data = DataFrame()
sp_names = Vector{String}()

for i in df[!, species]
    push!(sp_names, i)
end

sp_trait = Vector{Float64}()

for i in df[!, variable]
    push!(sp_trait, i)
end


data[!, :Species] = sp_names
data[!, :Variable] = sp_trait
v = []
children = ["NA", "NA"]
push!(v, children)
data[!, :Children] .= v
data[!,:Child1] .= 0.
data[!,:Child1Val] .= 0.
data[!, :Child2] .= 0.
data[!,:Child2Val] .= 0.
data[!,:BridgeMean] .= 0.
data[!, :BridgeVar] .= 0.
iter = 0
idx = 0
vals = start_vals
leaves = getleafnames(mars_tree)
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
            trait_sample = sample(bridgesim)
            Xhat = mean(bridgesim)
            Var = var(bridgesim)

        #bridgelen works - we can now do the diffusion bridge on them
    push!(data, (parent, trait_sample, children, bridgelen[1], val1, bridgelen[2], val2, Xhat, Var))
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





