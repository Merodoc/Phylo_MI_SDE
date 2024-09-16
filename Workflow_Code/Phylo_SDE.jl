using Phylo

#Phlo.jl Practice
simpletree = parsenewick("((,Tip:1.0)Internal,)Root;")
getbranches(simpletree)

tree = open(parsenewick, Phylo.path("H1N1.newick"))

ts = open(parsenexus, Phylo.path("H1N1.trees"))

gettreeinfo(ts)

ts["TREE1"]

out = Phylo.outputtree(simpletree, Newick())

Phylo.write("test.newick", simpletree)

using Phylo, Plots
default(linecolor = :black, size = (400, 400))
hummers = open(parsenewick, Phylo.path("hummingbirds.tree"))
plot(hummers, size = (400, 600), showtips = false)

#map_depthfirst to evolve over tree

evolve(tree) = map_depthfirst((val, node) -> val + randn(), 5., tree, Float64)
trait = evolve(hummers)
plot(hummers, line_z = trait, linecolor = :RdYlBu, linewidth = 5, showtips = false)

brownsampler = BrownianTrait(hummers, "Trait")
plot(hummers, showtips = false, marker_z = rand(brownsampler), linewidth = 2, markercolor = :RdYlBu, size = (400, 600))

using StochasticDiffEq, SciMLSensitivity, Plots
using DifferentialEquations

function OU_Wrap(val, α=1.0, σ = 0.1, μ = 5.0, u0 = 5.0)

    return val + α*(μ-val) + σ*val*randn()
end

OU_Evolve(tree) = map_depthfirst((val, node) -> OU_Wrap(val), 0., tree, Float64)

trait = OU_Evolve(hummers)

function map_depthfirst2(FUN, NOISE, start, tree, eltype = nothing)
    #This is an updated Map_depthfirst function that allows us to take SDE input
    #Currently the functionality is rough
    #importantly this is now scaling the functions on time, rather than treating each node as 1 time-step
    root = first(nodenamefilter(isroot, tree))
    eltype === nothing && (eltype = typeof(FUN(start, root)))
    ret = Vector{eltype}()
    function local!(val, node)
        push!(ret, val)
        for ch in getchildren(tree, node)
            branch = getinbound(tree, ch)
            #Need to figure out how to get branches from names, or nodes from names
            len = getlength(tree, branch)
            #this section here allows u to get time between nodes for when u come to this later
            tspan = (0.0, len)

            #Iterator seems to be working, now need to fix the SDE bit
            SDE = SDEProblem(FUN, NOISE, val, tspan, noise = WienerProcess(0.0,0.0,0.0))
            sol = solve(SDE, SOSRI())
            
            #can try to make the OU_Wrap function (Val, node, time)
            local!(last(sol), ch)
        end
    end
    local!(start, root)
    ret
end




using Phylo

#Import TimeTree of Species
#Different paths for each device
#C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code for PC
#C:/PhD/Phylo_MI_SDE/Workflow_Code/ for campus
mars_tree = open(parse(RootedTree), Phylo.path("C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/Mars_TimeTree.nwk"))

plot(mars_tree)
plot(hummers, marker_z = trait)
trait = OU_Evolve
using CSV
using DataFrames

# Read in Data
# Will need to do similar data matching to that in R
mars = CSV.read("C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/mars.csv", DataFrame)
mars_avg = CSV.read("C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/Imp_Mars_PVR25.csv", DataFrame)
using Statistics


#Convert Dataframe rows to numerics and remove missing data (should not be a problem later obviously)

#filter!(row -> row.dentary != "NA", mars_avg)
dentary = mars_avg[!,2]

#Define the SDE for map_depthfirst2
W = WienerProcess(0.0,0.0,0.0)
α = 1.5
σ = 0.1
μ = mean(dentary)
μ0 = μ
f(u, p, t) =  α*(μ-u)
g(u,p, t) = σ*u

#Model Traits using the OU Model (This is purely exploratory, no data involved)
OU_Evolve2(tree) = map_depthfirst2(f, g, 5.0, tree, Float64)


traits = OU_Evolve2(mars_tree)
plot(mars_tree, showtips = false, marker_z = traits, linewidth = 2, markercolor = :RdYlBu, size = (400, 600))

# Job from here is to be able to input the data, into the leaves and then evolve properly from those
# map_depthfirst 2 works currently, but evolves from the root to branches
# Need to modify it to evolve from the leaves
# Basically need to do independent contrasts and solve those SDEs 


# TODO Attach trait values to the leaves
# Map the tree backwards
# figure out how to get the pairs from the leaves - Pruning?


plot(mars_tree,
     size = (400, 800),
     linecolor = :orange, linewidth = 5,
     markersize = 10, markercolor = :steelblue, markerstrokecolor = :white,
     series_annotations = text.(1:nnodes(mars_tree), 5, :center, :center, :white),
     tipfont = (4,))

# Get children for each node?

#function Ancestor_Path(tree)
    #This should allow us to iterator down a tree
#    root = first(nodenamefilter(isroot, tree))
#    function local!(val, node)
#        for ch in getchildren(tree, node)
#            branch = getinbound(mars_tree, ch)
            #Need to figure out how to get branches from names, or nodes from names
#            len = getlength(tree, branch)
            #this section here allows u to get time between nodes for when u come to this later
#            tspan = (0.0, len)

            #Iterator seems to be working, now need to fix the SDE bit

            #can try to make the OU_Wrap function (Val, node, time)
#            local!(last(len), ch)
#        end
#    end
#    local!(start, root)
#end

#function Attach_branch(tree, node)
#    branch = getinbound(tree, node)
#    len = getlength(tree, branch)
#end

#get_Branches(tree) = map_depthfirst((val, node) -> Attach_branch(tree, node), 0., tree, Float64)

#get_Branches(mars_tree)

#root = first(nodenamefilter(isroot, mars_tree))

#for ch in getchildren(mars_tree, root)
    
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







    dent = Vector{Float64}()

for i in 1:length(getleafnames(mars_tree))
    idx = findall(x -> x == getleafnames(mars_tree)[i], mars_avg[:,1])
    #idx = parse.(Int, idx)
    if length(idx) >= 1
        idx = idx[1]
    end
    print(idx)
    val = mars_avg[idx,2]
    #print(val)
    push!(dent, collect(val))
end

dent

