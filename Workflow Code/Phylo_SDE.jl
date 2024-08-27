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
            branch = getinbound(mars_tree, ch)
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
mars_tree = open(parse(RootedTree), Phylo.path("H:/Code/Mars_TimeTree.nwk"))

plot(mars_tree)

using CSV
using DataFrames

# Read in Data
# Will need to do similar data matching to that in R
mars = CSV.read("H:/Code/mars.csv", DataFrame)
mars_avg = CSV.read("H:/Code/mars_avg.csv", DataFrame)
using Statistics


#Convert Dataframe rows to numerics and remove missing data (should not be a problem later obviously)
dentary = parse.(Float64,mars_avg.dentary)
filter!(row -> row.dentary != "NA", mars_avg)
dentary = parse.(Float64,mars_avg[!,2])

#Define the SDE for map_depthfirst2
W = WienerProcess(0.0,0.0,0.0)
α = 1.0
σ = 0.1
μ = mean(dentary)
μ0 = μ
f(u, p, t) =  α*(μ-u)
g(u,p, t) = σ*u

#Model Traits using the OU Model (This is purely exploratory, no data involved)
OU_Evolve2(tree) = map_depthfirst2(f, g, 5.0, tree, Float64)


traits = OU_Evolve2(mars_tree)
plot(mars_tree, marker_z = traits, linecolor = :RdYlBu, linewidth = 5, showtips = false)

# Job from here is to be able to input the data, into the leaves and then evolve properly from those
# map_depthfirst 2 works currently, but evolves from the root to branches
# Need to modify it to evolve from the leaves
# Basically need to do independent contrasts and solve those SDEs 
