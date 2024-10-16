include("Leaf_Prune.jl")


mars_tree = open(parse(RootedTree), Phylo.path("C:/PhD/Phylo_MI_SDE/Workflow_Code/newtree.nwk"))
mars_avg = CSV.read("C:/PhD/Phylo_MI_SDE/Workflow_Code/Imp_Mars_PVR25new5.csv", DataFrame, types = [String, Float64])

start_vals = mars_avg.mean
pruned = Leaf_Prune(mars_tree, start_vals)

