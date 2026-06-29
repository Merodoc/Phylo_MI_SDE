using Phylo
include("BridgeFuncs.jl")
nu = Nonultrametric(5);
tree = rand(nu)

plot(tree)
getroot(tree)

root = collect(nodenamefilter(isroot, tree))
branches = getoutbounds(tree, root[1])

getlength(tree, branches[1])

setnodedata!(tree, root[1], "Trait", 0.)

Iterator = collect(nodeiter(tree))

for node in Iterator
    getnodedata(tree, node)
end



nu = Nonultrametric(5);
tree = rand(nu)

iter = traversal(tree, preorder)

for node in iter 
    if isroot(tree, node)
        setnodedata!(tree, node, "Trait", 0.)
    elseif hasinbound(tree, node)
        inbound = getinbound(tree, node)
        len = getlength(tree, inbound)
        parent = getparent(tree, node)
        start = getnodedata(tree, parent)
        startval = start["Trait"]
        print(startval)
        W = sample(0:0.01:len, Wiener(), startval)
        setnodedata!(tree, node, "Trait", last(W.yy))
    end
end




function Treesim(Model, ntips, rootval = 0.)
    nu = Ultrametric(ntips);
    tree = rand(nu)
    
    iter = traversal(tree, preorder)
    
    for node in iter 
        if isroot(tree, node)
            setnodedata!(tree, node, "Trait", rootval)
        elseif hasinbound(tree, node)
            inbound = getinbound(tree, node)
            len = getlength(tree, inbound)
            parent = getparent(tree, node)
            start = getnodedata(tree, parent)
            startval = start["Trait"]
            W = sample(0:0.01:len, Wiener(), startval)
            X = solve(EulerMaruyama(), startval, W, Model)
            setnodedata!(tree, node, "Trait", last(X.yy))
        end
    end

    return tree
end

function TreesimBM(ntips)
    nu = Ultrametric(ntips);
    tree = rand(nu)
    
    iter = traversal(tree, preorder)
    
    for node in iter 
        if isroot(tree, node)
            setnodedata!(tree, node, "Trait", 0.)
        elseif hasinbound(tree, node)
            inbound = getinbound(tree, node)
            len = getlength(tree, inbound)
            parent = getparent(tree, node)
            start = getnodedata(tree, parent)
            startval = start["Trait"]
            W = sample(0:0.01:len, Wiener(), startval)
            setnodedata!(tree, node, "Trait", last(W.yy))
        end
    end

    return tree
end

OUtree = Treesim(OUMean(1.,0.,0.1), 32)

OUDict = Dict()

for nodename in getnodes(OUtree)
    dat = getnodedata(OUtree, nodename)
    val = dat["Trait"]
    OUDict[nodename] = val 
end

OUDict

plot(OUtree, line_z = collect(values(OUDict)))

values(OUDict)

CIRtree = Treesim(CIR(1.,0.,0.1), 32)
WFtree = Treesim(WF(1., 0., 0.1), 32)

# Need to put the leaf values into a .csv for the test simulations