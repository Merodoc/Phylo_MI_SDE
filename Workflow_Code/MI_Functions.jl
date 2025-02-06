#MI_Dir = "C:/Users/Rowan/OneDrive/Documents/GitHub/REG_PhD/Workflow_Code/MI_Data2lonly/"
3MIdict = Phybridge_Dict(MI_Dir)

#MIdata = MIdict["MarsMI"]

#MIdata.femur1


function Phybridge_MIDict(dir)
    try 
        readdir(dir)
    catch
        println("Invalid directory")
    end
    Files = readdir(dir)

    traitsource = CSV.read(string(dir, Files[1]), DataFrame)
    traits = names(select(traitsource, Not([:Species, :Individual])))

    df_dict = Dict{String, DataFrame}()
# Collates all the simulated data for each trait, can do stuff to the data frames in the dictionary
for i in traits 
    #file_idx = findall(x -> i == split(x, "_")[1], Files)
    trait_data = DataFrame()
    iter = 1
    for file in Files
        data = CSV.read(string(dir, file), DataFrame)
        if iter == 1
            trait_data[!, :Species] = data[!, 1]
        end
        iter = iter + 1
        trait_data[!, file] = data[!, i]
    end
    df_dict[i] = trait_data
end
return df_dict
end 

MI_dict = Phybridge_MIDict(MI_Dir)

MI_femur = MI_dict["femur"]

length(names(MI_femur))
p = plot(legend = false)
for i in 2:length(names(MI_femur))
    vals = collect(eachcol(MI_femur)[i])
    U = kde(vals)
    plot!(U.x, U.density)
end

display(p)

osp = filter(:Species => ==("Osphranter_robustus"), femur)
osp = select(osp, Not([:Species]))
osp = collect(eachrow(osp)[1])

U = kde(osp)

plot!(U.x, U.density, linewidth = 5, linecolor = :red)

osp = filter(:Species => ==("Onychogalea_unguifera"), femur)
osp = select(osp, Not([:Species]))
osp = collect(eachrow(osp)[1])

U = kde(osp)

plot(U.x, U.density, linewidth = 5, linecolor = :blue)
savefig("Onychogalea_femur")

std(osp)

