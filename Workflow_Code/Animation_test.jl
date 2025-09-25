using Bridge
using Plots
using DataFrames
# Define a diffusion process 

struct OrnsteinUhlenbeck <: ContinuousTimeProcess{Float64}
    β::Float64 #drift parameter 
    σ::Float64 #diffusion parameter
end

#define drift and diffusion coefficient
Bridge.b(t,x,P::OrnsteinUhlenbeck) = -P.β*x
Bridge.σ(t,x,P::OrnsteinUhlenbeck) = P.σ

#simulate ornstein uhlenbeck using Euler scheme

W = sample(0:0.01:3, Wiener())
X = solve(EulerMaruyama(), 0.1, W, OrnsteinUhlenbeck(2.0,1.0))
plot(X, level = "X")
X

p = plot(X, level = "X")

p = plot(1, legend = false)
title!("Standard Brownian Motion")

xlims!(0, 2)
ylims!(-2, 2)
X2t = reverse(X.tt)
X2y = reverse(X.yy)
anim = @animate for x = 1:length(X.tt)
    push!(p, 1, W.tt[x], W.yy[x])
    #push!(p, 2, X2t[x], X2y[x])
end

gif(anim, "BMAnim.gif", fps = 30)

function WienerAnim(t, dt, n)
    Y = Vector{Vector{Float64}}()
    T = Vector{Vector{Float64}}()
    p = plot(n)
    xlims!(0,t)
    ylims!(-5, 5)
    for i in 1:n
        W = sample(0:dt:t, Wiener())
        push!(Y, W.yy)
        push!(T, W.tt)
    end
    anim = @animate for x = 1:t*dt
        for i in 1:n
            push!(p, i, T[i][x], Y[i][x])
        end
    end
    println(length(Y[1]))
    println(length(T[1]))
return anim
end
    
test = WienerAnim(2, 0.1, 50)

gif(test)


function WienerAnim2(t, dt, n)
    Y = Vector{Vector{Float64}}()
    Y2 = Vector{Vector{Float64}}()
    T2 = Vector{Vector{Float64}}()
    T = Vector{Vector{Float64}}()
    p = plot(2*n, legend=false)
    xlims!(0,t)
    ylims!(-10,10)
    for i in 1:n
        W = sample(0:dt:t, Wiener())
        W2 = sample(0:dt:t, Wiener())
        push!(Y, W.yy)
        push!(T, W.tt)
        push!(Y2, W2.yy)
        push!(T2, reverse(W2.tt))
    end
    anim = @animate for x = 1:length(T[1])
        for i in 1:n
            push!(p, i, T[i][x], Y[i][x])
            push!(p, (2*n+1)-i, T2[i][x], Y2[i][x])
        end
    end
return anim
end
    
test = WienerAnim2(10, 0.1, 50)

gif(test)

function Phylo_Bridge_anim(start, fin, fin_time, anc, dt, samples)
    #start = start value
    #fin = final value
    #fin_time = total time
    #anc = time point of the ancestor node
    #samples = number of repeats
    N = 1:samples
    Xhat = Vector{Float64}()
    paths = DataFrame()
    for n in N
        B = sample(0:dt:fin_time, WienerBridge(fin_time,fin), start)
        idx = findall(x -> x == anc, B.tt)
        val = B.yy[idx][1]
        paths[!, string(n)] = Vector(B.yy)
        end

return paths
end

bridgeanim = Phylo_Bridge_anim(0., 0., 5., 2.5, 0.01, 10)

t = 0:0.01:5
p = plot(200, legend = false)
xlims!(0, 5)
ylims!(-5, 5)
anim = @animate for x = 1:250
    for i in 1:10
        y = collect(eachcol(bridgeanim)[i])
        if x == 250
            push!(p, i, t[x], y[x])
        end
        push!(p, i, t[x], y[x])
        push!(p, 201-i, reverse(t)[x], reverse(y)[x])
    end
end

gif(anim)


gif(anim, "BMAnimdasdas10.gif", fps = 60)


B = sample(0:0.1:2, WienerBridge(2,0), start)

bridgeanim

t = 0:0.01:5
p = plot(200, legend = false)
xlims!(0, 5)
ylims!(-5, 5)
    for i in 1:10
        y = collect(eachcol(bridgeanim)[i])
        if x == 250
            push!(p, i, t[x], y[x])
        end
        push!(p, i, t[x], y[x])
        push!(p, 201-i, reverse(t)[x], reverse(y)[x])
    end
end
