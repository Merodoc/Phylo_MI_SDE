using Bridge
using Plots

# Define a diffusion process 

struct OrnsteinUhlenbeck <: ContinuousTimeProcess{Float64}
    β::Float64 #drift parameter 
    σ::Float64 #diffusion parameter
end

#define drift and diffusion coefficient
Bridge.b(t,x,P::OrnsteinUhlenbeck) = -P.β*x
Bridge.σ(t,x,P::OrnsteinUhlenbeck) = P.σ

#simulate ornstein uhlenbeck using Euler scheme

W = sample(0:0.01:10, Wiener())
X = solve(EulerMaruyama(), 0.1, W, OrnsteinUhlenbeck(2.0,1.0))
plot(X, level = "X")
X

p = plot(X, level = "X")

p = plot(1)
X2t = reverse(X.tt)
X2y = reverse(X.yy)
anim = @animate for x = 1:length(X.tt)
    push!(p, 1, W.tt[x], W.yy[x])
    #push!(p, 2, X2t[x], X2y[x])
end

gif(anim)

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
    
test = WienerAnim(10, 1, 10)

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

function WienerAnim3(t, dt, n)
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
        Wyforward = W.yy
        Wybackward = reverse(W.tt)
        #Have to do a lot to fix this, as the y values are going same direction
        #need to account for the reverse of the time, not reverse of the value
        #Otherwise should be manageable 
        for i in 1:length(W.yy)
            if W.yy[i] > reverse(W2.yy)[i]-0.5 && W.yy[i] < reverse(W2.yy)[i] + 0.5
                println("Bridge Match")
                newW2 = reverse(W2.yy)[1:i]
                for j in i:length(W.yy)
                    push!(newW2, W.yy[j])
                end
                W2.yy = newW2
                break
            end
        end

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
    
test = WienerAnim3(10, 0.1, 50)