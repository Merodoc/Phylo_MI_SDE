#unsure on my pointers so right now we are being over the top

using Plots
mutable struct Particle
    rx
    ry
    vx
    vy
    radius
    mass
    count
end


function move(mol::Particle, dt)
    mol.rx += mol.vx*dt
    mol.ry += mol.vy*dt
    return mol.rx, mol.ry
end


#time to hits are specific to particles hence vy not V_y
function TTHV(mol, wall_length=1)
    #Time to hit vertical wall
    if (mol.vy > 0)
        return (wall_length-mol.ry-mol.radius)/mol.vy
    elseif (mol.vy < 0)
        return (mol.radius-mol.ry)/mol.vy
    else
        return Inf
    end
end

function TTHH(mol, wall_length =100)
    #Time to hit horizontal wall
    if (mol.vx > 0)
        return (wall_length-mol.rx-mol.radius)/mol.vx
    elseif (mol.vx < 0)
        return (mol.radius-mol.ry)/mol.vx
    else
        return Inf
    end
end

BounceV(mol) = mol.vx*-1
bounceH(mol) = mol.vy*-1

#need to determine our iterator
function TimeToHit(mol1, mol2, dt)
    #atm test with tuples put into the function
    #there is probably a matrix way to get a full matrix of dvs/drs
    # or we make a particle function that returns things about particles
    if (mol1 == mol2)
        return Inf
    else
    dx = mol1.rx-mol2.rx
    dy = mol1.ry-mol2.ry
    dvx = mol1.vx-mol2.vx
    dvy = mol1.vy - mol2.vy
    dvdr = dx*dvx + dy*dvy


    if (dvdr > 0)
        return Inf
    else
        dvdv = dvx*dvx + dy*dvy
        if (dvdv == 0)
            return Inf
        else
            drdr = dx*dx + dy*dy
            sigma = mol1.radius - mol2.radius
            d = (dvdr*dvdr) - dvdv*(drdr-sigma*sigma)
                if (d < 0)
                    return Inf
                else        
                    return -(dvdr + sqrt(d))/dvdv
                end
        end
    end
end
end

function ParticleBounce(mol1, mol2)
    #come back to this, we need to learn to use abstract types
dx = mol1.rx - mol2.rx
dy = mol1.ry - mol2.ry
dvx = mol1.vx - mol2.vx
dvy = mol1.vy - mol2.vy
dvdr = dx*dvx + dy*dvy
dist = mol1.radius + mol2.radius
J = (2*mol1.mass*mol2.mass*dvdr)/((mol1.mass + mol2.mass) * dist)
Jx = J*dx/dist
Jy = J*dy/dist
mol1.vx += Jx/mol1.mass
mol1.vy += Jy/mol1.mass
mol2.vx += Jx/mol2.mass
mol2.vy += Jy/mol2.mass
end

function ParticleBounceSurf!(dvdr = 0.5)
    mass = range(0,10,step =0.1)
    dist = range(0,10,step =0.1)
    f(x,y) = (2 .*x .*x .* dvdr/(2 .*x .*y))/y
    surface(mass, dist, f)
end
# Event Class Stuff

function dvdrSurf(wall = 100, vlim = 10, mass = 1, Rad = 1)
    dx = range(-wall/2, wall/2, length = 1000)
    dy = range(-wall/2, wall/2, length = 1000)
    dvx = range(0, vlim, length =1000)
    dvy = range(0, vlim, length =1000)
    dvdr = dx .* dvx + dy.*dvy
    J = 2 .* dvdr/(4)
    Jx = J .* dx/2
    surface(dx, dvx, Jx)
end

mutable struct Event
    time
    thisa::Particle
    thisb::Particle
end

function isValid(Event, a, b)
    if (a.count != Event.thisa.count)
        return false
    elseif (b.count != Event.thisb.count)
        return false
    else
        return true
    end
end

#Initialize random particles

function init(N, mass = 1, radius = 1)
    ret = []
    for i in 1:N
        particle = Particle(rand()*100, rand()*100, rand(1:10), rand(1:10), radius, mass, 0)
        push!(ret, particle)
    end


    #right now we are playing with constant mass/radius
    #there are options to change this but it will get complicated
    return ret
end

particles = init(100)

function particleSim(particles, time=10, dt = 0.1, vwall = 100, hwall = 100)
    p1x = []
    p1y = []
    for i in 0:dt:time
        pdΔt = deepcopy(particles)
        for j in 1:length(particles)
            move(pdΔt[j],dt)
        end
        for j in 1:length(pdΔt)
            bounces = 0
            if pdΔt[j].ry >= vwall || pdΔt[j].ry <= 0
                particles[j].vy = BounceV(particles[j])
                bounces +=1 #this sucks fix it
            end
            if pdΔt[j].rx >= hwall || pdΔt[j].rx <= 0
                particles[j].vx = bounceH(particles[j])
                bounces += 1#see above
            end
            #= Elastic Collisions are being real fucky
                for k in j:length(pdΔt)
                    if abs(pdΔt[j].ry - pdΔt[k].ry) <= (pdΔt[j].radius + pdΔt[k].radius) && abs(pdΔt[j].rx - pdΔt[k].rx) <= (pdΔt[j].radius + pdΔt[k].radius)
                        if bounces < 1
                            ParticleBounce(particles[j], particles[k])
                        end
                    end
                end
=#
            end
        end
        for j in 1:length(particles)
        move(particles[j], dt)
        end
    end


test = init(10)
pddt = deepcopy(test)
dt = 0.1

for j in 1:length(pddt)
    move(pddt[j], dt)
end

for j in 1:length(pddt)
    if pddt[j].ry >= 100 || pddt[j].ry <= 0
        test[j].vy = BounceV(test[j])
    end
    if pddt[j].rx >=100 || pddt[j].rx <= 0
        test[j].vx = bounceH(test[j])
    end
    for k in j:length(pddt)
        if abs(pddt[j].ry - pddt[k].ry) <= (pddt[j].radius + pddt[k].radius) && abs(pddt[j].rx - pddt[k].rx) <= (pddt[j].radius + pddt[k].radius)
            ParticleBounce(test[j], test[k])
        end
    end
end

for j in 1:length(test)
    move(test[j], dt)
end

test

ptest = particleSim(particles)

