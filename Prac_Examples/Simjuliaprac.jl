
using ResumableFunctions
using ConcurrentSim

@resumable function car(env::Environment)
    while true
        println("Start Parking and charging at ", now(env))
        charge_duration = 5.0
        charge_process = @process charge(sim, charge_duration)
        try
            @yield charge_process
        catch
            println("Was interrupted. Hopefully, the battery is full enough ...")
        end
        println("Start driving at ", now(env))
        trip_duration = 2.0
        @yield timeout(sim, trip_duration)
    end
end

@resumable function charge(env::Environment, duration::Number)
    @yield timeout(env, duration)
end

@resumable function driver(env::Environment, car_process::Process)
    @yield timeout(env, 3)
    @yield interrupt(car_process)
end

sim = Simulation()
car_process = @process car(sim)
@process driver(sim, car_process)

run(sim, 15)

@resumable function car(env::Environment, name::Int, bcs::Resource, driving_time::Number, charge_duration::Number)
    @yield timeout(sim, driving_time)
    println(name, " arriving at ", now(env))
    @yield request(bcs)
    println(name, " starting to charge at ", now(env))
    @yield timeout(sim, charge_duration)
    println(name, " leaving the bcs at ", now(env))
    @yield unlock(bcs)
end

sim = Simulation()
bcs = Resource(sim, 2)

for i in 1:4
    @process car(sim, i, bcs, 2i, 5)
end

run(sim)