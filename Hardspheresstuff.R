#Hard-Spheres

#Set up hard sphers/disks on a lattice

#assign random initial velocities to N particles
# E_k = 3NkT/2 with some desired T
# Total Momentum = 0

#For each pair of particles (i,j) in the system get t_ij
#Time for pair to meet

t_ij <- (-b - sqrt(b^2-v^2*(r^2-d^2))/v^2)

# d is the sphere diameter,
# r is the distance between centers of i and j

b = (r_j -r_i).(v_j-v_i)
v = abs(v_j-v_i)

# solutions to t_ij

q = -b[b +sign(b)*sqrt(b^2-v^2*(r^2-d^2))]

t_ij = min(q/v^2, (r^2-d^2)/q)

#set up 2 arrays that contain
#for each particle, the smallest positie collision time
# t(i) = min(t_ij)
# and the next collision partner j(i)

#if no collision partner then j(i) = 0 and t(i) = inf

#this loop over all N(N-1)/2 pairs need only be performed once

#find the smallest element in table of free flight times
#t(i_0)
#this is the time until the next encounter between two particles
# indices are i_0 and j_0

#during the time t(i_0) all particles perform a free flight

r_i -> r_i + v_i.t(i_0)
t_i -> t(i) - t(i_0)

#now an elastic collision between i = i_0 and j = j_0 occurs
dela_v = b*r_ij/d^2

vprime_i = v_i + deltav
vprime_j = v_j - deltav

# now repeat
