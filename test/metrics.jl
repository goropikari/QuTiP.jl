x = fock_dm(5,3)
y = coherent_dm(5,1)
@test abs(fidelity(x,y) - 0.24104350624628332) < 1e-5
@test abs(tracedist(x,y) - 0.9705143161472971) < 1e-5

