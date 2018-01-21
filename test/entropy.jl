psi = ghz_state(2)
rho = ptrace(psi * psi', 0)
@test abs(1.0 - concurrence(psi)) < 1.0e-10
@test entropy_vn(rho, 2) == 1.0

rho=0.5*fock_dm(2,0)+0.5*fock_dm(2,1)
@test entropy_linear(rho) == 0.5
