using QuTiP.Nonmarkov
kappa = 1.0 # cavity decay rate
wc = 0.0*kappa # cavity frequency
wa = 0.0*kappa # qubit frequency
g = 10.0*kappa # coupling strength
N = 3 # size of cavity basis

# intial state
psi0c = basis(N,0)
rho0c = ket2dm(psi0c)
rho0a = ket2dm(basis(2,0))
rho0 = tensor(rho0a,rho0c)
rho0avec = operator_to_vector(rho0a)

# identity superoperator
Id = tensor(qeye(2), qeye(N))
E0 = sprepost(Id,Id)

# partial trace over the cavity, reprsented as a superoperator
ptracesuper = tensor_contract(E0,(1,3))

# intial state of the cavity, represented as a superoperator
superrho0cav = sprepost(tensor(qeye(2),psi0c),
                           tensor(qeye(2),dag(psi0c)))

# operators
a  = tensor(qeye(2), destroy(N))
sm = tensor(sigmam(), qeye(N))
sz = tensor(sigmaz(), qeye(N))

# Hamiltonian
H = wc * dag(a) * a + wa * dag(sm) * sm + g * (dag(a) * sm + a * dag(sm))
c_ops = [sqrt(kappa)*a]

function dynmap(t)
    # reduced dynamical map for the qubit at time t
    Et = mesolve(H, E0, [0.,t], c_ops, [])[:states][end]
    return ptracesuper*(Et*superrho0cav)
end

exacttimes = collect(0:0.01:4.99)
exactsol = mesolve(H, rho0, exacttimes, c_ops, [])

times = 0:0.1:4.9 # total extrapolation time
ttmsols = []
maxlearningtimes = [0.5, 2.0] # maximal learning times
for T in maxlearningtimes
    learningtimes = 0:0.1:T-0.1
    learningmaps = [dynmap(t) for t in learningtimes] # generate exact dynamical maps to learn from
    push!(ttmsols, ttmsolve(learningmaps, rho0a, times)) # extrapolate using TTM
end

@test !(isempty(ttmsols))
