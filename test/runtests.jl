using QuTiP, PyCall
using Base.Test
@pyimport qutip as qt

# write your own tests here
# @test 1 == 1

include("qobj.jl")
include("eseries.jl")


@test  sigmax() * basis(2,0) == qt.sigmax() * qt.basis(2,0)

@test sigmay() * basis(2,0) == qt.sigmay() * qt.basis(2,0)

@test sigmaz() * basis(2,0) == qt.sigmaz() * qt.basis(2,0)

@test qt.sigmay() == im * sigmax() * sigmaz()

q = Qobj([1,0])
@test (q[:isherm], q[:type]) == (false, "ket")

@test cnot() == qt.cnot()

@test hadamard_transform() == qt.hadamard_transform()

X = sigmax()
Y = 1 + dag(X) + dag(X) + X + 1
@test full(Y) == [2.+0im 3.; 3. 2.]

es3 = eseries([0.5*sigmaz(), 0.5*sigmaz()], [1im, -1im]) + eseries([-0.5im*sigmax(), 0.5im*sigmax()], [1im, -1im]) 
rho = fock_dm(2, 1)
es3_expect = expect(rho, es3)
@test es3_expect[:shape] == Any[1,1]

q = tensor(qidentity(2), basis(2))
s_prep = sprepost(q, dag(q))
@test full(tensor_contract(to_super(cnot()), (1, 3)) * s_prep) == Complex{Float64}[1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]


@test full(qnum(3, offset=10)) == Complex{Float64}[10 0 0; 0 11 0; 0 0 12];
