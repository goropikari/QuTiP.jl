using QuTiP
using Base.Test

# write your own tests here
# @test 1 == 1

x = sigmax() * basis(2,0)
@test vec(x[:full]()) == Complex{Float64}[0,1]

x = sigmay() * basis(2,0)   
@test vec(x[:full]()) == Complex{Float64}[0,im]

x = sigmaz() * basis(2,0)   
@test vec(x[:full]()) == Complex{Float64}[1,0]

@test sigmay() == im * sigmax() * sigmaz()

q = Qobj([1,0])
@test (q[:isherm], q[:type]) == (false, "ket")

@test cnot()[:full]() == Complex{Float64}[1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]

@test hadamard_transform() == Qobj(1/sqrt(2) * [1 1; 1 -1])

X = sigmax()
Y = 1 + dag(X) + dag(X) + X + 1
@test Y[:full]() == [2.+0im 3.; 3. 2.]

es3 = eseries([0.5*sigmaz(), 0.5*sigmaz()], [1im, -1im]) + eseries([-0.5im*sigmax(), 0.5im*sigmax()], [1im, -1im]) 
rho = fock_dm(2, 1)
es3_expect = expect(rho, es3)
@test es3_expect[:shape] == Any[1,1]

q = tensor(qidentity(2), basis(2))
s_prep = sprepost(q, dag(q))
@test (tensor_contract(to_super(cnot()), (1, 3)) * s_prep)[:full]() == Complex{Float64}[1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]


@test qnum(3, offset=10)[:full]() == Complex{Float64}[10 0 0; 0 11 0; 0 0 12];