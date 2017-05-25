using QuTiP
using Base.Test

# write your own tests here
# @test 1 == 1

x = sigmax() * basis(2,0)
@test vec(x[:full]()) == Complex{Float64}[0,1]

q = Qobj([1,0])
@test (q[:isherm], q[:type]) == (false, "ket")

X = sigmax()
Y = 1 + dag(X) + dag(X) + X + 1
@test Y[:full]() == [2.+0im 3.; 3. 2.]
