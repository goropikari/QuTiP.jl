@test basis(2,0) == qt.basis(2,0)
@test 1 + sigmax() + 1 == PyObject(1) + qt.sigmax() + PyObject(1)
@test dag(basis(2,0)) == qt.basis(2,0)[:dag]()

A = sigmax() + sigmay()
B = qt.sigmax() + qt.sigmay()
@test A == B
@test dims(A) == B[:dims]
@test shape(A) == B[:shape]
@test superrep(A) == B[:superrep]
@test isherm(A) == B[:isherm]
@test iscp(A) == B[:iscp]
@test ishp(A) == B[:ishp]
@test istp(A) == B[:istp]
@test iscptp(A) == B[:iscptp]
@test isket(A) == B[:isket]
@test isbra(A) == B[:isbra]
@test isoper(A) == B[:isoper]
@test issuper(A) == B[:issuper]
@test isoperket(A) == B[:isoperket]
@test isoperbra(A) == B[:isoperbra]

@test conj(A) == B[:conj]()

@test haskey(A, "isherm")
@test A - 1 == B - 1
@test 1 - A == B * (-1) + 1
@test -1 + A == B + (-1)
@test - A + 1 == B * (-1) + 1
@test A * 9 == B * 9
@test A / 9 == B / 9
@test A.o + A == B * 2
@test A.o - A == B - B
@test A.o * A == B * B
@test A + A.o == B * 2
@test A - A.o == B - B
@test A * A.o == B * B
@test A + A == B + B
@test A - A == B - B
@test A * A == B * B
@test +A == A
@test -A == A * (-1)
