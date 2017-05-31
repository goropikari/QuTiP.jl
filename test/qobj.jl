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

