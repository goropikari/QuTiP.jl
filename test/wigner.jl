w0 = 1.0  * 2 * pi  # cavity frequency
wa = 1.0  * 2 * pi  # atom frequency
g  = 0.05 * 2 * pi  # coupling strength

kappa = 0.04        # cavity dissipation rate
gamma = 0.00        # atom dissipation rate
Gamma = 0.35        # atom pump rate

N = 50              # number of cavity fock states
n_th_a = 0.0        # avg number of thermal bath excitation

a  = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))
sx = tensor(qeye(N), sigmax())
H = w0 * dag(a) * a + wa * dag(sm) * sm + g * (dag(a) + a) * sx
c_ops = []
rate = kappa * (1 + n_th_a)
if rate > 0.0
    push!(c_ops, sqrt(rate) * a)
end
rate = kappa * n_th_a
if rate > 0.0
    push!(c_ops, sqrt(rate) * dag(a))
end
rate = gamma
if rate > 0.0
    push!(c_ops, sqrt(rate) * sm)
end
rate = Gamma
if rate > 0.0
    push!(c_ops, sqrt(rate) * dag(sm))
end
rho_ss = steadystate(H, c_ops)
xvec = linspace(-5,5,200)
rho_cavity = ptrace(rho_ss, 0)
W = wigner(rho_cavity, xvec, xvec)
@test W == qt.wigner(rho_cavity, xvec, xvec)
