delta = 0.0 * 2 * pi
epsilon = 0.5 * 2 * pi
gamma = 0.25
H = delta/2 * sigmay() + epsilon/2 * sigmaz()
a_ops = [sigmax()]
S_w = [w -> gamma * (w >= 0)]
Hqt = qt.sigmay() * delta/2 + qt.sigmaz() * epsilon/2
a_opsqt = [qt.sigmax()]
@test bloch_redfield_tensor(H, a_ops, S_w) == qt.bloch_redfield_tensor(Hqt, a_opsqt, S_w)
