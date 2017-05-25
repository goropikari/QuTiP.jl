export jmat, spin_Jx, spin_Jy, spin_Jz, spin_Jm, spin_Jp, spin_J_set, sigmap, sigmam, sigmax, sigmay, sigmaz, destroy, create, qeye, momentum, squeezing, displace, commutator, qutrit_ops, qdiags, phase, zero_oper, enr_destroy, enr_idenitty
const operators_class = (:jmat, :spin_Jx, :spin_Jy, :spin_Jz, :spin_Jm, :spin_Jp, :spin_J_set, :sigmap, :sigmam, :sigmax, :sigmay, :sigmaz, :destroy, :create, :qeye, :momentum, :squeezing, :displace, :commutator, :qutrit_ops, :qdiags, :phase, :zero_oper, :enr_destroy, :enr_idenitty)

# To avoid conflict with Base module
export qidentity, qnum, qposition, qsqueeze
const renamedfn = (:identity, :num, :position, :squeeze)
for f in renamedfn
    sf = string(f)
    nf = Symbol("q", f)
    @eval @doc LazyHelp(qutip,$sf) function $nf(args...; kws...)
        if !haskey(qutip, $sf)
            error("qutip ", version, " does not have qutip.", $sf)
        end
        return pycall(qutip[$sf], Quantum, args...; kws...)
    end
end
