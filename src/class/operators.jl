export jmat, spin_Jx, spin_Jy, spin_Jz, spin_Jm, spin_Jp, spin_J_set, sigmap, sigmam, sigmax, sigmay, sigmaz, destroy, create, qeye, position, momentum, num, squeeze, squeezing, displace, commutator, qutrit_ops, qdiags, phase, zero_oper, enr_destroy, enr_idenitty
const operators_class = (:jmat, :spin_Jx, :spin_Jy, :spin_Jz, :spin_Jm, :spin_Jp, :spin_J_set, :sigmap, :sigmam, :sigmax, :sigmay, :sigmaz, :destroy, :create, :qeye, :position, :momentum, :num, :squeeze, :squeezing, :displace, :commutator, :qutrit_ops, :qdiags, :phase, :zero_oper, :enr_destroy, :enr_idenitty)


export qidentity
sf = string(:identity)
# @eval @doc LazyHelp(qutip,$sf) function qidentity(args...; kws...)
# 	if !haskey(qutip, $sf)
# 		error("qutip ", version, " does not have qutip", $sf)
# 	end
# 	return pycall(qutip[$sf], PyAny, args...; kws...)
# end

@eval @doc LazyHelp(qutip,$sf) function qidentity(args)
	if !haskey(qutip, $sf)
		error("qutip ", version, " does not have qutip", $sf)
	end
    if  typeof(args) <: Integer
        return pycall(qutip[$sf], PyAny, args)
    elseif typeof(args) <: Array
        # return pycall(qutip[$sf], PyAny, [[i] for i in args])
        return pycall(qutip[$sf], PyAny, map(x -> [x], args))
    end
end
