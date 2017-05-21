__precompile__()

module QuTiP 
using PyCall
import PyCall: PyNULL, pyimport_conda, pycall
export qutip

using Compat
@compat import Base.show

###########################################################################
# quoted from PyPlot.jl
# Julia 0.4 help system: define a documentation object
# that lazily looks up help from a PyObject via zero or more keys.
# This saves us time when loading PyPlot, since we don't have
# to load up all of the documentation strings right away.
immutable LazyHelp
    o::PyObject
    keys::Tuple{Vararg{Compat.String}}
    LazyHelp(o::PyObject) = new(o, ())
    LazyHelp(o::PyObject, k::AbstractString) = new(o, (k,))
    LazyHelp(o::PyObject, k1::AbstractString, k2::AbstractString) = new(o, (k1,k2))
    LazyHelp(o::PyObject, k::Tuple{Vararg{AbstractString}}) = new(o, k)
end
@compat function show(io::IO, ::MIME"text/plain", h::LazyHelp)
    o = h.o
    for k in h.keys
        o = o[k]
    end
    if haskey(o, "__doc__")
        print(io, convert(AbstractString, o["__doc__"]))
    else
        print(io, "no Python docstring found for ", h.k)
    end
end
Base.show(io::IO, h::LazyHelp) = @compat show(io, "text/plain", h)
function Base.Docs.catdoc(hs::LazyHelp...)
    Base.Docs.Text() do io
        for h in hs
            @compat show(io, MIME"text/plain"(), h)
        end
    end
end

###########################################################################


const qutip = PyNULL()
# copy!(qutip, pyimport_conda("qutip", "qutip"))

# ref
# https://github.com/JuliaPy/PyPlot.jl/blob/master/src/PyPlot.jl#L166
# for f in plt_funcs
#     sf = string(f)
#     @eval @doc LazyHelp(plt,$sf) function $f(args...; kws...)
#         if !haskey(plt, $sf)
#             error("matplotlib ", version, " does not have pyplot.", $sf)
#         end
#         return pycall(plt[$sf], PyAny, args...; kws...)
#     end
# end


# export ducumented qutip API
# states
export ghz_state, w_state, triplet_states, singlet_state, bell_state, spin_coherent, spin_state, zero_ket, phase_basis, enr_thermal_dm, enr_fock, enr_state_dictionaries, state_number_qobj, state_index_number, state_number_index, state_number_enumerate, bra, ket, qstate, projection, ket2dm, maximally_mixed_dm, thermal_dm, fock, fock_dm, coherent_dm, coherent, qutrit_basis, basis 

# random_objects
export rand_stochastic, rand_super_bcsz, rand_super, rand_kraus_map, rand_dm_hs, rand_dm_ginibre, rand_dm, ket_haar, rand_ket, rand_unitary_haar, rand_unitary, rand_herm, randnz, rand_jacobi_rotation

# continuous_variables
export logarithmic_negativity, wigner_covariance_matrix, correlation_matrix_quadrature, correlation_matrix_field, covariance_matrix, correlation_matrix

# superoperator
export sprepost, spre, spost, mat2vec_index, vec2mat_index, vec2mat, mat2vec, vector_to_operator, operator_to_vector, lindblad_dissipator, liouvillian_ref, liouvillian

# superop_reps
export to_stinespring, to_kraus, to_super, to_chi, to_choi, choi_to_stinespring, chi_to_choi, choi_to_chi, kraus_to_super, kraus_to_choi, choi_to_kraus, choi_to_super, super_to_choi

# subsystem_apply
export subsystem_apply

# operatos
export jmat, spin_Jx, spin_Jy, spin_Jz, spin_Jm, spin_Jp, spin_J_set, sigmap, sigmam, sigmax, sigmay, sigmaz, destroy, create, qeye, identity, position, momentum, num, squeeze, squeezing, displace, commutator, qutrit_ops, qdiags, phase, zero_oper, enr_destroy, enr_idenitty


const states = (:ghz_state, :w_state, :triplet_states, :singlet_state, :bell_state, :spin_coherent, :spin_state, :zero_ket, :phase_basis, :enr_thermal_dm, :enr_fock, :enr_state_dictionaries, :state_number_qobj, :state_index_number, :state_number_index, :state_number_enumerate, :bra, :ket, :qstate, :projection, :ket2dm, :maximally_mixed_dm, :thermal_dm, :fock, :fock_dm, :coherent_dm, :coherent, :qutrit_basis, :basis)

const random_objects = (:rand_stochastic, :rand_super_bcsz, :rand_super, :rand_kraus_map, :rand_dm_hs, :rand_dm_ginibre, :rand_dm, :ket_haar, :rand_ket, :rand_unitary_haar, :rand_unitary, :rand_herm, :randnz, :rand_jacobi_rotation)

const continuous_variables = (:logarithmic_negativity, :wigner_covariance_matrix, :correlation_matrix_quadrature, :correlation_matrix_field, :covariance_matrix, :correlation_matrix)

const superoperator = (:sprepost, :spre, :spost, :mat2vec_index, :vec2mat_index, :vec2mat, :mat2vec, :vector_to_operator, :operator_to_vector, :lindblad_dissipator, :liouvillian_ref, :liouvillian
)

const superop_reps = (:to_stinespring, :to_kraus, :to_super, :to_chi, :to_choi, :choi_to_stinespring, :chi_to_choi, :choi_to_chi, :kraus_to_super, :kraus_to_choi, :choi_to_kraus, :choi_to_super, :super_to_choi)

const subsystem_apply = (subsystem_apply, )

const operators = (:jmat, :spin_Jx, :spin_Jy, :spin_Jz, :spin_Jm, :spin_Jp, :spin_J_set, :sigmap, :sigmam, :sigmax, :sigmay, :sigmaz, :destroy, :create, :qeye, :identity, :position, :momentum, :num, :squeeze, :squeezing, :displace, :commutator, :qutrit_ops, :qdiags, :phase, :zero_oper, :enr_destroy, :enr_idenitty)

const qutipfn = (states..., 
				random_objects...,
                continuous_variables...,
                superoperator..., 
                superop_reps..., 
                subsystem_apply..., 
                operators...)
for f in qutipfn
    sf = string(f)
    @eval @doc LazyHelp(qutip,$sf) function $f(args...; kws...)
        if !haskey(qutip, $sf)
            error("qutip ", version, " does not have qutip.", $sf)
        end
        return pycall(qutip[$sf], PyAny, args...; kws...)
    end
end


function __init__()
    copy!(qutip, pyimport_conda("qutip", "qutip"))
    global const version = try
        convert(VersionNumber, qutip[:__version__])
    catch
        v"0.0" # fallback
    end
end



end # module

