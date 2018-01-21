# From qutip 4.2.0, https://github.com/qutip/qutip/tree/v4.2.0

# qip/gate.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/qip/gates.py
export rx, ry, rz, sqrtnot, snot, phasegate, cphase, cnot,
csign, berkeley, swapalpha, swap, iswap, sqrtswap,
sqrtiswap, fredkin, toffoli, rotation, controlled_gate,
globalphase, hadamard_transform, gate_sequence_product,
gate_expand_1toN, gate_expand_2toN, gate_expand_3toN,
qubit_clifford_group
const gate_module = (:rx, :ry, :rz, :sqrtnot, :snot, :phasegate, :cphase, :cnot,
                     :csign, :berkeley, :swapalpha, :swap, :iswap, :sqrtswap,
                     :sqrtiswap, :fredkin, :toffoli, :rotation, :controlled_gate,
                     :globalphase, :hadamard_transform, :gate_sequence_product,
                     :gate_expand_1toN, :gate_expand_2toN, :gate_expand_3toN,
                     :qubit_clifford_group)


# about.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/about.py
## export about
const about_module = (:about, )

# bloch.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/about.py
export Bloch
const bloch_module = (:Bloch, )

# bloch3d.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/bloch3d.py
export Bloch3d
const bloch3d_module = (:Bloch3d, )

# bloch_redfield.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/bloch_redfield.py
# bloch_redfield_tensor is defined individually.
export brmesolve, bloch_redfield_solve
const bloch_redfield_module = (:brmesolve, :bloch_redfield_solve)

# configrc.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/configrc.py
# There is no __all__.

# continuous_variables.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/configrc.py
export correlation_matrix, covariance_matrix,  correlation_matrix_field, correlation_matrix_quadrature, wigner_covariance_matrix, logarithmic_negativity
const continuous_variables_module = (:correlation_matrix, :covariance_matrix, :correlation_matrix_field, :correlation_matrix_quadrature, :wigner_covariance_matrix, :logarithmic_negativity)

# correlation.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/correlation.py
export correlation_2op_1t, correlation_2op_2t, correlation_3op_1t, correlation_3op_2t, coherence_function_g1, coherence_function_g2, spectrum, spectrum_correlation_fft, correlation_ss, correlation, correlation_4op_1t, correlation_4op_2t, spectrum_ss, spectrum_pi
const correlation_module = (:correlation_2op_1t, :correlation_2op_2t, :correlation_3op_1t, :correlation_3op_2t, :coherence_function_g1, :coherence_function_g2, :spectrum, :spectrum_correlation_fft, :correlation_ss, :correlation, :correlation_4op_1t, :correlation_4op_2t, :spectrum_ss, :spectrum_pi)

# countstat.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/countstat.py
export countstat_current, countstat_current_noise
const countstat_module = (:countstat_current, :countstat_current_noise)

# dimensions.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/dimensions.py
# Everything are explicitly imported.


# distributions.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/distributions.py
export Distribution, WignerDistribution, QDistribution, TwoModeQuadratureCorrelation, HarmonicOscillatorWaveFunction, HarmonicOscillatorProbabilityFunction
const distributions_module = (:Distribution, :WignerDistribution, :QDistribution, :TwoModeQuadratureCorrelation, :HarmonicOscillatorWaveFunction, :HarmonicOscillatorProbabilityFunction)

# entropy.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/entropy.py
export entropy_vn, entropy_linear, entropy_mutual, negativity, concurrence, entropy_conditional, entangling_power
const entropy_module = (:entropy_vn, :entropy_linear, :entropy_mutual, :negativity, :concurrence, :entropy_conditional, :entangling_power)

# eseries.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/eseries.py
# esspec and esval are defined individually.
export eseries, estidy
const eseries_module = (:eseries, :estidy)

# essolve.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/essolve.py
# essolve is defined individually.
export ode2es
const essolve_module = (:ode2es, )

# expect.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/expect.py
# expect is defined individually.
export variance
const expect_module = (:variance, )

# fastsparse.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/fastsparse.py
# There is no __all__.

# fileio.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/fastsparse.py
export file_data_store, file_data_read, qsave, qload
const fileio_module = (:file_data_store, :file_data_read, :qsave, :qload)

# floquet.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/floquet.py
export floquet_modes, floquet_modes_t, floquet_modes_table, floquet_modes_t_lookup, floquet_states, floquet_states_t, floquet_wavefunction, floquet_wavefunction_t, floquet_state_decomposition, fsesolve, floquet_master_equation_rates, floquet_collapse_operators, floquet_master_equation_tensor, floquet_master_equation_steadystate, floquet_basis_transform, floquet_markov_mesolve, fmmesolve
const floquet_module = (:floquet_modes, :floquet_modes_t, :floquet_modes_table, :floquet_modes_t_lookup, :floquet_states, :floquet_states_t, :floquet_wavefunction, :floquet_wavefunction_t, :floquet_state_decomposition, :fsesolve, :floquet_master_equation_rates, :floquet_collapse_operators, :floquet_master_equation_tensor, :floquet_master_equation_steadystate, :floquet_basis_transform, :floquet_markov_mesolve, :fmmesolve)

# graph.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/graph.py
export graph_degree, column_permutation, breadth_first_search,reverse_cuthill_mckee, maximum_bipartite_matching, weighted_bipartite_matching
const graph_module = (:graph_degree, :column_permutation, :breadth_first_search, :reverse_cuthill_mckee, :maximum_bipartite_matching, :weighted_bipartite_matching)

# hardware_info.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/hardware_info.py
const hardware_info_module = (:hardware_info, )

# interpolate.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/interpolate.py
export Cubic_Spline
const interpolate_module = (:Cubic_Spline, )

# ipynbtools.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/ipynbtools.py
# export version_table, plot_animation, HTMLProgressBar
# ipynbtools_module = (:version_table, :plot_animation, :HTMLProgressBar)

# logging_utils.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/logging_utils.py
export get_logger
const logging_utils_module = (:get_logger, )

# matplotlib_utilities.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/matplotlib_utilities.py
export wigner_cmap, MidpointNorm, complex_phase_cmap
const matplotlib_utilities_module = (:wigner_cmap, :MidpointNorm, :complex_phase_cmap)

# mcsolve.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/mcsolve.py
export mcsolve
const mcsolve_module = (:mcsolve, )

# mesolve.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/mesolve.py
export mesolve, odesolve
const mesolve_module = (:mesolve, :odesolve)

# metrics.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/metrics.py
export fidelity, tracedist, bures_dist, bures_angle, hilbert_dist, average_gate_fidelity, process_fidelity, unitarity, dnorm
const metrics_module = (:fidelity, :tracedist, :bures_dist, :bures_angle, :hilbert_dist, :average_gate_fidelity, :process_fidelity, :unitarity, :dnorm)

# operators.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/operators.py
# identity, position, num and squeeze are renamed in order to avoid name conflict with Base module functions.
export jmat, spin_Jx, spin_Jy, spin_Jz, spin_Jm, spin_Jp, spin_J_set, sigmap, sigmam, sigmax, sigmay, sigmaz, destroy, create, qeye, momentum, squeezing, displace, commutator, qutrit_ops, qdiags, phase, qzero, enr_destroy, enr_identity, charge, tunneling
const operators_module = (:jmat, :spin_Jx, :spin_Jy, :spin_Jz, :spin_Jm, :spin_Jp, :spin_J_set, :sigmap, :sigmam, :sigmax, :sigmay, :sigmaz, :destroy, :create, :qeye, :identity, :position, :momentum, :num, :squeeze, :squeezing, :displace, :commutator, :qutrit_ops, :qdiags, :phase, :qzero, :enr_destroy, :enr_identity, :charge, :tunneling)

# orbital.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/orbital.py
export orbital
const orbital_module = (:orbital, )

# parallel.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/parallel.py
export parfor, parallel_map, serial_map
const parallel_module = (:parfor, :parallel_map, :serial_map)

# partial_transpose.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/partial_transpose.py
export partial_transpose
const partial_transpose_module = (:partial_transpose, )

# permute.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/permute.py
export reshuffle
const permute_module = (:reshuffle, )

# propagator.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/propagator.py
export propagator, propagator_steadystate
const propagator_module = (:propagator, :propagator_steadystate)

# qobj.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/qobj.py
export Qobj, qobj_list_evaluate, ptrace, dag, isequal, issuper, isoper, isoperket, isoperbra, isket, isbra, isherm, shape, dims
const qobj_module = (:Qobj, :qobj_list_evaluate, :ptrace, :dag, :isequal, :issuper, :isoper, :isoperket, :isoperbra, :isket, :isbra, :isherm, :shape, :dims)

# random_objects.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/random_objects.py
export rand_herm, rand_unitary, rand_ket, rand_dm, rand_unitary_haar, rand_ket_haar, rand_dm_ginibre, rand_dm_hs, rand_super_bcsz, rand_stochastic, rand_super
const random_objects_module = (:rand_herm, :rand_unitary, :rand_ket, :rand_dm, :rand_unitary_haar, :rand_ket_haar, :rand_dm_ginibre, :rand_dm_hs, :rand_super_bcsz, :rand_stochastic, :rand_super)

# rcsolve.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/rcsolve.py
export rcsolve
const rcsolve_module = (:rcsolve, )

# rhs_generate.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/rhs_generate.py
export rhs_generate, rhs_clear
const rhs_generate_module = (:rhs_generate, :rhs_clear)

# semidefinite.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/semidefinite.py
# There is no __all__.

# sesolve.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/sesolve.py
export sesolve
const sesolve_module = (:sesolve, )

# settings.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/settings.py
# There is no __all__.

# simdiag.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/simdiag.py
export simdiag
const simdiag_module = (:simdiag, )

# solver.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/solver.py
export Options, Odeoptions, Odedata
const solver_module = (:Options, :Odeoptions, :Odedata)

# sparse.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/sparse.py
export sp_fro_norm, sp_inf_norm, sp_L2_norm, sp_max_norm, sp_one_norm, sp_reshape, sp_eigs, sp_expm, sp_permute, sp_reverse_permute, sp_bandwidth, sp_profile
const sparse_module = (:sp_fro_norm, :sp_inf_norm, :sp_L2_norm, :sp_max_norm, :sp_one_norm, :sp_reshape, :sp_eigs, :sp_expm, :sp_permute, :sp_reverse_permute, :sp_bandwidth, :sp_profile)

# states.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/states.py
export basis, qutrit_basis, coherent, coherent_dm, fock_dm, fock, thermal_dm, maximally_mixed_dm, ket2dm, projection, qstate, ket, bra, state_number_enumerate, state_number_index, state_index_number, state_number_qobj, phase_basis, zero_ket, spin_state, spin_coherent, bell_state, singlet_state, triplet_states, w_state, ghz_state, enr_state_dictionaries, enr_fock, enr_thermal_dm
const states_module = (:basis, :qutrit_basis, :coherent, :coherent_dm, :fock_dm, :fock, :thermal_dm, :maximally_mixed_dm, :ket2dm, :projection, :qstate, :ket, :bra, :state_number_enumerate, :state_number_index, :state_index_number, :state_number_qobj, :phase_basis, :zero_ket, :spin_state, :spin_coherent, :bell_state, :singlet_state, :triplet_states, :w_state, :ghz_state, :enr_state_dictionaries, :enr_fock, :enr_thermal_dm)

# steadystate.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/steadystate.py
export steadystate, steady, build_preconditioner, pseudo_inverse
const steadystate_module = (:steadystate, :steady, :build_preconditioner, :pseudo_inverse)

# stochastic.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/stochastic.py
export ssesolve, ssepdpsolve, smesolve, smepdpsolve
const stochastic_module = (:ssesolve, :ssepdpsolve, :smesolve, :smepdpsolve)

# subsystem_apply.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/subsystem_apply.py
export subsystem_apply
const subsystem_apply_module = (:subsystem_apply, )

# superop_reps.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/superop_reps.py
export super_to_choi, choi_to_super, choi_to_kraus, kraus_to_choi, kraus_to_super, choi_to_chi, chi_to_choi, to_choi, to_chi, to_super, to_kraus, to_stinespring
const superop_reps_module = (:super_to_choi, :choi_to_super, :choi_to_kraus, :kraus_to_choi, :kraus_to_super, :choi_to_chi, :chi_to_choi, :to_choi, :to_chi, :to_super, :to_kraus, :to_stinespring)

# superoperator.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/superoperator.py
export liouvillian, liouvillian_ref, lindblad_dissipator, operator_to_vector, vector_to_operator, mat2vec, vec2mat, vec2mat_index, mat2vec_index, spost, spre, sprepost
const superoperator_module = (:liouvillian, :liouvillian_ref, :lindblad_dissipator, :operator_to_vector, :vector_to_operator, :mat2vec, :vec2mat, :vec2mat_index, :mat2vec_index, :spost, :spre, :sprepost)

# tensor.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/tensor.py
export tensor, super_tensor, composite, tensor_swap, tensor_contract
const tensor_module = (:tensor, :super_tensor, :composite, :tensor_swap, :tensor_contract)

# testing.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/testing.py
# There is no __all__.

# three_level_atom.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/three_level_atom.py
export three_level_basis, three_level_ops
const three_level_atom_module = (:three_level_basis, :three_level_ops)

# tomography.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/tomography.py
export qpt_plot, qpt_plot_combined, qpt
const tomography_module = (:qpt_plot, :qpt_plot_combined, :qpt)

# utilities.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/utilities.py
export n_thermal, linspace_with, clebsch, convert_unit, view_methods
const utilities_module = (:n_thermal, :linspace_with, :clebsch, :convert_unit, :view_methods)

# visualization.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/visualization.py
export hinton, sphereplot, energy_level_diagram, plot_fock_distribution, wigner_fock_distribution, plot_wigner_fock_distribution, plot_wigner, plot_expectation_values, plot_spin_distribution_2d, plot_spin_distribution_3d, plot_qubism, plot_schmidt, complex_array_to_rgb, matrix_histogram, matrix_histogram_complex, sphereplot
const visualization_module = (:hinton, :sphereplot, :energy_level_diagram, :plot_fock_distribution, :wigner_fock_distribution, :plot_wigner_fock_distribution, :plot_wigner, :plot_expectation_values, :plot_spin_distribution_2d, :plot_spin_distribution_3d, :plot_qubism, :plot_schmidt, :complex_array_to_rgb, :matrix_histogram, :matrix_histogram_complex, :sphereplot)

# wigner.py https://github.com/qutip/qutip/blob/v4.2.0/qutip/wigner.py
# wigner is defined individually.
export qfunc, spin_q_function, spin_wigner
const wigner_module = (:qfunc, :spin_q_function, :spin_wigner)


###############################################################
# Function
###############################################################
const qutipfn = (
                gate_module...,
                about_module...,
                bloch_module...,
                bloch3d_module...,
                bloch_redfield_module...,
                continuous_variables_module...,
                correlation_module...,
                countstat_module...,
                distributions_module...,
                # entropy_module...,
                eseries_module...,
                essolve_module...,
                expect_module...,
                fileio_module...,
                floquet_module...,
                graph_module...,
                hardware_info_module...,
                interpolate_module...,
                # ipynbtools_module...,
                logging_utils_module...,
                matplotlib_utilities_module...,
                mcsolve_module...,
                mesolve_module...,
                # metrics_module...,
                operators_module...,
                orbital_module...,
                parallel_module...,
                partial_transpose_module...,
                permute_module...,
                propagator_module...,
                qobj_module...,
                random_objects_module...,
                rcsolve_module...,
                rhs_generate_module...,
                sesolve_module...,
                simdiag_module...,
                solver_module...,
                sparse_module...,
                states_module...,
                steadystate_module...,
                stochastic_module...,
                subsystem_apply_module...,
                superop_reps_module...,
                superoperator_module...,
                tensor_module...,
                three_level_atom_module...,
                tomography_module...,
                # utilities_module...,
                # visualization_module...,
                wigner_module...
               )

for f in qutipfn
    sf = string(f)
    @eval @doc LazyHelp(qutip,$sf) function $f(args...; kws...)
        if !haskey(qutip, $sf)
            error("qutip ", version, " does not have qutip.", $sf)
        end
        return pycall(qutip[$sf], Quantum, args...; kws...)
    end
end

# for f in ipynbtools_module
#     sf = string(f)
#     @eval @doc LazyHelp(ipynbtools,$sf) function $f(args...; kws...)
#         if !haskey(ipynbtools, $sf)
#             error("qutip ", version, " does not have qutip.ipynbtools", $sf)
#         end
#         return pycall(ipynbtools[$sf], PyAny, args...; kws...)
#     end
# end

for f in metrics_module
    sf = string(f)
    @eval @doc LazyHelp(qutip,$sf) function $f(args...; kws...)
        if !haskey(qutip, $sf)
            error("qutip ", version, " does not have qutip.", $sf)
        end
        return pycall(qutip[$sf], PyAny, args...; kws...)
    end
end

for f in entropy_module
    sf = string(f)
    @eval @doc LazyHelp(qutip,$sf) function $f(args...; kws...)
        if !haskey(qutip, $sf)
            error("qutip ", version, " does not have qutip.", $sf)
        end
        return pycall(qutip[$sf], PyAny, args...; kws...)
    end
end

for f in visualization_module
    sf = string(f)
    @eval @doc LazyHelp(visualization,$sf) function $f(args...; kws...)
        if !haskey(visualization, $sf)
            error("qutip ", version, " does not have qutip.visualization", $sf)
        end
        return pycall(visualization[$sf], PyAny, args...; kws...)
    end
end

for f in utilities_module
    sf = string(f)
    @eval @doc LazyHelp(qutip,$sf) function $f(args...; kws...)
        if !haskey(qutip, $sf)
            error("qutip ", version, " does not have qutip.utilities. ", $sf)
        end
        return pycall(qutip[$sf], PyAny, args...; kws...)
    end
end


# Functions whose type of return value is not Qobj.
export expect
export esspec, esval
export essolve
export wigner
for f in (:expect, :esspec, :esval, :essolve, :wigner,
         correlation_module...)
    sf = string(f)
    @eval @doc LazyHelp(qutip,$sf) function $f(args...; kws...)
        if !haskey(qutip, $sf)
            error("qutip ", version, " does not have qutip.", $sf)
        end
        return pycall(qutip[$sf], PyAny, args...; kws...)
    end
end

export bloch_redfield_tensor
f = (:bloch_redfield_tensor)
sf = string(f)
@eval @doc LazyHelp(qutip,$sf) function $f(args...; kws...)
    if !haskey(qutip, $sf)
        error("qutip ", version, " does not have qutip.", $sf)
    end
    return pycall(qutip[$sf], Tuple{Quantum, Vector{Quantum}}, args...; kws...)
end

#############################
# Nonmarkov module
#############################
module Nonmarkov
using ..QuTiP
import PyCall: PyNULL, pyimport, pycall
import QuTiP: LazyHelp, Quantum, qutip, nonmarkov, transfertensor, memorycascade

# nonmarkov/memorycascade.py and nonmarkov/transfertensor.py
export ttmsolve, MemoryCascade

f = :ttmsolve
sf = string(f)
@eval @doc LazyHelp(transfertensor, $sf) function $f(args...; kws...)
    if !haskey(transfertensor, $sf)
        error("qutip.nonmarkov.transfertensor ", version, " does not have qutip.nonmarkov.transfertensor", $sf)
    end
    return pycall(transfertensor[$sf], Quantum, args...; kws...)
end

f = :MemoryCascade
sf = string(f)
@eval @doc LazyHelp(memorycascade,$sf) function $f(args...; kws...)
    if !haskey(memorycascade, $sf)
        error("qutip.nonmarkov.memorycascade ", version, " does not have qutip.nonmarkov.memorycascade", $sf)
    end
    return pycall(memorycascade[$sf], Quantum, args...; kws...)
end

end # module nonmarkov

