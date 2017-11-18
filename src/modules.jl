# From qutip 4.2.0

# about.py
# __all__ = [about]
## export about
about_module = (:about, )

# bloch.py
# __all__ = [Bloch]
export Bloch
bloch_module = (:Bloch, )

# bloch3d.py
# __all__ = [Bloch3d]
export Bloch3d
bloch3d_module = (:Bloch3d, )

# bloch_redfield.py
# __all__ = [brmesolve, bloch_redfield_solve, bloch_redfield_tensor]
# bloch_redfield_tensor is defined individually.
export brmesolve, bloch_redfield_solve
bloch_redfield_module = (:brmesolve, :bloch_redfield_solve)

# configrc.py
# There is no __all__.

# continuous_variables.py
# __all__ = [correlation_matrix, covariance_matrix, correlation_matrix_field, correlation_matrix_quadrature, wigner_covariance_matrix, logarithmic_negativity]
export correlation_matrix, covariance_matrix,  correlation_matrix_field, correlation_matrix_quadrature, wigner_covariance_matrix, logarithmic_negativity
const continuous_variables_module = (:correlation_matrix, :covariance_matrix, :correlation_matrix_field, :correlation_matrix_quadrature, :wigner_covariance_matrix, :logarithmic_negativity)

# correlation.py
# __all__ = [correlation_2op_1t, correlation_2op_2t, correlation_3op_1t, correlation_3op_2t, coherence_function_g1, coherence_function_g2, spectrum, spectrum_correlation_fft, correlation_ss, correlation, correlation_4op_1t, correlation_4op_2t, spectrum_ss, spectrum_pi]
export correlation_2op_1t, correlation_2op_2t, correlation_3op_1t, correlation_3op_2t, coherence_function_g1, coherence_function_g2, spectrum, spectrum_correlation_fft, correlation_ss, correlation, correlation_4op_1t, correlation_4op_2t, spectrum_ss, spectrum_pi
correlation_module = (:correlation_2op_1t, :correlation_2op_2t, :correlation_3op_1t, :correlation_3op_2t, :coherence_function_g1, :coherence_function_g2, :spectrum, :spectrum_correlation_fft, :correlation_ss, :correlation, :correlation_4op_1t, :correlation_4op_2t, :spectrum_ss, :spectrum_pi)

# countstat.py
# __all__ = [countstat_current, countstat_current_noise]
export countstat_current, countstat_current_noise
countstat_module = (:countstat_current, :countstat_current_noise)

# dimensions.py
# __all__ = [] # Everything should be explicitly imported, not made available
#              # by default.

# distributions.py
# __all__ = [Distribution, WignerDistribution, QDistribution, TwoModeQuadratureCorrelation, HarmonicOscillatorWaveFunction, HarmonicOscillatorProbabilityFunction]
export Distribution, WignerDistribution, QDistribution, TwoModeQuadratureCorrelation, HarmonicOscillatorWaveFunction, HarmonicOscillatorProbabilityFunction
distributions_module = (:Distribution, :WignerDistribution, :QDistribution, :TwoModeQuadratureCorrelation, :HarmonicOscillatorWaveFunction, :HarmonicOscillatorProbabilityFunction)

# entropy.py
# __all__ = [entropy_vn, entropy_linear, entropy_mutual, negativity, concurrence, entropy_conditional, entangling_power]
export entropy_vn, entropy_linear, entropy_mutual, negativity, concurrence, entropy_conditional, entangling_power
entropy_module = (:entropy_vn, :entropy_linear, :entropy_mutual, :negativity, :concurrence, :entropy_conditional, :entangling_power)

# eseries.py
# __all__ = [eseries, esval, esspec, estidy]
# esspec and esval are defined individually.
export eseries, estidy
eseries_module = (:eseries, :estidy)

# essolve.py
# __all__ = [essolve, ode2es]
# essolve is defined individually.
export ode2es
essolve_module = (:ode2es, )

# expect.py
# __all__ = [expect, variance]
# expect is defined individually.
export variance
expect_module = (:variance, )

# fastsparse.py
# There is no __all__.

# fileio.py
# __all__ = [file_data_store, file_data_read, qsave, qload]
export file_data_store, file_data_read, qsave, qload
fileio_module = (:file_data_store, :file_data_read, :qsave, :qload)

# floquet.py
# __all__ = [floquet_modes, floquet_modes_t, floquet_modes_table, floquet_modes_t_lookup, floquet_states, floquet_states_t, floquet_wavefunction, floquet_wavefunction_t, floquet_state_decomposition, fsesolve, floquet_master_equation_rates, floquet_collapse_operators, floquet_master_equation_tensor, floquet_master_equation_steadystate, floquet_basis_transform, floquet_markov_mesolve, fmmesolve]
export floquet_modes, floquet_modes_t, floquet_modes_table, floquet_modes_t_lookup, floquet_states, floquet_states_t, floquet_wavefunction, floquet_wavefunction_t, floquet_state_decomposition, fsesolve, floquet_master_equation_rates, floquet_collapse_operators, floquet_master_equation_tensor, floquet_master_equation_steadystate, floquet_basis_transform, floquet_markov_mesolve, fmmesolve
floquet_module = (:floquet_modes, :floquet_modes_t, :floquet_modes_table, :floquet_modes_t_lookup, :floquet_states, :floquet_states_t, :floquet_wavefunction, :floquet_wavefunction_t, :floquet_state_decomposition, :fsesolve, :floquet_master_equation_rates, :floquet_collapse_operators, :floquet_master_equation_tensor, :floquet_master_equation_steadystate, :floquet_basis_transform, :floquet_markov_mesolve, :fmmesolve)

# graph.py
# __all__ = [graph_degree, column_permutation, breadth_first_search,reverse_cuthill_mckee, maximum_bipartite_matching, weighted_bipartite_matching]
export graph_degree, column_permutation, breadth_first_search,reverse_cuthill_mckee, maximum_bipartite_matching, weighted_bipartite_matching
graph_module = (:graph_degree, :column_permutation, :breadth_first_search, :reverse_cuthill_mckee, :maximum_bipartite_matching, :weighted_bipartite_matching)

# hardware_info.py
# __all__ = [hardware_info]
hardware_info_module = (:hardware_info, )

# interpolate.py
# __all__ = [Cubic_Spline]
export Cubic_Spline
interpolate_module = (:Cubic_Spline, )

# ipynbtools.py
#            if IPython.version_info[0] >= 4:
#                try:
#                    from ipyparallel import Client
#                    __all__ = [version_table, parfor, plot_animation,
#                                parallel_map, HTMLProgressBar]
#                except:
#                     __all__ = [version_table, plot_animation, HTMLProgressBar]
#            else:
#                try:
#                    from IPython.parallel import Client
#                    __all__ = [version_table, parfor, plot_animation,
#                                parallel_map, HTMLProgressBar]
#                except:
#                     __all__ = [version_table, plot_animation, HTMLProgressBar]
# export version_table, plot_animation, HTMLProgressBar
# ipynbtools_module = (:version_table, :plot_animation, :HTMLProgressBar)

# logging_utils.py
# __all__ = [get_logger]
export get_logger
logging_utils_module = (:get_logger, )

# matplotlib_utilities.py
# __all__ = [wigner_cmap, MidpointNorm, complex_phase_cmap]
export wigner_cmap, MidpointNorm, complex_phase_cmap
matplotlib_utilities_module = (:wigner_cmap, :MidpointNorm, :complex_phase_cmap)

# mcsolve.py
# __all__ = [mcsolve]
export mcsolve
mcsolve_module = (:mcsolve, )

# mesolve.py
# __all__ = [mesolve, odesolve]
export mesolve, odesolve
mesolve_module = (:mesolve, :odesolve)

# metrics.py
# __all__ = [fidelity, tracedist, bures_dist, bures_angle, hilbert_dist, average_gate_fidelity, process_fidelity, unitarity, dnorm]
export fidelity, tracedist, bures_dist, bures_angle, hilbert_dist, average_gate_fidelity, process_fidelity, unitarity, dnorm
metrics_module = (:fidelity, :tracedist, :bures_dist, :bures_angle, :hilbert_dist, :average_gate_fidelity, :process_fidelity, :unitarity, :dnorm)

# operators.py
# __all__ = [jmat, spin_Jx, spin_Jy, spin_Jz, spin_Jm, spin_Jp, spin_J_set, sigmap, sigmam, sigmax, sigmay, sigmaz, destroy, create, qeye, identity, position, momentum, num, squeeze, squeezing, displace, commutator, qutrit_ops, qdiags, phase, qzero, enr_destroy, enr_identity, charge, tunneling]
# identity, position, num and squeeze are renamed in order to avoid name conflict with Base module functions.
export jmat, spin_Jx, spin_Jy, spin_Jz, spin_Jm, spin_Jp, spin_J_set, sigmap, sigmam, sigmax, sigmay, sigmaz, destroy, create, qeye, momentum, squeezing, displace, commutator, qutrit_ops, qdiags, phase, qzero, enr_destroy, enr_identity, charge, tunneling
operators_module = (:jmat, :spin_Jx, :spin_Jy, :spin_Jz, :spin_Jm, :spin_Jp, :spin_J_set, :sigmap, :sigmam, :sigmax, :sigmay, :sigmaz, :destroy, :create, :qeye, :identity, :position, :momentum, :num, :squeeze, :squeezing, :displace, :commutator, :qutrit_ops, :qdiags, :phase, :qzero, :enr_destroy, :enr_identity, :charge, :tunneling)

# orbital.py
# __all__ = [orbital]
export orbital
orbital_module = (:orbital, )

# parallel.py
# __all__ = [parfor, parallel_map, serial_map]
export parfor, parallel_map, serial_map
parallel_module = (:parfor, :parallel_map, :serial_map)

# partial_transpose.py
# __all__ = [partial_transpose]
export partial_transpose
partial_transpose_module = (:partial_transpose, )

# permute.py
# __all__ = [reshuffle]
export reshuffle
permute_module = (:reshuffle, )

# propagator.py
# __all__ = [propagator, propagator_steadystate]
export propagator, propagator_steadystate
propagator_module = (:propagator, :propagator_steadystate)

# qobj.py
# __all__ = [Qobj, qobj_list_evaluate, ptrace, dag, isequal, issuper, isoper, isoperket, isoperbra, isket, isbra, isherm, shape, dims]
export Qobj, qobj_list_evaluate, ptrace, dag, isequal, issuper, isoper, isoperket, isoperbra, isket, isbra, isherm, shape, dims
qobj_module = (:Qobj, :qobj_list_evaluate, :ptrace, :dag, :isequal, :issuper, :isoper, :isoperket, :isoperbra, :isket, :isbra, :isherm, :shape, :dims)

# random_objects.py
# __all__ = [rand_herm, rand_unitary, rand_ket, rand_dm, rand_unitary_haar, rand_ket_haar, rand_dm_ginibre, rand_dm_hs, rand_super_bcsz, rand_stochastic, rand_super]
export rand_herm, rand_unitary, rand_ket, rand_dm, rand_unitary_haar, rand_ket_haar, rand_dm_ginibre, rand_dm_hs, rand_super_bcsz, rand_stochastic, rand_super
random_objects_module = (:rand_herm, :rand_unitary, :rand_ket, :rand_dm, :rand_unitary_haar, :rand_ket_haar, :rand_dm_ginibre, :rand_dm_hs, :rand_super_bcsz, :rand_stochastic, :rand_super)

# rcsolve.py
# __all__ = [rcsolve]
export rcsolve
rcsolve_module = (:rcsolve, )

# rhs_generate.py
# __all__ = [rhs_generate, rhs_clear]
export rhs_generate, rhs_clear
rhs_generate_module = (:rhs_generate, :rhs_clear)

# semidefinite.py
# There is no __all__.

# sesolve.py
# __all__ = [sesolve]
export sesolve
sesolve_module = (:sesolve, )

# settings.py
# There is no __all__.

# simdiag.py
# __all__ = [simdiag]
export simdiag
simdiag_module = (:simdiag, )

# solver.py
# __all__ = [Options, Odeoptions, Odedata]
export Options, Odeoptions, Odedata
solver_module = (:Options, :Odeoptions, :Odedata)

# sparse.py
# __all__ = [sp_fro_norm, sp_inf_norm, sp_L2_norm, sp_max_norm, sp_one_norm, sp_reshape, sp_eigs, sp_expm, sp_permute, sp_reverse_permute, sp_bandwidth, sp_profile]
export sp_fro_norm, sp_inf_norm, sp_L2_norm, sp_max_norm, sp_one_norm, sp_reshape, sp_eigs, sp_expm, sp_permute, sp_reverse_permute, sp_bandwidth, sp_profile
sparse_module = (:sp_fro_norm, :sp_inf_norm, :sp_L2_norm, :sp_max_norm, :sp_one_norm, :sp_reshape, :sp_eigs, :sp_expm, :sp_permute, :sp_reverse_permute, :sp_bandwidth, :sp_profile)

# states.py
# __all__ = [basis, qutrit_basis, coherent, coherent_dm, fock_dm, fock, thermal_dm, maximally_mixed_dm, ket2dm, projection, qstate, ket, bra, state_number_enumerate, state_number_index, state_index_number, state_number_qobj, phase_basis, zero_ket, spin_state, spin_coherent, bell_state, singlet_state, triplet_states, w_state, ghz_state, enr_state_dictionaries, enr_fock, enr_thermal_dm]
export basis, qutrit_basis, coherent, coherent_dm, fock_dm, fock, thermal_dm, maximally_mixed_dm, ket2dm, projection, qstate, ket, bra, state_number_enumerate, state_number_index, state_index_number, state_number_qobj, phase_basis, zero_ket, spin_state, spin_coherent, bell_state, singlet_state, triplet_states, w_state, ghz_state, enr_state_dictionaries, enr_fock, enr_thermal_dm
states_module = (:basis, :qutrit_basis, :coherent, :coherent_dm, :fock_dm, :fock, :thermal_dm, :maximally_mixed_dm, :ket2dm, :projection, :qstate, :ket, :bra, :state_number_enumerate, :state_number_index, :state_index_number, :state_number_qobj, :phase_basis, :zero_ket, :spin_state, :spin_coherent, :bell_state, :singlet_state, :triplet_states, :w_state, :ghz_state, :enr_state_dictionaries, :enr_fock, :enr_thermal_dm)

# steadystate.py
# __all__ = [steadystate, steady, build_preconditioner, pseudo_inverse]
export steadystate, steady, build_preconditioner, pseudo_inverse
steadystate_module = (:steadystate, :steady, :build_preconditioner, :pseudo_inverse)

# stochastic.py
# __all__ = [ssesolve, ssepdpsolve, smesolve, smepdpsolve]
export ssesolve, ssepdpsolve, smesolve, smepdpsolve
stochastic_module = (:ssesolve, :ssepdpsolve, :smesolve, :smepdpsolve)

# subsystem_apply.py
# __all__ = [subsystem_apply]
export subsystem_apply
subsystem_apply_module = (:subsystem_apply, )

# superop_reps.py
# __all__ = [super_to_choi, choi_to_super, choi_to_kraus, kraus_to_choi, kraus_to_super, choi_to_chi, chi_to_choi, to_choi, to_chi, to_super, to_kraus, to_stinespring]
export super_to_choi, choi_to_super, choi_to_kraus, kraus_to_choi, kraus_to_super, choi_to_chi, chi_to_choi, to_choi, to_chi, to_super, to_kraus, to_stinespring
superop_reps_module = (:super_to_choi, :choi_to_super, :choi_to_kraus, :kraus_to_choi, :kraus_to_super, :choi_to_chi, :chi_to_choi, :to_choi, :to_chi, :to_super, :to_kraus, :to_stinespring)

# superoperator.py
# __all__ = [liouvillian, liouvillian_ref, lindblad_dissipator, operator_to_vector, vector_to_operator, mat2vec, vec2mat, vec2mat_index, mat2vec_index, spost, spre, sprepost]
export liouvillian, liouvillian_ref, lindblad_dissipator, operator_to_vector, vector_to_operator, mat2vec, vec2mat, vec2mat_index, mat2vec_index, spost, spre, sprepost
superoperator_module = (:liouvillian, :liouvillian_ref, :lindblad_dissipator, :operator_to_vector, :vector_to_operator, :mat2vec, :vec2mat, :vec2mat_index, :mat2vec_index, :spost, :spre, :sprepost)

# tensor.py
# __all__ = [tensor, super_tensor, composite, tensor_swap, tensor_contract]
export tensor, super_tensor, composite, tensor_swap, tensor_contract
tensor_module = (:tensor, :super_tensor, :composite, :tensor_swap, :tensor_contract)

# testing.py
# There is no __all__.

# three_level_atom.py
# __all__ = [three_level_basis, three_level_ops]
export three_level_basis, three_level_ops
three_level_atom_module = (:three_level_basis, :three_level_ops)

# tomography.py
# __all__ = [qpt_plot, qpt_plot_combined, qpt]
export qpt_plot, qpt_plot_combined, qpt
tomography_module = (:qpt_plot, :qpt_plot_combined, :qpt)

# utilities.py
# __all__ = [n_thermal, linspace_with, clebsch, convert_unit, view_methods]
export n_thermal, linspace_with, clebsch, convert_unit, view_methods
utilities_module = (:n_thermal, :linspace_with, :clebsch, :convert_unit, :view_methods)

# visualization.py
# __all__ = [hinton, sphereplot, energy_level_diagram, plot_fock_distribution, wigner_fock_distribution, plot_wigner_fock_distribution, plot_wigner, plot_expectation_values, plot_spin_distribution_2d, plot_spin_distribution_3d, plot_qubism, plot_schmidt, complex_array_to_rgb, matrix_histogram, matrix_histogram_complex, sphereplot]
export hinton, sphereplot, energy_level_diagram, plot_fock_distribution, wigner_fock_distribution, plot_wigner_fock_distribution, plot_wigner, plot_expectation_values, plot_spin_distribution_2d, plot_spin_distribution_3d, plot_qubism, plot_schmidt, complex_array_to_rgb, matrix_histogram, matrix_histogram_complex, sphereplot
visualization_module = (:hinton, :sphereplot, :energy_level_diagram, :plot_fock_distribution, :wigner_fock_distribution, :plot_wigner_fock_distribution, :plot_wigner, :plot_expectation_values, :plot_spin_distribution_2d, :plot_spin_distribution_3d, :plot_qubism, :plot_schmidt, :complex_array_to_rgb, :matrix_histogram, :matrix_histogram_complex, :sphereplot)

# wigner.py
# __all__ = ['wigner', 'qfunc', 'spin_q_function', 'spin_wigner']
# wigner is defined individually.
export qfunc, spin_q_function, spin_wigner
wigner_module = (:qfunc, :spin_q_function, :spin_wigner)


###############################################################
# Function
###############################################################
const qutipfn = (
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
#     @eval @doc LazyHelp(qutip[:ipynbtools],$sf) function $f(args...; kws...)
#         if !haskey(ipynbtools, $sf)
#             error("qutip ", version, " does not have qutip.ipynbtools", $sf)
#         end
#         return pycall(qutip[:ipynbtools][$sf], PyAny, args...; kws...)
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
