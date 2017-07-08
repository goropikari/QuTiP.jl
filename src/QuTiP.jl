__precompile__()

module QuTiP
using PyCall
import PyCall: PyNULL, pyimport_conda, pycall, PyObject
import Base: +, -, *, /, ==, hash, getindex, setindex!, haskey, keys, show, convert, collect
import Base: conj, expm, sqrtm, full, norm, diag

export qutip

import Base.show

type Quantum
    o::PyObject
end


###########################################################################
# quoted from PyPlot.jl
# Julia 0.4 help system: define a documentation object
# that lazily looks up help from a PyObject via zero or more keys.
# This saves us time when loading PyPlot, since we don't have
# to load up all of the documentation strings right away.
immutable LazyHelp
    o::PyObject
    keys::Tuple{Vararg{String}}
    LazyHelp(o::PyObject) = new(o, ())
    LazyHelp(o::PyObject, k::AbstractString) = new(o, (k,))
    LazyHelp(o::PyObject, k1::AbstractString, k2::AbstractString) = new(o, (k1,k2))
    LazyHelp(o::PyObject, k::Tuple{Vararg{AbstractString}}) = new(o, k)
end
function show(io::IO, ::MIME"text/plain", h::LazyHelp)
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
Base.show(io::IO, h::LazyHelp) = show(io, "text/plain", h)
function Base.Docs.catdoc(hs::LazyHelp...)
    Base.Docs.Text() do io
        for h in hs
            show(io, MIME"text/plain"(), h)
        end
    end
end

###########################################################################

const qutip = PyNULL()
const ipynbtools = PyNULL()
const visualization = PyNULL()

function __init__()
    pyimport_conda("IPython", "IPython")
    pyimport_conda("matplotlib", "matplotlib")
    copy!(qutip, pyimport_conda("qutip", "qutip", "conda-forge"))
    copy!(ipynbtools, pyimport("qutip.ipynbtools"))
    copy!(visualization, pyimport("qutip.visualization"))
    global const version = try
        convert(VersionNumber, qutip[:__version__])
    catch
        v"0.0" # fallback
    end
end


PyObject(f::Quantum) = f.o
convert(::Type{Quantum}, o::PyObject) = Quantum(o)
collect(x::Quantum) = collect(PyObject(x))
==(f::Quantum, g::Quantum) = f.o == g.o
==(f::Quantum, g::PyObject) = f.o == g
==(f::PyObject, g::Quantum) = f == g.o
hash(f::Quantum) = hash(f.o)
pycall(f::Quantum, args...; kws...) = pycall(f.o, args...; kws...)
(f::Quantum)(args...; kws...) = pycall(f.o, Quantum, args...; kws...)

getindex(f::Quantum, x) = getindex(f.o, x)
# getindex(f::Quantum, x) = convert(Quantum, getindex(f.o, x)) # error when basis(2,0)[:isherm]
# getindex(f::Quantum, x) = try # inefficient
#         convert(Quantum, getindex(f.o, x))
#     catch
#         getindex(f.o, x)
#     end
setindex!(f::Quantum, v, x) = setindex!(f.o, v, x)
haskey(f::Quantum, x) = haskey(f.o, x)
keys(f::Quantum) = keys(f.o)



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
include("modules.jl")
include("attributes_methods.jl")

###############################################################
# Function
###############################################################
const qutipfn = (#utilities_module...,
                sparse_module...,
                simdiag_module...,
                permute_module...,
                parallel_module...,
                # ipynbtools_module...,
                hardware_info_module...,
                graph_module...,
                fileio_module...,
                about_module...,
                tensor_module...,
                qobj_module...,
                partial_transpose_module...,
                expect_module...,
                metrics_module...,
                entropy_module...,
                countstat_module...,
                three_level_atom_module...,
                states_module...,
                random_objects_module...,
                continuous_variables_module...,
                superoperator_module...,
                superop_reps_module...,
                subsystem_apply_module...,
                operators_module...,
                bloch_redfield_module...,
                # correlation_module...,
                eseries_module...,
                essolve_module...,
                floquet_module...,
                hsolve_module...,
                mcsolve_module...,
                mesolve_module...,
                propagator_module...,
                rcsolve_module...,
                rhs_generate_module...,
                sesolve_module...,
                solver_module...,
                steadystate_module...,
                stochastic_module...,
                memorycascade_module...,
                transfertensor_module...,
                settings_module...,
                bloch_module...,
                bloch3d_module...,
                distributions_module...,
                orbital_module...,
                tomography_module...,
                # visualization_module...,
                wigner_module...,
                gate_module...
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

export ⊗
⊗(a::Quantum, b::Quantum) = tensor(a,b)


for f in ipynbtools_module
    sf = string(f)
    @eval @doc LazyHelp(ipynbtools,$sf) function $f(args...; kws...)
        if !haskey(ipynbtools, $sf)
            error("qutip.ipynbtools ", version, " does not have qutip.ipynbtools", $sf)
        end
        return pycall(ipynbtools[$sf], PyAny, args...; kws...)
    end
end

for f in visualization_module
    sf = string(f)
    @eval @doc LazyHelp(visualization,$sf) function $f(args...; kws...)
        if !haskey(visualization, $sf)
            error("qutip.visualization ", version, " does not have qutip.visualization", $sf)
        end
        return pycall(visualization[$sf], PyAny, args...; kws...)
    end
end

for f in utilities_module
    sf = string(f)
    @eval @doc LazyHelp(qutip,$sf) function $f(args...; kws...)
        if !haskey(qutip, $sf)
            error("qutip.utilities ", version, " does not have qutip.utilities. ", $sf)
        end
        return pycall(qutip[$sf], PyAny, args...; kws...)
    end
end

# Functions whose type of return value is not Qobj.
export expect
export esspec, esval
export essolve
for f in (:expect, :esspec, :esval, :essolve,
         correlation_module..., )
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


###############################################################
# attributes and methods
###############################################################
for m in attributes
    sm = string(m)
    @eval function $m(x::Quantum)
        if !haskey(x, $sm)
            error("KeyError: key $sm not found")
        end
        return getindex(x, Symbol($sm))
    end
end

for m in methods_qobj
    sm = string(m)
    @eval function $m(x::Quantum, args...; kws...)
        if !haskey(x, $sm)
            error("KeyError: key $sm not found")
        end
        return convert(Quantum, x[$sm](args...; kws...))
    end
end

for m in methods
    sm = string(m)
    @eval function $m(x::Quantum, args...; kws...)
        if !haskey(x, $sm)
            error("KeyError: key $sm not found")
        end
        return x[$sm](args...; kws...)
    end
end


export ampl
function ampl(x::Quantum)
    if !haskey(x, "ampl")
        error("KeyError: key 'ampl' not found")
    end
    return convert(Vector{Quantum}, x[:ampl])
end

export conj, expm, sqrtm, sinm, cosm
for m in (:conj, :expm, :sqrtm, :sinm, :cosm)
    sm = string(m)
    @eval function $m(x::Quantum, args...; kws...)
        if !haskey(x, $sm)
            error("KeyError: key $sm not found")
        end
        return convert(Quantum, x[$sm](args...; kws...))
    end
end

export eigenstates, groundstate
function eigenstates(x::Quantum, args...; kws...)
    if !haskey(x, "eigenstates")
        error("KeyError: key 'eigenstates' not found")
    end
    return convert(Tuple{Vector{Float64},Vector{Quantum}}, x[:eigenstates](args...; kws...))
end

function groundstate(x::Quantum, args...; kws...)
    if !haskey(x, "groundstate")
        error("KeyError: key 'groundstate' not found")
    end
    return convert(Tuple{Float64, Quantum}, x[:groundstate](args...; kws...))
end

export value, spec
for m in (:value, :spec)
    sm = string(m)
    @eval function $m(x::Quantum, args...; kws...)
        if !haskey(x, $sm)
            error("KeyError: key $sm not found")
        end
        return convert(Vector{Quantum}, x[$sm](args...; kws...))
    end
end

#######################################################################
# To avoid name conflict with Base module functions, add prefix 'q'.
#######################################################################
export qidentity, qnum, qposition, qsqueeze # operators module
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

export qtype # attributes
function qtype(x::Quantum)
    sm = :type
    if !haskey(x, sm)
        error("KeyError: key $sm not found")
    end
    return x[sm]
end

export qpermute #methods
function qpermute(x::Quantum, args...; kws...)
    sm = :permute
    if !haskey(x, sm)
        error("KeyError: key $sm not found")
    end
    return convert(Quantum, x[sm](args...; kws...))
end


###################################################
# arithmetic
#
# Why I define arithmetic?
# julia> @pyimport qutip as qt
# julia> qt.sigmax() + 1
# PyObject Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
# Qobj data =
# [[ 1.  1.]
#  [ 1.  1.]]
#
# julia> 1 + qt.sigmax()
# ERROR: MethodError: no method matching +(::Int64, ::PyCall.PyObject)
# Closest candidates are:
#   +(::Any, ::Any, ::Any, ::Any...) at operators.jl:424
#   +(::PyCall.PyObject, ::Any) at /home/tk/.julia/v0.6/PyCall/src/PyCall.jl:702
#   +(::T<:Union{Int128, Int16, Int32, Int64, Int8, UInt128, UInt16, UInt32, UInt64, UInt8}, ::T<:Union{Int128, Int16, Int32, Int64, Int8, UInt128, UInt16, UInt32, UInt64, UInt8}) where T<:Union{Int128, Int16, Int32, Int64, Int8, UInt128, UInt16, UInt32, UInt64, UInt8} at int.jl:32
#   ...
###################################################

+(a::Number, b::Quantum) = convert(Quantum, PyObject(b) + a)
-(a::Number, b::Quantum) = convert(Quantum, PyObject(b) * (-1) + a)
*(a::Number, b::Quantum) = convert(Quantum, PyObject(b) * a)

+(a::Quantum, b::Number) = convert(Quantum, PyObject(a) + b)
-(a::Quantum, b::Number) = convert(Quantum, PyObject(a) - b)
*(a::Quantum, b::Number) = convert(Quantum, PyObject(a) * b)
/(a::Quantum, b::Number) = convert(Quantum, PyObject(a) / b)

+(a::PyObject, b::Quantum) = convert(Quantum, PyObject(a) + PyObject(b))
-(a::PyObject, b::Quantum) = convert(Quantum, PyObject(a) - PyObject(b))
*(a::PyObject, b::Quantum) = convert(Quantum, PyObject(a) * PyObject(b))

+(a::Quantum, b::PyObject) = convert(Quantum, PyObject(a) + PyObject(b))
-(a::Quantum, b::PyObject) = convert(Quantum, PyObject(a) - PyObject(b))
*(a::Quantum, b::PyObject) = convert(Quantum, PyObject(a) * PyObject(b))

+(a::Quantum, b::Quantum) = convert(Quantum, PyObject(a) + PyObject(b))
-(a::Quantum, b::Quantum) = convert(Quantum, PyObject(a) - PyObject(b))
*(a::Quantum, b::Quantum) = convert(Quantum, PyObject(a) * PyObject(b))

(+)(a::Quantum) = a
(-)(a::Quantum) = -1.0 * a

end # module

