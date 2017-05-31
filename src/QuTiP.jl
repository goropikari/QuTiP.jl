__precompile__()

module QuTiP 
using PyCall
import PyCall: PyNULL, pyimport_conda, pycall, PyObject
import Base: +, -, *, /, ==, hash, getindex, setindex!, haskey, keys, show, convert
import Base: conj, expm, sqrtm

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
const utilities = PyNULL()

function __init__()
    pyimport_conda("IPython", "IPython")
    pyimport_conda("matplotlib", "matplotlib")
    copy!(qutip, pyimport_conda("qutip", "qutip", "conda-forge"))
    copy!(ipynbtools, pyimport("qutip.ipynbtools"))
    copy!(visualization, pyimport("qutip.visualization"))
    copy!(utilities, pyimport("qutip.utilities"))
    global const version = try
        convert(VersionNumber, qutip[:__version__])
    catch
        v"0.0" # fallback
    end
end


PyObject(f::Quantum) = f.o
convert(::Type{Quantum}, o::PyObject) = Quantum(o)
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
include("class/utilities.jl")
include("class/sparse.jl")
include("class/simdiag.jl")
include("class/permute.jl")
include("class/parallel.jl")
include("class/ipynbtools.jl")
include("class/hardware_info.jl")
include("class/graph.jl")
include("class/fileio.jl")
include("class/about.jl")

include("class/tensor.jl")
include("class/qobj.jl")
include("class/partial_transpose.jl")
include("class/expect.jl")

include("class/metrics.jl")
include("class/entropy.jl")
include("class/countstat.jl")

include("class/three_level_atom.jl")
include("class/states.jl")
include("class/random_objects.jl")
include("class/continuous_variables.jl")
include("class/superoperator.jl")
include("class/superop_reps.jl")
include("class/subsystem_apply.jl")
include("class/operators.jl")

include("class/bloch_redfield.jl")
include("class/correlation.jl")
include("class/eseries.jl")
include("class/essolve.jl")
include("class/floquet.jl")
include("class/hsolve.jl")
include("class/mcsolve.jl")
include("class/mesolve.jl")
include("class/propagator.jl")
include("class/rcsolve.jl")
include("class/rhs_generate.jl")
include("class/sesolve.jl")
include("class/solver.jl")
include("class/steadystate.jl")
include("class/stochastic.jl")
include("class/memorycascade.jl")
include("class/transfertensor.jl")

include("class/settings.jl")

include("class/bloch.jl")
include("class/bloch3d.jl")
include("class/distributions.jl")
include("class/orbital.jl")
include("class/tomography.jl")
include("class/visualization.jl")
include("class/wigner.jl")

include("class/gate.jl")

include("attributes_methods.jl")

###############################################################
# Function
###############################################################
const qutipfn = (#utilities_class...,
                sparse_class...,
                simdiag_class...,
                permute_class...,
                parallel_class...,
                # ipynbtools_class...,
                hardware_info_class...,
                graph_class...,
                fileio_class...,
                about_class...,
                tensor_class..., 
                qobj_class...,
                partial_transpose_class...,
                expect_class...,
                metrics_class...,
                entropy_class..., 
                countstat_class..., 
                three_level_atom_class...,
                states_class..., 
                random_objects_class...,
                continuous_variables_class...,
                superoperator_class..., 
                superop_reps_class..., 
                subsystem_apply_class..., 
                operators_class...,
                bloch_redfield_class...,
                # correlation_class...,
                eseries_class...,
                essolve_class...,
                floquet_class...,
                hsolve_class...,
                mcsolve_class...,
                mesolve_class...,
                propagator_class...,
                rcsolve_class...,
                rhs_generate_class...,
                sesolve_class...,
                solver_class...,
                steadystate_class...,
                stochastic_class...,
                memorycascade_class...,
                transfertensor_class...,
                settings_class...,
                bloch_class...,
                bloch3d_class...,
                distributions_class...,
                orbital_class...,
                tomography_class...,
                # visualization_class...,
                wigner_class...,
                gate_class...
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

for f in ipynbtools_class
    sf = string(f)
    @eval @doc LazyHelp(ipynbtools,$sf) function $f(args...; kws...)
        if !haskey(ipynbtools, $sf)
            error("qutip.ipynbtools ", version, " does not have qutip.ipynbtools", $sf)
        end
        return pycall(ipynbtools[$sf], PyAny, args...; kws...)
    end
end

for f in visualization_class
    sf = string(f)
    @eval @doc LazyHelp(visualization,$sf) function $f(args...; kws...)
        if !haskey(visualization, $sf)
            error("qutip.visualization ", version, " does not have qutip.visualization", $sf)
        end
        return pycall(visualization[$sf], PyAny, args...; kws...)
    end
end

for f in utilities_class
    sf = string(f)
    @eval @doc LazyHelp(utilities,$sf) function $f(args...; kws...)
        if !haskey(utilities, $sf)
            error("qutip.utilities ", version, " does not have qutip.utilities", $sf)
        end
        return pycall(utilities[$sf], PyAny, args...; kws...)
    end
end

# Functions whose type of return value is not Qobj.
export  expect
export esspec, esval
export  essolve
for f in (:expect, :esspec, :esval, :essolve, 
         correlation_class..., )
    sf = string(f)
    @eval @doc LazyHelp(qutip,$sf) function $f(args...; kws...)
        if !haskey(qutip, $sf)
            error("qutip ", version, " does not have qutip.", $sf)
        end
        return pycall(qutip[$sf], PyAny, args...; kws...)
    end
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



#######################################################################
# To avoid name conflict with Base module functions, add prefix 'q'.
#######################################################################
export qidentity, qnum, qposition, qsqueeze # operators class
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

export qfull, qnorm, qpermute # methdos
const renamedmethods = (:full, :norm, :permute,)
for m in renamedmethods
    sm = string(m)
    nm = Symbol("q", m)
    @eval function $nm(x, args...; kws...)
        if !haskey(x, $sm)
            error("KeyError: key $sm not found")
        end
        return x[$sm](args...; kws...)
    end
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
-(a::Number, b::Quantum) = convert(Quantum, - PyObject(b) + a)
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


end # module

