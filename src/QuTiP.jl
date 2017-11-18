__precompile__()

module QuTiP
using PyCall
import PyCall: PyNULL, pyimport_conda, pycall, PyObject
import Base: +, -, *, /, ==, hash, getindex, setindex!, haskey, keys, show, convert, collect
import Base: conj, expm, sqrtm, full, norm, diag, ctranspose

export qutip, Quantum

import Base.show

type Quantum
    o::PyObject
end

include("display.jl")

###########################################################################
# quoted from PyPlot.jl https://github.com/JuliaPy/PyPlot.jl/blob/2476177334dd0640711c2ca5d2c5c6bb0dc3c317/src/PyPlot.jl#L18
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
const visualization = PyNULL()
function __init__()
    pyimport_conda("IPython", "IPython")
    pyimport_conda("matplotlib", "matplotlib")
    copy!(qutip, pyimport_conda("qutip", "qutip", "conda-forge"))
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

###############################################################################
# export ducumented qutip API
###############################################################################
include("modules.jl")
include("attributes_methods.jl")

    ############################################################################
    # In order to avoid name conflict with Base module functions, add prefix 'q'.
    ############################################################################
    export qidentity, qposition, qnum, qsqueeze# operators module
    const renamedfn = (:identity, :position, :num, :squeeze)
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

# define functions for convenience
export ⊗, ctranspose
⊗(a::Quantum, b::Quantum) = tensor(a,b)
ctranspose(x::Quantum) = dag(x::Quantum)

###################################################
# arithmetic
#
# Why do we define arithmetic?
# Example
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
