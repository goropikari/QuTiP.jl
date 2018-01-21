import Base: conj, expm, sqrtm, full, norm, diag

###################
# Qobj
###################
export data, dims, shape, superrep, isherm, iscp, ishp, istp, iscptp, isket, isbra, isoper, issuper, isoperket, isoperbra # Qob
export rates # eseries
const attributes = (:data, :dims, :shape, :type, :superrep, :isherm, :iscp, :ishp, :istp, :iscptp, :isket, :isbra, :isoper, :issuper, :isoperket, :isoperbra, # Qobj
:rates, # eseries
)

export dual_chan, tidyup, trans, transform, trunc_neg, unit, eliminate_states, evaluate, extract_states # Qob
const methods_qobj =  (:dual_chan, :tidyup, :trans, :transform, :trunc_neg, :unit, :eliminate_states, :evaluate, :extract_states, # QObj
)

export full, norm, dnorm, eigenenergies, matrix_element, tr, diag, overlap
const methods = (:full, :norm, :dnorm, :eigenenergies, :matrix_element, :tr, :diag, :overlap)

# renamed
# type

##########################
# individualy define
##########################
#=
    eigenstates, groundstate
    permute

    ampl
    value, spec
=#


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
