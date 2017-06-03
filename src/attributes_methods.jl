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
