export data, dims, shape, superrep, isherm, iscp, ishp, istp, iscptp, isket, isbra, isoper, issuper, isoperket, isoperbra
const attributes = (:data, :dims, :shape, :type, :superrep, :isherm, :iscp, :ishp, :istp, :iscptp, :isket, :isbra, :isoper, :issuper, :isoperket, :isoperbra)

export dual_chan, tidyup, trans, transform, trunc_neg, unit, eliminate_states, evaluate, extract_states
const methods_qobj =  (:dual_chan, :tidyup, :trans, :transform, :trunc_neg, :unit, :eliminate_states, :evaluate, :extract_states)

export full, norm, dnorm, eigenenergies, matrix_element, tr, diag, overlap
const methods = (:full, :norm, :dnorm, :eigenenergies, :matrix_element, :tr, :diag, :overlap)

# renamed
# type

# individualy define
# eigenstates, groundstate
# permute
