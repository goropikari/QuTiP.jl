export data, dims, shape, superrep, isherm, iscp, ishp, istp, iscptp, isket, isbra, isoper, issuper, isoperket, isoperbra
const attributes = (:data, :dims, :shape, :type, :superrep, :isherm, :iscp, :ishp, :istp, :iscptp, :isket, :isbra, :isoper, :issuper, :isoperket, :isoperbra)

export dnorm, dual_chan, tidyup, trans, transform, trunc_neg, unit, eliminate_states, evaluate, extract_states
const methods_qobj =  (:dnorm, :dual_chan, :tidyup, :trans, :transform, :trunc_neg, :unit, :eliminate_states, :evaluate, :extract_states)

export full, norm, permute, eigenenergies, eigenstates, groundstate, matrix_element, tr, diag, overlap
const methods = (:full, :norm, :permute, :eigenenergies, :eigenstates, :groundstate, :matrix_element, :tr, :diag, :overlap)

# renamed
# type 
