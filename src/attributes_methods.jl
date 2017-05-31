export data, dims, shape, superrep, isherm, iscp, ishp, istp, iscptp, isket, isbra, isoper, issuper, isoperket, isoperbra
const attributes = (:data, :dims, :shape, :type, :superrep, :isherm, :iscp, :ishp, :istp, :iscptp, :isket, :isbra, :isoper, :issuper, :isoperket, :isoperbra)

export dnorm, dual_chan, tidyup, trans, transform, trunc_neg, unit
const methods_qobj =  (:dnorm, :dual_chan, :tidyup, :trans, :transform, :trunc_neg, :unit)

export eigenenergies, eigenstates, groundstate, matrix_element, tr
const methods = (:eigenenergies, :eigenstates, :groundstate, :matrix_element, :tr)

# renamed
# conj, expm, cosm, sinm, sqrtm, full, norm, permute, type 
