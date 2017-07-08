omega = 1.0
es1 = eseries(0.5 * sigmax(), im * omega) + eseries(0.5 * sigmax(), -im * omega)
es1qt = qt.eseries(qt.sigmax() * 0.5, im * omega) + qt.eseries(qt.sigmax() * 0.5, -im * omega)
@test ampl(es1) == es1qt[:ampl]
@test rates(es1) == es1qt[:rates]
@test dims(es1) == es1qt[:dims]
@test shape(es1) == es1qt[:shape]
@test value(es1, 0:0.5:1) == es1qt[:value](0:0.5:1)
@test tidyup(es1)[:ampl] == es1qt[:tidyup]()[:ampl]
