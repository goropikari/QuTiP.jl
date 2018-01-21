@test try
    show(STDOUT, "text/latex", basis(2,0))
    true
catch
    false 
end

# @test try
#     show(STDOUT, "text/latex", [basis(2,0), basis(2,1)])
#     true
# catch
#     false 
# end


omega = 1.0
es = eseries(0.5 * sigmax(), im * omega) + eseries(0.5 * sigmax(), -im * omega)
@test try
    show(STDOUT, "text/latex", es)
    true
catch
    false
end

