using QuTiP, PyCall
using Base.Test
@pyimport qutip as qt

include("display.jl")
include("bloch_redfield.jl")
include("entropy.jl")
include("eseries.jl")
include("metrics.jl")
include("qobj.jl")
include("utilities.jl")
include("wigner.jl")

@test try
   QuTiP.about()
   true
catch
   false
end
