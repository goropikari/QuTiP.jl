# QuTiP.jl

[![Build Status](https://travis-ci.org/goropikari/QuTiP.jl.svg?branch=master)](https://travis-ci.org/goropikari/QuTiP.jl)

[![Coverage Status](https://coveralls.io/repos/goropikari/QuTiP.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/goropikari/QuTiP.jl?branch=master)

[![codecov.io](http://codecov.io/github/goropikari/QuTiP.jl/coverage.svg?branch=master)](http://codecov.io/github/goropikari/QuTiP.jl?branch=master)

*This package is under development. Do not use it for anything important*

This package is a wrapper of [QuTiP](http://qutip.org/) using [PyCall](https://github.com/stevengj/PyCall.jl).



# Install

```julia
Pkg.add("QuTiP")
# Pkg.checkout("QuTiP") # checkout master branch
```

# Translate Python code into Julia code
## Basic usage
Almost all syntax is same as original QuTiP, but some features are different.
As [PyCall](https://github.com/JuliaPy/PyCall.jl)'s troubleshooting says, use `foo[:bar]` and `foo[:bar](...)` rather than `foo.bar` and `foo.bar(...)`, respectively.

```python
# python
# quoted from [QuTiP tutorial](http://nbviewer.jupyter.org/github/qutip/qutip-notebooks/blob/master/examples/superop-contract.ipynb)
import numpy as np
import qutip as qt
%matplotlib inline
qt.settings.colorblind_safe = True

q = basis(2,0)
q.dag()

qt.visualization.hinton(qt.identity([2, 3]).unit());
```

```julia
# julia
using QuTiP, PyPlot

q = basis(2,0)
dag(q) # or q'

hinton(qidentity([[2], [3]])[:unit]()) # instead of identity use qidentity
```

## If arguments of QuTiP function contain Julia functions
Julia functions cannot be passed to python as native python functions.
So we have to use `py"..."`.

For examples,
```python
# python
# from http://qutip.org/docs/latest/guide/dynamics/dynamics-time.html
def H1_coeff(t, args):
    return 9 * np.exp(-(t / 5.) ** 2)

H = [H0,[H1,H1_coeff]]
output = mesolve(H, psi0, t, c_ops, [ada, sigma_UU, sigma_GG])
```



```julia
# Julia
using PyCall

fnj2py = py"lambda f: lambda t, args: f(t, args)"
H1_coeff(t, args) = 9 * exp(-(t / 5.) ^ 2)

H = [H0,[H1, fnj2py(H1_coeff)]]
output = mesolve(H, psi0, t, c_ops, [ada, sigma_UU, sigma_GG])
```

# Examples
In order to test this package and compare python and julia, I translate some Jupyter notebooks about qutip into Julia.
All original python codes are left as comment.  

From [jrjohansson/qutip-lectures](https://github.com/jrjohansson/qutip-lectures)
- [Lecture 0 Introduction to QuTiP - The Quantum Toolbox in Python](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-0-Introduction-to-QuTiP.ipynb)
- [Lecture 1 QuTiP lecture: Vacuum Rabi oscillations in the Jaynes-Cummings model](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-1-Jaynes-Cumming-model.ipynb)
- [Lecture 2A Cavity Qubit Gates](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-2A-Cavity-Qubit-Gates.ipynb)
- [Lecture 2B QuTiP lecture: Single-Atom-Lasing](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-2B-Single-Atom-Lasing.ipynb)
- [Lecture 3A Dicke model](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-3A-Dicke-model.ipynb)
- [Lecture 3B Jaynes Cumming model with ultrastrong coupling](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-3B-Jaynes-Cumming-model-with-ultrastrong-coupling.ipynb)
- [Lecture 4 QuTiP lecture: Correlation functions](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-4-Correlation-Functions.ipynb)
- [Lecture 5 QuTiP lecture: Evolution and quantum statistics of a quantum parameter amplifier](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-5-Parametric-Amplifier.ipynb)
- [Lecture 6 QuTiP lecture: Quantum Monte-Carlo Trajectories](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-6-Quantum-Monte-Carlo-Trajectories.ipynb)
- [Lecture 7 Two-qubit iSWAP gate and process tomography](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-7-iSWAP-gate.ipynb)
- [Lecture 8 Adiabatic quantum computing](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-8-Adiabatic-quantum-computing.ipynb)
- [Lecture 13 Resonance flourescence](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-13-Resonance-flourescence.ipynb)
- [Lecture 15: Nonclassically driven atoms (cascaded quantum systems)](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-15-Nonclassically-driven-atoms.ipynb)
- [Lecture 16: Gallery of Wigner functions](http://nbviewer.jupyter.org/github/goropikari/qutip-lectures/blob/For_qutip_v4.2.0/julia/Lecture-16-Gallery-of-Wigner-functions.ipynb)

# convert Qobj to Julia array
To convert Oobj to julia array, use `full`.
```julia
julia> x = basis(2,0)
QuTiP.Quantum(PyObject Quantum object: dims = [[2], [1]], shape = (2, 1), type = ket
Qobj data =
[[ 1.]
 [ 0.]]

julia> full(x)
2×1 Array{Complex{Float64},2}:
 1.0+0.0im
 0.0+0.0im

julia> full(sigmax())
2×2 Array{Complex{Float64},2}:
 0.0+0.0im  1.0+0.0im
 1.0+0.0im  0.0+0.0im
```

# Renamed functions
In order to avoid name conflict, some functions are renamed.  
original name --> renamed
- position --> qposition
- identity --> qidentity
- num      --> qnum
- squeeze  --> qsqueeze

# Contributing
I appreciate all kinds of help. Please feel free to open an issue.
