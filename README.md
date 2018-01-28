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
