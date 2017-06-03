# QuTiP.jl

[![Build Status](https://travis-ci.org/goropikari/QuTiP.jl.svg?branch=master)](https://travis-ci.org/goropikari/QuTiP.jl)

[![Coverage Status](https://coveralls.io/repos/goropikari/QuTiP.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/goropikari/QuTiP.jl?branch=master)

[![codecov.io](http://codecov.io/github/goropikari/QuTiP.jl/coverage.svg?branch=master)](http://codecov.io/github/goropikari/QuTiP.jl?branch=master)


This package is a wrapper of [QuTiP](http://qutip.org/) using [PyCall](https://github.com/stevengj/PyCall.jl).


# Install

```julia
Pkg.add("QuTiP")
```

# Translate Python code into Julia code
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
dag(q)

hinton(qidentity([[2], [3]])[:unit]()) # instead of identity use qidentity
```

To test this package and compare python and julia, I translate some Jupyter notebooks about qutip into Julia. 
All original python codes are left as comment.  

From [jrjohansson/qutip-lectures](https://github.com/jrjohansson/qutip-lectures)
- [Lecture 0 - Introduction to QuTiP - The Quantum Toolbox in Python](https://github.com/goropikari/qutip-lectures/blob/master/Lecture-0-Introduction-to-QuTiP.ipynb)

From [qutip/qutip-notebooks](https://github.com/qutip/qutip-notebooks)
- [QuTiP example: Correlation functions and spectrum of a atom-cavity system](https://github.com/goropikari/qutip-notebooks/blob/master/examples/atom-cavity-correlation-function.ipynb)
- [QuTiP example: Dynamics of an atom-cavity system using three different solvers](https://github.com/goropikari/qutip-notebooks/blob/master/examples/atom-cavity-dynamics.ipynb)
- [QuTiP example: Bloch-Redfield Master Equation (bloch-redfield.ipynb)](https://github.com/goropikari/qutip-notebooks/blob/master/examples/bloch-redfield.ipynb)
- [QuTiP example: Bloch-Redfield Master Equation (brmesolve.ipynb)](https://github.com/goropikari/qutip-notebooks/blob/master/examples/brmesolve.ipynb)
- [QuTiP example: Single-Qubit Dynamics](https://github.com/goropikari/qutip-notebooks/blob/master/examples/qubit-dynamics.ipynb)
- [QuTiP example: Energy-levels of a quantum systems as a function of a single parameter](https://github.com/goropikari/qutip-notebooks/blob/master/examples/energy-levels.ipynb)
- [QuTiP example: eseries](https://github.com/goropikari/qutip-notebooks/blob/master/examples/eseries.ipynb)
- [QuTiP example: Dynamics of a Spin Chain](https://github.com/goropikari/qutip-notebooks/blob/master/examples/spin-chain.ipynb)

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
In order to avoid name conflict, some functions are rennamed.  
original name --> renamed
- position --> qposition
- identity --> qidentity
- num      --> qnum
- squeeze  --> qsqueeze
