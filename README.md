# QuTiP.jl

<!-- 
[![Build Status](https://travis-ci.org/goropikari/QuTiP.jl.svg?branch=master)](https://travis-ci.org/goropikari/QuTiP.jl)

[![Coverage Status](https://coveralls.io/repos/goropikari/QuTiP.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/goropikari/QuTiP.jl?branch=master)

[![codecov.io](http://codecov.io/github/goropikari/QuTiP.jl/coverage.svg?branch=master)](http://codecov.io/github/goropikari/QuTiP.jl?branch=master)
-->

This package is a wrapper of [QuTiP](http://qutip.org/) using [PyCall](https://github.com/stevengj/PyCall.jl).

The [QuTiP](http://qutip.org/) package is a Python library for quantum information.

# Install

```julia
 Pkg.clone("https://github.com/goropikari/QuTiP.jl")
```

# Translate Python code to Julia code
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
    q[:dag]() # or dag(q)

    hinton(qidentity([2, 3])[:unit]()) # instead of identity use qidentity
```

# Renamed functions
original name --> renamed
- position --> qposition
- identity --> qidentity
