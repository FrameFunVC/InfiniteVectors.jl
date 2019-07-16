# InfiniteVectors.jl Documentation

*A Julia package for infinite vectors operations such as convolution, z-transform.*

For installation instructions, see [Installation](@ref).

For a  full description of the functionality use the manual:
```@contents
Pages = ["man/indices.md",
            "man/doublyinfinite.md"]
```

## Installation

InfiniteVectors.jl is not added to `General.jl`.
The package can easily be installed by cloning its git repository. From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add https://github.com/FrameFunVC/InfiniteVectors.jl
```

or in a file you could use

```julia
using Pkg
pkg"add https://github.com/FrameFunVC/InfiniteVectors.jl"
```


# Development

For development instructions see the [Development](@ref development)
