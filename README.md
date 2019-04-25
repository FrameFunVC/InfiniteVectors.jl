# InfiniteVectors.jl

*A Julia package for infinite vectors operations such as convolution and z-transform.*



| **Documentation** | **Build Status** | **Coverage** |
|-------------------|------------------|--------------|
| [![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://vincentcp.github.io/InfiniteVectors.jl/stable)  [![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://vincentcp.github.io/InfiniteVectors.jl/dev) | [![Build Status](https://travis-ci.org/vincentcp/InfiniteVectors.jl.png)](https://travis-ci.org/vincentcp/InfiniteVectors.jl) [![Build status](https://ci.appveyor.com/api/projects/status/gh4ka7m9a7qekqu8?svg=true)](https://ci.appveyor.com/project/vincentcp/InfiniteVectors-jl) | [![Coverage](https://codecov.io/gh/vincentcp/InfiniteVectors.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/vincentcp/InfiniteVectors.jl)  [![Coverage Status](https://coveralls.io/repos/github/vincentcp/InfiniteVectors.jl/badge.svg)](https://coveralls.io/github/vincentcp/InfiniteVectors.jl) |


Installation instructions and full functionality can be found in the [documentation](https://vincentcp.github.io/InfiniteVectors.jl/dev)

## Examples 

```julia
julia> using InfiniteVectors, PGFPlotsX
julia> Plot(CompactInfiniteVector(1:4,-1))
```
![svg](test/example1.svg)
```julia
julia> Plot(PeriodicInfiniteVector(1:10))
```
![svg](test/example2.svg)
