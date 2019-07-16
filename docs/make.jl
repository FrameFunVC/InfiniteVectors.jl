using Pkg
Pkg.develop(PackageSpec(path=splitdir(@__DIR__)[1]))
pkg"instantiate"
using Documenter, InfiniteVectors

const render_pdf = "pdf" in ARGS
let r = r"buildroot=(.+)", i = findfirst(x -> occursin(r, x), ARGS)
    global const buildroot = i === nothing ? (@__DIR__) : first(match(r, ARGS[i]).captures)
end

const format = if render_pdf
    LaTeX(
        platform = "texplatform=docker" in ARGS ? "docker" : "native"
    )
else
    Documenter.HTML(
        prettyurls = ("deploy" in ARGS),
    )
end

makedocs(sitename="InfiniteVectors.jl",
    modules = [InfiniteVectors],
    authors = "vincentcp",
    format = format,
    pages = [
        "Home" => "index.md",
        "Manual" =>Any["man/indices.md","man/vectors.md","man/vectoroperations.md"],
        "Developer information" =>["development.md"]
        ],
    doctest=true
)

if "deploy" in ARGS && Sys.ARCH === :x86_64 && Sys.KERNEL === :Linux
    deploydocs(
        repo = "github.com/FrameFunVC/InfiniteVectors.jl.git",
    )
end
