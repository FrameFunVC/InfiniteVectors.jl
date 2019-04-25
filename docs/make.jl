using Pkg
pkg"rm Sequences"
Pkg.develop(PackageSpec(path=splitdir(@__DIR__)[1]))
pkg"instantiate"
using Documenter, Sequences

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

makedocs(sitename="Sequences.jl",
    modules = [Sequences],
    authors = "vincentcp",
    format = format,
    pages = [
        "Home" => "index.md",
        ],
    doctest=true
)

if "deploy" in ARGS && Sys.ARCH === :x86_64 && Sys.KERNEL === :Linux
    deploydocs(
        repo = "github.com/vincentcp/Sequences.jl.git",
    )
end
