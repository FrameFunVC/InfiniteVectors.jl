Plot(data::InfiniteVector, trailing...) =
    Plot(Options(), data, trailing...)

function Plot(options::Options, data::InfiniteVector{T}, trailing...) where T<:Number
    local s
    if haskey(options.dict, "samples_at")
        samples_at = options["samples_at"]
        if samples_at isa String
            s = parse_samples_at(samples_at)
        elseif samples_at isa UnitRange
            s = samples_at
        else
            throw(ArgumentError("samples_at option unknown."))
        end
    else
        s = -10:10
        options = @pgf {options..., samples_at =string(s)}
    end
    options = @pgf {options..., ycomb, mark="*"}
    if T <: Real
        Plot(options, Table([s, data[s]]), trailing...)
    else
        Plot(options, Table([vcat(s,s), vcat(real.(data[s]), imag.(data[s]))]), trailing...)
    end
end

function parse_samples_at(str)
    elements = split(str, ",")
    elements = map(x->parse(Int,x), elements[[1,2,4]])
    elements[1]:(elements[2]-elements[1]):elements[3]
end
