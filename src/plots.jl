

for plot in (:Plot, :PlotInc)
    @eval begin
        $(plot)(data::InfiniteVector, trailing...) =
            $(plot)(Options(), data, trailing...)

        function $(plot)(options::Options, data::InfiniteVector{T}, trailing...) where T<:Number
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
                s = default_range(data)
                options = @pgf {options..., samples_at =string(s)}
            end
            options = @pgf {options..., ycomb, mark="*"}
            if T <: Real
                $(plot)(options, Table([s, data[s]]), trailing...)
            else
                $(plot)(options, Table([vcat(s,s), vcat(real.(data[s]), imag.(data[s]))]), trailing...)
            end
        end
    end
end

default_range(::InfiniteVector) = -10:10
default_range(vec::PeriodicInfiniteVector) = -period(vec)>>1:period(vec)>>1+1

function parse_samples_at(str)
    elements = split(str, ",")
    elements = map(x->parse(Int,x), elements[[1,2,4]])
    elements[1]:(elements[2]-elements[1]):elements[3]
end
