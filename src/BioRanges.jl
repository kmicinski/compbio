# TODO:doc

module BioRanges
        using IntervalDef
	using ITree
	using Option
        include("IntervalDef.jl")
	include("IRanges.jl")
        include("GRange.jl")
	include("GRanges.jl")
        include("GRangesList.jl")
end
