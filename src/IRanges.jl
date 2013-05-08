# IRanges Class for Julia BioSeq
#   Modeled after original BioConductor documentation
#   :http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
#
#   CMSC702, Spring 2013, Term Project (University of Maryland)
#   James Parker
#   Kristopher Micinski
#   Matthew Mauriello
#

module IRanges
using Option

export # TODO: make sure everything we need is exported.

# IRange data structure. Represents an integer interval with an optional name.
type IRange
	start::Int
	finish::Int
	width::Int
	name::option(String)

	# Constructor for IRange. It checks that finish >= start. The name argument is optional and defaults to None.
	IRange( start, finish, name=None) = begin
		if ( finish < start)
			error( "Finish must be greater than or equal to start.")
		end

		new( start, finish, finish - start + 1, name)
	end
end

# IRanges data structure.
type IRanges
	iranges::Vector{IRange}
	nameDict::option(Dict{String,Int}) # Maps names to an IRange index. Names must be unique and present in the iranges vector.

	IRanges( iranges::Vector{IRange}, index::Bool) = begin
		self = new( iranges, None)

		if index
			indexNames!( self)
		end

		self
	end
end

# Many operations like merge, disjoin, etc obliterate this index so only create index when finished manipulating iranges.
function indexNames!( iranges:IRanges)
	if !isNone( iranges.nameDict)
		error( "name index already exists")
	end

	dict = Dict{String,Int}()
	ranges = iranges.iranges
	sizehint( dict, length( ranges))

	for i = 1:lenth( ranges)
		range = ranges[i]
		name = range.name
		if !isNone( name)
			# Check that dict does not already contain name.
			if has( dict, name)
				error( "name, $name, is not unique")
			end

			# Add index to dictionary.
			merge!( dict, {name => i})
		end
	end
		
	iranges.nameDict = dict
end

end
