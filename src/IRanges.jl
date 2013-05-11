# IRanges Class for Julia BioSeq
#   Modeled after original BioConductor documentation
#   :http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
#
#   CMSC702, Spring 2013, Term Project (University of Maryland)
#   James Parker
#   Kristopher Micinski
#   Matthew Mauriello
#

using Option
export IRange, IRanges, indexNames!, start, finish, width # TODO: make sure everything we need is exported.

# IRange data structure. Represents an integer interval with an optional name.
type IRange
	start::Int
	finish::Int
	width::Int
	name::option(String)

	# Constructor for IRange. It checks that finish >= start. The name argument is optional and defaults to None.
	IRange( start, finish, name) = begin
		if ( finish < start)
			error( "Finish must be greater than or equal to start.")
		end

		new( start, finish, finish - start + 1, name)
	end
end

# Define equality for IRange as equal when both endpoints are equal.
#function isequal( r1::IRange, r2::IRange)
#	(r1.start == r2.start) && (r1.finish == r2.finish)
#end

# A stupid comparison that just checks if r2.start is less than r2.start.
#function isless( r1::IRange, r2::IRange)
#	r1.start < r2.start
#end

# IRanges data structure.
type IRanges
        iranges::Vector{IRange} # Should this be changed to a btree???
	nameDict::option(Dict{String,Int}) # Maps names to an IRange index. Names must be unique and present in the iranges vector.

        IRanges( iranges::Vector{IRange}, index::Bool) = begin
                self = new(iranges, Nothing)

                if index
                       indexNames!(self)
                end

		self
        end
end

# Many operations like merge, disjoin, etc obliterate this index so only create index when finished manipulating iranges.
function indexNames!( iranges::IRanges)

        #if !isNone( iranges.nameDict)
        #	error( "name index already exists")
        #end

	dict = Dict{String,Int}()
	ranges = iranges.iranges
	sizehint( dict, length( ranges))

        for i = 1:length( ranges)
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

function start( ranges::IRanges)
	[ range.start for range in ranges.iranges ]
end

function finish( ranges::IRanges)
	[ range.finish for range in ranges.iranges ]
end

function width( ranges::IRanges)
	[ range.width for range in ranges.iranges ]
end

