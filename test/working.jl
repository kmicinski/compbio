# IRanges Class for Julia BioSeq
#   Modeled after original BioConductor documentation
#   :http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
#
#   CMSC702, Spring 2013, Term Project (University of Maryland)
#   James Parker
#   Kristopher Micinski
#   Matthew Mauriello
#

#export IRange, IRanges, start, finish, width, isless, isequal # TODO: make sure everything we need is exported.

using Option

abstract Interval

function start( i::Interval) # Not sure if these methods actually do anything..
	error( "`$(typeof(i))` does not implement `start`")
end
function finish( i::Interval) # Not sure if these methods actually do anything..
	error( "`$(typeof(i))` does not implement `finish`")
end


# IRange data structure. Represents an integer interval with an optional name.
type IRange <: Interval # TODO: change to immutable once everyone has 2.0
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

# Getters for IRange.
function start(r :: IRange)
	r.start
end

function finish(r::IRange)
	r.finish
end

import Base.isless
import Base.isequal

# Define equality for IRange as equal when the ranges don't overlap at all.
function isless( r1::IRange, r2::IRange)
	r1.finish < r2.start
end

function isequal( r1::IRange, r2::IRange)
	!(r1 < r2 || r2 < r1)
end

# Define comparisons relative to a point.
function isless( x::Float64, r::IRange)
	x < r.start
end

function isless( r::IRange, x::Float64)
	r.finish < x
end

function isequal( x::Float64, r::IRange)
	!(x < r || r < x)
end

function isequal( r::IRange, x::Float64)
	isequal( x, r)
end

# An alternative definition for comparisons.
## Define equality for IRange as equal when both endpoints are equal.
#function isequal( r1::IRange, r2::IRange)
#	(r1.start == r2.start) && (r1.finish == r2.finish)
#end

## A stupid comparison that just checks if r2.start is less than r2.start.
#function isless( r1::IRange, r2::IRange)
#	r1.start < r2.start
#end

# IRanges data structure.
type IRanges # TODO: change to immutable once everyone has 2.0
	iranges::Vector{IRange} # Should this be changed to a btree???
	nameDict::option(Dict{String,Int}) # Maps names to an IRange index. Names must be unique and present in the iranges vector.

	IRanges( iranges::Vector{IRange}, index::Bool) = begin
		if index
			dict = indexNames( iranges)
			new( iranges, dict)
		else
			new( iranges, nothing)
		end
	end
end

# Create a dictionary to map names to ranges.
function indexNames!( ranges::Vector{IRange})

	#if !isNone( iranges.nameDict)
	#	error( "name index already exists")
	#end

	dict = Dict{String,Int}()
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

# TODO ...
# Inspired by: https://github.com/JuliaLang/julia/blob/d7efe83cbbd9312efa14620d857c30f8c90ec289/examples/tree.jl

#module ITree

#using IntervalDef

#include("IRanges.jl")

#export IntervalTree, intervalTree, findOverlaps, countOverlaps #, EmptyNode, IntervalNode

abstract IntervalTree

type EmptyNode <: IntervalTree # TODO: can i make this singleton?
end

type IntervalNode <: IntervalTree
	center::Float64
	forwardIntervals::Vector{IRange}
	backwardIntervals::Vector{IRange}
	left::IntervalTree
	right::IntervalTree

	# TODO: Can this be made tail recursive?
	function IntervalNode( intervals::Vector)
		# Check if the interval is empty.
		if length( intervals) == 0
			return EmptyNode() # TODO: HERE
		end
		center = middle( intervals)

		# TODO: can I prealloc a length? Or use a list or something..
		leftIntervals = Array( IRange, 0)
		overlapIntervals = Array( IRange, 0)
		rightIntervals = Array( IRange, 0)

		# Add intervals to their appropriate sets.
		for i in intervals
			if i < center
				push!( leftIntervals, i)
			elseif i > center
				push!( rightIntervals, i)
			else
				push!( overlapIntervals, i)
			end
		end

		forwardIntervals = sortby( overlapIntervals, i->start(i))
		sortby!( overlapIntervals, i->finish(i))
		backwardIntervals = overlapIntervals

		leftNode = IntervalNode( leftIntervals)
		rightNode = IntervalNode( rightIntervals)

		new( center, forwardIntervals, backwardIntervals, leftNode, rightNode)
	end
end

# Public constructor.
function intervalTree( intervals::Vector)
	IntervalNode( intervals) # TODO: HERE
end

#function insert!{T <: Interval}( tree::IntervalTree{T}, interval::T)
#	# TODO ...
#
#	#empty = EmptyNode()
#	#node = IntervalNode( interval, empty, empty)
#	#insert!( tree, node)
#end
#
#function insert!{T <: Interval}( tree::EmptyNode{T}, node::IntervalNode{T})
#	node
#end
#
#function insert!{T <: Interval}( tree::IntervalNode{T}, node::IntervalNode{T})
#	# TODO ...
#end

# Finds the middle of the intervals. Assumes that intervals is not empty.
# Somewhat naive. Could improve to split more evenly.
function middle( intervals::Vector)
#function middle( intervals::Vector{Interval})
	minn = Inf
	maxx = -Inf

	for i in intervals
		minn = min( minn, start(i))
		maxx = max( maxx, finish(i))
	end

	(minn + maxx) / 2.0
end

# Checks if the tree is the empty node.
function isEmpty( tree::IntervalTree)
	isa( tree, EmptyNode) # TODO: HERE
end

# Find and return a vector of all intervals in tree that overlap query.
# I didn't exactly follow the textbook on this one but I think it should work...
# TODO: TEST!!!
function findOverlaps( root, query)
	overlaps = Array( IRange, 0)
	queue = Array(IntervalTree, 1)
	queue[1] = root # TODO: Change this to Queue?

	while length( queue) > 0
		# Pop the next tree node.
		tree = pop!( queue)

		# Check if the node is empty.
		if isEmpty( tree) # Maybe dispatch would be cleaner here...
			break
		end

		# Check if center overlaps with query.
		if query == tree.center
			append!( overlaps, tree.forwardIntervals)

			push!( queue, tree.left)
			push!( queue, tree.right)
		# Check if query is left of center.
		elseif query < tree.center
			for i in 1:length(tree.forwardIntervals)
				interval = tree.forwardIntervals[i]
				if start(interval) <= finish(query)
					push!( overlaps, interval)
				else
					break
				end
			end

			push!( queue, tree.left)
		# Check if query is right of center.
		elseif tree.center < query
			for i in reverse(1:length(tree.backwardIntervals))
				interval = tree.backwardIntervals[i]
				if finish(interval) >= start(query)
					push!( overlaps, interval)
				else
					break
				end
			end

			push!( queue, tree.right)
		else
			println( tree.center)
			println( tree.center == query)
			println( query)
			error( "Interval Tree implementation error. This should never get called.")
		end
	end

	overlaps
end

# Count the intervals in tree that overlap query.
# TODO: Wastes memory but works for now...
function countOverlaps( tree::IntervalTree, query::IRange)
	length( findOverlaps( tree, query))
end

#end
# Abstract representation of intervals

#module IntervalDef

#export Interval, start, finish

#end
# TODO: ...

#require("../src/ITree.jl") # Do you have to do this?
#require("../src/Option.jl") # Do you have to do this?
#require("../src/BioRanges.jl") # Do you have to do this?

using Base.Test
#using BioRanges
#using ITree

#import ITree.*

# Set seed.
srand(63)

# Construct random intervals.
n = 1000000
starts = rand( 1:1000, n)
finishes = starts + rand( 190:210, n)
intervals = Array( IRange, n)
for i in 1:n
	intervals[i] = IRange( starts[i], finishes[i], nothing)
end

# Test middle.
#middle(intervals) == (min( starts) + max(finishes)) / 2.0
#@test middle(intervals) == (min( starts) + max(finishes)) / 2.0

# Test countOverlaps.

tree = intervalTree( intervals)

function construct()
	tree = intervalTree( intervals)
end

function query1()
	q = IRange( 0, 1000, nothing)
	#findOverlaps( tree, i)#IRange( 0, 1000, nothing))
	countOverlaps( tree, q)
end

function query2()
	q = IRange( 100, 300, nothing)
	#findOverlaps( tree, i)#IRange( 0, 1000, nothing))
	countOverlaps( tree, q)
end

function query3()
	q = IRange( 650, 850, nothing)
	#findOverlaps( tree, i)#IRange( 0, 1000, nothing))
	countOverlaps( tree, q)
end

#using Benchmark
function timer(f::Function, n::Integer)
    # Call once to force JIT compilation
    f()

    times = Array(Float64, n)
    for itr in 1:n
        times[itr] = @elapsed f()
    end

    return times
end

println(timer(construct,5))
println(timer(query1,5))
println(timer(query2,5))
println(timer(query3,5))
#benchmark( construct, "construction", 1)
#benchmark( query1, "query1", 1)
#benchmark( query2, "query2", 1)
#benchmark( query3, "query3", 1)
