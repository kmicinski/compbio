# TODO ...
# Inspired by: https://github.com/JuliaLang/julia/blob/d7efe83cbbd9312efa14620d857c30f8c90ec289/examples/tree.jl

module ITree

export Interval, IntervalTree, intervalTree, findOverlaps, countOverlaps #, EmptyNode, IntervalNode

abstract Interval
#function start( i::Interval) # Not sure if these methods actually do anything..
#	error( "`$(typeof(i))` does not implement `start`")
#end
#function finish( i::Interval) # Not sure if these methods actually do anything..
#	error( "`$(typeof(i))` does not implement `finish`")
#end
# comparisons...?

abstract IntervalTree

type EmptyNode <: IntervalTree # TODO: can i make this singleton?
end

type IntervalNode <: IntervalTree
	center::Float64
	forwardIntervals::Vector{Interval}
	backwardIntervals::Vector{Interval}
	left::IntervalTree
	right::IntervalTree

	# TODO: Can this be made tail recursive?
	function IntervalNode{T <: Interval}( intervals::Vector{T})
		# Check if the interval is empty.
		if length( intervals) == 0
			EmptyNode() # TODO: HERE
		end
		center = middle( intervals)

		# TODO: can I prealloc a length? Or use a list or something..
		leftIntervals = Array( Interval, 0)
		overlapIntervals = Array( Interval, 0)
		rightIntervals = Array( Interval, 0)

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
function intervalTree{T <: Interval}( intervals::Vector{T})
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
function middle{T <: Interval}( intervals::Vector{T})
#function middle( intervals::Vector{Interval})
	minn = Inf
	maxx = -Inf

	for i in intervals
		println( typeof(i))
		minn = min( minn, start(i))
		maxx = max( maxx, finish(i))
	end

	(minn + maxx) / 2.0
end

# Checks if the tree is the empty node.
function isEmpty( tree::IntervalTree)
	isa( tree, EmptyNode{T}) # TODO: HERE
end

# Find and return a vector of all intervals in tree that overlap query.
# I didn't exactly follow the textbook on this one but I think it should work...
# TODO: TEST!!!
function findOverlaps{T <: Interval}( root::IntervalTree, query::T)
	overlaps = Array( T, 0)
	queue = [root] # TODO: Change this to Queue?

	while length( queue) > 0
		# Pop the next tree node.
		tree = pop!( queue)

		# Check if the node is empty.
		if isEmpty( Tree) # Maybe dispatch would be cleaner here...
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
					append!( overlaps, interval)
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
					append!( overlaps, interval)
				else
					break
				end
			end

			push!( queue, tree.right)
		else
			error( "Interval Tree implementation error. This should never get called.")
		end
	end

	overlaps
end

# Count the intervals in tree that overlap query.
# TODO: Wastes memory but works for now...
function countOverlaps{T <: Interval}( tree::IntervalTree, query::T)
	length( findOverlaps( tree, query))
end

end
