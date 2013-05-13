# TODO ...
# Inspired by: https://github.com/JuliaLang/julia/blob/d7efe83cbbd9312efa14620d857c30f8c90ec289/examples/tree.jl

module ITree

export Interval, IntervalTree, intervalTree #, EmptyNode, IntervalNode

abstract Interval
abstract IntervalTree{T <: Interval}

type EmptyNode{T} <: IntervalTree{T} # TODO: can i make this singleton?
end

type IntervalNode{T} <: IntervalTree{T}
	center::Float64
	forwardIntervals::Vector{T}
	backwardIntervals::Vector{T}
	left::IntervalTree{T}
	right::IntervalTree{T}

	# TODO: Can this be made tail recursive?
	function IntervalNode( intervals::Vector{T})
		center = middle( intervals)

		# TODO: can I prealloc a length? Or use a list or something..
		leftIntervals = Array( Interval, 0)
		overlapIntervals = Array( Interval, 0)
		rightIntervals = Array( Interval, 0)

		# Add intervals to their appropriate sets.
		for i in intervals
			if i < center
				push!( leftIntervals, i)
			else if i > center
				push!( rightINtervals, i)
			else
				push!( overlapIntervals, i)
			end
		end

		forwardIntervals = sortby( overlapIntervals, i->i.start)
		sortby!( overlapIntervals, i->i.finish)
		backwardIntervals = overlapIntervals

		leftNode = if length( leftIntervals) == 0
			EmptyNode()
		else
			IntervalNode( leftIntervals)
		end

		rightNode = if length( rightIntervals) == 0
			EmptyNode()
		else
			IntervalNode( rightIntervals)
		end

		new( center, forwardIntervals, backwardIntervals, leftNode, rightNode)
	end
end

function intervalTree{T <: Interval}( intervals::Vector{T})
	if length( intervals) == 0
		EmptyNode()
	else
		IntervalNode( intervals)
	end
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

# Finds the middle of the intervals.
# Somewhat naive. Could improve to split more evenly.
function middle{T <: Interval}( intervals::Vector{T})
	minn = Inf
	maxx = -Inf

	for i in intervals
		minn = min( minn, i.start)
		maxx = max( maxx, i.finish)
	end

	(minn + maxx) / 2.0
end

end
