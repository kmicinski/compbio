# TODO ...
# Inspired by: https://github.com/JuliaLang/julia/blob/d7efe83cbbd9312efa14620d857c30f8c90ec289/examples/tree.jl

module ITree

export IntervalTree, EmptyNode, IntervalNode, intervalTree

abstract IntervalTree{T <: Interval}

type EmptyNode{T} <: IntervalTree{T} # TODO: can i make this singleton?
end

type IntervalNode{T} <: IntervalTree{T}
	center::Float64
	forwardIntervals::Vector{T}
	backwardIntervals::Vector{T}
	left::IntervalTree{T}
	right::IntervalTree{T}
end

function intervalTree{T <: Interval}( intervals::Vector{T})
	tree = EmptyNode()
	for i in intervals
		insert!( tree, i)
	end
	tree
end

function insert!{T <: Interval}( tree::IntervalTree{T}, interval::T)
	# TODO ...

	#empty = EmptyNode()
	#node = IntervalNode( interval, empty, empty)
	#insert!( tree, node)
end

function insert!{T <: Interval}( tree::EmptyNode{T}, node::IntervalNode{T})
	node
end

function insert!{T <: Interval}( tree::IntervalNode{T}, node::IntervalNode{T})
	# TODO ...
end

end
