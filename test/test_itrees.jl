# TODO: ...

require("../src/ITree.jl") # Do you have to do this?
require("../src/Option.jl") # Do you have to do this?
require("../src/BioRanges.jl") # Do you have to do this?

using Base.Test
using BioRanges
using ITree

#import ITree.*

# Set seed.
srand(63)

# Construct random intervals.
n = 100
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
i = IRange( 0, 1000, nothing)
println("here")
findOverlaps( tree, i)#IRange( 0, 1000, nothing))
#@test countOverlaps( tree, IRange( 0, 1000, nothing)) == n

# TODO: test boundaries, empty query, empty tree, etc...
