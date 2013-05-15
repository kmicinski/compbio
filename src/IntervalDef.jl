# Abstract representation of intervals

module IntervalDef

export Interval, start, finish

abstract Interval

function start( i::Interval) # Not sure if these methods actually do anything..
	error( "`$(typeof(i))` does not implement `start`")
end
function finish( i::Interval) # Not sure if these methods actually do anything..
	error( "`$(typeof(i))` does not implement `finish`")
end

end
