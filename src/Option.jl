# Module to emulate options in OCaml (although not that rigorously).

module Option
export isNone, option

isNone(x) = x == nothing
option(T) = Union(T,Nothing)

end
