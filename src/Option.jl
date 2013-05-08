# Module to emulate options in OCaml (although not that rigorously).

module Option

isNone(x) = x == None
option(T) = Union(T,UnionKind)

end
