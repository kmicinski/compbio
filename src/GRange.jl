# GRanges Class for Julia BioSeq
#   Modeled after original BioConductor documentation
#   :http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
#
#   CMSC702, Spring 2013, Term Project (University of Maryland)
#   James Parker
#   Kristopher Micinski
#   Matthew Mauriello
#
#   See GRanges.jl for notes.

export GRange, start, finish, width, name, isless, isequal, isless, isequal

type GRange
    seqname::String
    range::IRange
    strand::Char
    score::Int #This is suppose to be metadata
    GC::Float64 #This is suppoe to be meatadata

    # Constructor
    GRange(seqname, range, strand, score, GC) = begin
        new (seqname, range, strand, score, GC)
    end

end

# [Utilized by countOverlap()]
# Getters for IRange data contained in GRange
function start(gr::GRange)
    return gr.range.start
end

function finish(gr::GRange)
    return gr.range.finish
end

function width(gr::GRange)
    return gr.range.width
end

function name(gr::GRange)
    return gr.range.name
end

# Define equality for IRange as equal when the ranges don't overlap at all.
function isless( r1::GRange, r2::GRange)
        finish(r1) < start(r2)
end

function isequal( r1::GRange, r2::GRange)
        !(r1 < r2 || r2 < r1)
end

# Define comparisons relative to a point.
function isless( x::Float64, r::GRange)
        x < start(r)
end

function isless( r::GRange, x::Float64)
        finish(r) < x
end

function isequal( x::Float64, r::GRange)
        !(x < r || r < x)
end

function isequal( r::GRange, x::Float64)
        isequal( x, r)
end

