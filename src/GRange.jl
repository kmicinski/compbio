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

export GRange, start, finish, width, range

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
