# GRanges Class for Julia BioSeq
#   Modeled after original BioConductor documentation
#   :http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
#
#   CMSC702, Spring 2013, Term Project (University of Maryland)
#   James Parker
#   Kristopher Micinski
#   Matthew Mauriello
#

export GRange, GRanges


# GRanges type; poorly approximated
# TO DO: Write simple constructor? Seperate into metadata and GRange data using matrix structure.
# TO DO: Change 'type' to  'Immutable'?
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

type GRanges
    granges::Vector{GRange}

    #Constructor
    GRanges( granges::Vector{GRange}) = begin
        new(granges)
    end

end

# Show the structure and data contained within the GRanges object; decently approximated
# This is all assigned directly at initialization; however, you could imagine a constructor that
# automated some of these fields, maybe?
function show(gr::GRanges)
    println(string("\tGRanges with ", length(gr.seqnames)," ranges and ", 2, " metadata columns:"))
    println("\t\tseqnames \tranges \t\tstrand \t | \tscore \tGC")
    println("\t\t<ASCIIString> \t<IRanges> \t<Char> \t | \t<Int> \t<Float64>")
    for i = 1:length(gr.seqnames)
        println(string("\t", gr.names[i],"\t", gr.seqnames[i], "\t\t", "[",gr.ranges[i].start, " ", gr.ranges[i].finish,"]", "\t\t", gr.strand[i], "\t | \t", gr.score[i], "\t", gr.GC[i]))
    end
    println("\t---")
    print("\tseqlengths:")
    print("\n\t")
    for i = 1:length(gr.seqnames)
        print(string(gr.seqnames[i], "\t"))
    end
    print("\n\t")
    for i = 1:length(gr.seqnames)
        # TO DO: this should be the total length of the sequence?
        print(string("NA", "\t"))
    end
    print("\n")
end

# seqnames appears to be just another view of the data
# TO DO: figure out what runs are. Sounds silly, I know.
function seqnames(gr::GRanges)
    println(string("\t", typeof(gr.seqnames), " of length ", length(gr.seqnames), " with X runs"))

    print("\tLengths:\t")
    for i = 1:length(gr.seqnames)
        print(string(gr.ranges[i].width,"\t"))
    end
    print("\n")

    print("\tValues:\t\t")
    for i = 1:length(gr.seqnames)
        print(string(gr.seqnames[i],"\t"))
    end
    print("\n")
    # TO DO: THIS PART IS MOSTLY NON-FUNCTIONAL
    # Count of unique sequences?
    print("\tLevels(X):\t")
    for i = 1:length(gr.seqnames)
        # List of unique sequences?
        print(string("X","\t"))
    end
    print("\n")
end

function ranges(gr::GRanges)
    print("\tIRanges of length: ");
    println(length(gr.seqnames));
    println("\t\tstart\tend\twidth\tnames");
    for i = 1:length(gr.seqnames)
        println(string("\t[", i, "]\t",gr.ranges[i].start, "\t", gr.ranges[i].finish, "\t", ((gr.ranges[i].finish - gr.ranges[i].start) + 1), "\t", gr.names[i]));
    end
end

function strand(gr::GRanges)
end

# Right now I have poorly inserted what is termed metadata at the top; when writing this function we'll have to define some kind of structure and matrix
# for holding metadata. (Impact: GRanges, show, and probably most other functions if this isn't handled soon
function mcols(gr::GRanges)
end


function seqlengths(gr::GRanges)
end

function names(gr::GRanges)
end

# length renamed to len as to not conflict with Julia default length function? Perhaps, just use the Length()
function len(gr::GRanges)
end

function splitGRanges(gr::GRanges, each::Int)
end

function c(gr::GRanges, br::GRanges)
end

#TO DO: Define Subsetting method

#TO DO: Define seqselect method

function head(gr::GRanges, n::Int)
end

function rep(gr::GRanges, times::Int)
end

function rev(gr::GRanges)
end

function tail(gr::GRanges, n::Int)
end

function window(gr::GRanges, start::Int, finish::Int)
end

function start(gr::GRanges)
end

function finish(gr::GRanges)
end

function width(gr::GRanges)
end

function range(gr::GRanges)
end

function flank(gr::GRanges, n::Int, start::Bool)
end

function shift(gr::GRanges, n::Int)
end

function resize(gr::GRanges, n::Int)
end

function reduce(gr::GRanges)
end

function gaps(gr::GRanges)
end

function disjoin(gr::GRanges)
end

function converage(gr::GRanges)
end

function union(gr::GRanges, br::GRanges)
end

function intersect(gr::GRanges, br::GRanges)
end

function setdiff(gr::GRanges, br::GRanges)
end

# TO DO: Investigate "p" methods.

# TO DO: Investigate GRangesList.

# TO DO: Investigate Gapped Alignments.

# Example usage:
#gr = GRanges(["a", "b", "c"],["chr1", "chr2", "chr3"], [IRanges(1,7,7-1,0), IRanges(0,0,0,0), IRanges(0,0,0,0)] , ['+','-','*'] , [1,2,3],  [0.0,0.0,0.0])

