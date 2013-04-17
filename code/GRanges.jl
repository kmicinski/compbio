# GRanges/IRanges Class for Julia BioSeq
#   Modeled after original BioConductor documentation
#   :http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
#
#   CMSC702, Spring 2013, Term Project (University of Maryland)
#   James Parker
#   Kristopher Micinski
#   Matthew Mauriello
#


# IRanges type; poorly approximated - not sure what needs to be in here.
type IRanges
    start::Int
    finish::Int
    width::Int
    names::Int
end

# GRanges type; poorly approximated
# TO DO: Write simple constructor? Seperate into metadata and GRange data using matrix structure.
type GRanges
    seqnames::Array{ASCIIString}
    ranges::Array{IRanges}
    strand::Array{Char}
    score::Array{Int} #This is suppose to be metadata
    GC::Array{Float64} #This is suppoe to be meatadata
end

# Show the structure and data contained within the GRanges object; decently approximated
# This is all assigned directly at initialization; however, you could imagine a constructor that
# automated some of these fields, maybe?
function show(gr::GRanges)
    println(string("\tGRanges with ", length(gr.seqnames)," ranges and ", 2, " metadata columns:"))
    println("\tseqnames \tranges \t\tstrand \t | \tscore \tGC")
    println("\t<ASCIIString> \t<IRanges> \t<Char> \t | \t<Int> \t<Float64>")
    for i = 1:length(gr.seqnames)
        println(string("\t", gr.seqnames[i], "\t\t", "[",gr.ranges[i].start, " ", gr.ranges[i].finish,"]", "\t\t", gr.strand[i], "\t | \t", gr.score[i], "\t", gr.GC[i]))
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
end

function strand(gr::GRanges)
end

# Right now I have poorly inserted what is termed metadata at the top; when writing this function we'll have to define some kind of structure and matrix
# for holding metadata. (Impact: GRanges, show, and probably most other functions if this isn't handled soon
function mcols(gr::GRanges)
end

function names(gr::GRanges)
end

# length renamed to len as to not conflict with Julia default length function? Perhaps, just use the Length()
function len(gr::GRanges)
end

# Example usage:
gr = GRanges(["chr1", "chr2", "chr3"], [IRanges(0,0,0,0), IRanges(0,0,0,0), IRanges(0,0,0,0)] , ['+','-','*'] , [1,2,3],  [0.0,0.0,0.0])

