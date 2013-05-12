# GRanges Class for Julia BioSeq
#   Modeled after original BioConductor documentation
#   :http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
#
#   CMSC702, Spring 2013, Term Project (University of Maryland)
#   James Parker
#   Kristopher Micinski
#   Matthew Mauriello
#

export GRange, GRanges, show, seqnames


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
    println(string("\tGRanges with ", length(gr.granges)," ranges and ", 2, " metadata columns:"))
    println("\t\tseqnames \tranges \t\tstrand \t | \tscore \tGC")
    println("\t\t<ASCIIString> \t<IRanges> \t<Char> \t | \t<Int> \t<Float64>")
    for i = 1:length(gr.granges)
        println(string("\t", "\t", gr.granges[i].seqname, "\t\t", "[",gr.granges[i].range.start, " ", gr.granges[i].range.finish,"]", "\t\t", gr.granges[i].strand, "\t | \t", gr.granges[i].score, "\t", gr.granges[i].GC))
    end
    println("\t---")
    print("\tseqlengths:")
    print("\n\t")
    for i = 1:length(gr.granges)
        print(string(gr.granges[i].seqname, "\t"))
    end
    print("\n\t")
    for i = 1:length(gr.granges)
        # TO DO: this should be the total length of the sequence?
        print(string("NA", "\t"))
    end
    print("\n")
end

# seqnames appears to be just another view of the data
# TO DO: figure out what runs are. Sounds silly, I know.
function seqnames(gr::GRanges)

    runs = 0
    cur = ""
    seq = Array(String, length(gr.granges))
    for i = 1:length(gr.granges)
        if seq == ""
            cur = gr.granges[i].seqname
            runs = runs + 1
        elseif cur != gr.granges[i].seqname
            cur = gr.granges[i].seqname
            runs = runs + 1
        end
        seq[i] = gr.granges[i].seqname
    end
    uniqueElems = sort(unique(seq))


    println(string("\t", typeof(gr.granges), " of length ", length(gr.granges), " with ", runs," runs"))

    print("\tLengths:\t")
    count = 0
    seq = ""
    for i = 1:length(gr.granges)
        if seq == ""
            seq = gr.granges[i].seqname
            count = 1
        elseif seq != gr.granges[i].seqname
            print(string(count,"\t"))
            count = 1
            seq = gr.granges[i].seqname
        else
            count = count + 1
        end
    end
    println(count)


    print("\tValues:\t\t")
    seq = ""
    for i = 1:length(gr.granges)
        if seq == ""
            seq = gr.granges[i].seqname
        elseif seq != gr.granges[i].seqname
            print(string(seq,"\t"))
            seq = gr.granges[i].seqname
        end
    end
    println(seq)

    print(string("\tLevels(", length(uniqueElems) ,"):\t"))
    for i = 1:length(uniqueElems)
        # List of unique sequences?
        print(string(uniqueElems[i],"\t"))
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

