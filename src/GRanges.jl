
# GRanges Class for Julia BioSeq
#   Modeled after original BioConductor documentation
#   :http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
#
#   CMSC702, Spring 2013, Term Project (University of Maryland)
#   James Parker
#   Kristopher Micinski
#   Matthew Mauriello
#
#   Example Construction:
#   gr = GRanges(["a", "b", "c"],["chr1", "chr2", "chr3"], [IRanges(1,7,7-1,0), IRanges(0,0,0,0), IRanges(0,0,0,0)] , ['+','-','*'] , [1,2,3],  [0.0,0.0,0.0])
#
#   Limitations:
#   Parallel Operations in R are not supported by this package. See "TO DO"s.
#
#   TO DO: Investigate Gapped Alignments.

export GRanges,
       show,
       seqnames,
       ranges,
       strand,
       mcols,
       names,
       len,
       mergeGRanges,
       head,
       rep,
       rev,
       tail,
       window,
       seqselect,
       start,
       finish,
       width,
       range,
       flank,
       shift,
       resize,
       union

# GRanges type; poorly approximated TO DO: Write simple constructor?
# Seperate into metadata and GRange data using matrix structure.  TO
# DO: Change 'type' to 'Immutable'?

type GRanges
    granges::Vector{GRange}

    #Constructor
    GRanges(granges::Vector{GRange}) = begin
        new(granges)
    end
end

# Show the structure and data contained within the GRanges object;
# decently approximated This is all assigned directly at
# initialization; however, you could imagine a constructor that
# automated some of these fields, maybe?
function show(gr::GRanges)
    println(string("\tGRanges with ", length(gr.granges)," ranges and ", 2, " metadata columns:"))
    println("\t\tseqnames \tranges \t\tstrand \t | \tscore \tGC")
    println("\t\t<String> \t<IRange> \t<Char> \t | \t<Int> \t<Float64>")
    for i = 1:length(gr.granges)
        println(string("\t", "\t", gr.granges[i].seqname, "\t\t"
                       , "[",gr.granges[i].range.start, " "
                       , gr.granges[i].range.finish,"]", "\t\t"
                       , gr.granges[i].strand, "\t | \t", gr.granges[i].score
                       , "\t", gr.granges[i].GC))
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
    (uniqueElems, runs) = call_unique_seq(gr)

    println(string("\t", typeof(gr.granges), " of length "
                   , length(gr.granges), " with ", runs," runs"))

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
    println(length(gr.granges));
    println("\t\tstart\tend\twidth\tnames");
    sample = Array(IRange, len(gr))
    for i = 1:length(gr.granges)
        println(string("\t[", i, "]\t",gr.granges[i].range.start, "\t", gr.granges[i].range.finish, "\t", gr.granges[i].range.width, "\t", gr.granges[i].range.name));
        sample[i] = gr.granges[i].range
    end
    return sample
end

function strand(gr::GRanges)
    (uniqueElems, runs) = call_unique_strand(gr)
    print("\tfactor-Rle of length: ")
    println(string(length(gr.granges), " with ", runs , " runs"))

    print("\tLengths:\t")
    count = 0
    seq = ""
    for i = 1:length(gr.granges)
        if seq == ""
            seq = gr.granges[i].strand
            count = 1
        elseif seq != gr.granges[i].strand
            print(string(count,"\t"))
            count = 1
            seq = gr.granges[i].strand
        else
            count = count + 1
        end
    end
    println(count)

    print("\tValues:\t\t")
    seq = ""
    for i = 1:length(gr.granges)
        if seq == ""
            seq = gr.granges[i].strand
        elseif seq != gr.granges[i].strand
            print(string(seq,"\t"))
            seq = gr.granges[i].strand
        end
    end
    println(seq)

    print(string("\tLevels(", length(uniqueElems) ,"):\t"))
    for i = 1:length(uniqueElems)
        # List of unique sequences
        print(string(uniqueElems[i],"\t"))
    end
    print("\n")

end

# Right now I have poorly inserted what is termed metadata at the top;
# when writing this function we'll have to define some kind of
# structure and matrix for holding metadata. (Impact: GRanges, show,
# and probably most other functions if this isn't handled soon
function mcols(gr::GRanges)
    print("\tDataFrame with X rows and Y columns")
    print("\n\t\tscore\t\t")
    println("GC")
    println("\t\t<integer>\t<numeric>")

    for i = 1:length(gr.granges)
        println(string("\t", i, "\t\t", gr.granges[i].score, "\t", gr.granges[i].GC))
    end
end


function seqlengths(gr::GRanges)
end

function names(gr::GRanges)

    print("\t[1] ")
    for i = 1:length(gr.granges)
        print(string(gr.granges[i].range.name, " "))
    end
    print("\n")

end

# length renamed to len as to not conflict with Julia default length function? Perhaps, just use the Length()
function len(gr::GRanges)
    return length(gr.granges)
end

function mergeGRanges(gr::GRanges, br::GRanges)
    array = Array(GRange, (len(gr) + len(br)))
    for i = 1:(len(gr) + len(br))
        if (i > len(gr))
            array[i] = br.granges[(i-len(gr))]
        else
            array[i] = gr.granges[i]
        end
    end
    return GRanges(array)

end

#TO DO: Define Subsetting method

#TO DO: Define seqselect method

function head(gr::GRanges, n::Int)
    if (n > len(gr))
        n = len(gr)
    end

    sample = gr.granges[1:n]
    ret = GRanges(sample)
    show(ret)
    return ret
end

function rep(gr::GRanges, times::Int)
    if len(gr) <= 0 || times <= 0
        print("Input contains zero")
        return
    end

    sample = Array(GRange, len(gr)*times)

    for i = 1:times
        for j = 1:len(gr)
            sample[(j + ((i-1) * len(gr)))] = gr.granges[j]
        end
    end

    range = GRanges(sample)
    return range
end

function rev(gr::GRanges)
    sample = GRanges(reverse(gr.granges))
    return sample
end

function tail(gr::GRanges, n::Int)
    if (n > len(gr))
        n = len(gr)
    end

    sample = gr.granges[(len(gr)-n+1):len(gr)]
    show(GRanges(sample))
end

function window(gr::GRanges, start::Int, finish::Int)
    if (start > finish || start <= 0 || finish >= len(gr) || start >= len(gr) || finish <= 0)
        print("Invalid Range")
    end

    sample = gr.granges[start:finish]
    show(GRanges(sample))
end

function seqselect(gr::GRanges, start::Array{Int32}, finish::Array{Int32})
    if (size(start)[1] == 1 
        && size(start)[2] == 2
        && size(finish)[1] == 1
        && size(finish)[2] == 2)
        if (
            start[1, 1] < start[1,2]
            && finish[1,1] < finish[1,2]
            && start[1,2] <= len(gr)
            && start[1,1] > 0
            && finish[1,2] <= len(gr)
            && finish[1,2] > 0)
            count = 0
            for i = 1:len(gr)
                if ((gr.granges[i].range.start >= start[1,1] 
                     && gr.granges[i].range.start <= start[1,2])
                    || 
                    (gr.granges[i].range.start >= finish[1,1]
                     && gr.granges[i].range.start <= finish[1,2]))
                    count = count + 1
                end
            end

            if (count > 0)
                sample = Array(GRange, count)
                count = 1
                for i = 1:len(gr)
                    if (
                        (gr.granges[i].range.start >= start[1,1]
                         && gr.granges[i].range.start <= start[1,2])
                        || (gr.granges[i].range.start >= finish[1,1]
                            && gr.granges[i].range.start <= finish[1,2]))
                        sample[count] = gr.granges[i]
                        count = count + 1
                    end
                end
                found = GRanges(sample)
                show(found)
                return found
            else
                print("No Ranges Found")
            end

        else
            println("Ranges in start or finish are out of bounds.")
        end
    else
        println("Start/Finish should be 1x2 tuples representing range of display.")
    end

    return
end

function start(gr::GRanges)
    starts = Array(Int, len(gr))
    for i = 1:len(gr)
        starts[i] = gr.granges[i].range.start
    end
    return starts
end

function finish(gr::GRanges)
    finishes = Array(Int, len(gr))
    for i = 1:len(gr)
        finishes[i] = gr.granges[i].range.finish
    end
    return finishes
end

function width(gr::GRanges)
    widths = Array(Int, len(gr))
    for i = 1:len(gr)
        widths[i] = gr.granges[i].range.width
    end
    return widths
end

function range(gr::GRanges)
    println(string("\tGRanges with ", length(gr.granges)," ranges and ", 0, " metadata columns:"))
    println("\t\tseqnames \tranges \t\tstrand")
    println("\t\t<String> \t<IRange> \t<Char>")
    for i = 1:length(gr.granges)
        println(string("\t", "\t", gr.granges[i].seqname, "\t\t", "[",gr.granges[i].range.start, " ", gr.granges[i].range.finish,"]", "\t\t", gr.granges[i].strand))
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

function flank(gr::GRanges, n::Int, start::Bool)
    print("Function Incomplete")
end

function shift(gr::GRanges, n::Int)
    print("Function Incomplete")
end

# Resize a GRanges object.
function resize(gr::GRanges, n::Int)
    print("Function Incomplete")
end

function reduce(gr::GRanges)
    print("Function Incomplete")
end

function gaps(gr::GRanges)
    print("Function Incomplete")
end

function disjoin(gr::GRanges)
    print("Function Incomplete")
end

function converage(gr::GRanges)
    print("Function Incomplete")
end

# Create a union of two GRanges objects.
function union(gr::GRanges, br::GRanges)

    for i = 1:len(br)
        add = true
        for j = 1:len(gr)
            if (isequal(gr.granges[j].range, br.granges[i].range))
                add = false
                break
            end
        end
        if (add)
            gr = mergeGRanges(gr, GRanges([br.granges[i]]))
        end
    end

    return gr
end


function intersect(gr::GRanges, br::GRanges)

end

function setdiff(gr::GRanges, br::GRanges)
    print("Function Incomplete")
end

# Internal functions -- not exported.

function call_unique_seq(gr::GRanges)
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
    return uniqueElems, runs
end

function call_unique_strand(gr::GRanges)
    runs = 0
    cur = ""
    seq = Array(Char, length(gr.granges))
    for i = 1:length(gr.granges)
        if seq == ""
            cur = gr.granges[i].strand
            runs = runs + 1
        elseif cur != gr.granges[i].strand
            cur = gr.granges[i].strand
            runs = runs + 1
        end
        seq[i] = gr.granges[i].strand
    end
    uniqueElems = sort(unique(seq))
    return uniqueElems, runs
end
