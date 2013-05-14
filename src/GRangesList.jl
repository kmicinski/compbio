# GRangesList Object for Julia BioSeq
#   Modeled after original BioConductor documentation
#   :http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
#
#   CMSC702, Spring 2013, Term Project (University of Maryland)
#   James Parker
#   Kristopher Micinski
#   Matthew Mauriello
#
#   This file adds the GRangesList Object, currently the only supported function allows GRanges to split into a GRangesList object.

export GRangesList, splitGRanges

type GRangesList
    grangeslist::Vector{GRanges}

    #Constructor
    GRangesList(grangeslist::Vector{GRanges}) = begin
        new(grangeslist)
    end
end

# [Utilized by GRanges]
# Split GRanges Object into several GRanges;
#   return GRangesList Object
function splitGRanges(gr::GRanges, each::Int32)
    if len(gr) <= each || each <= 0
        print("Split would have no effect")
        return
    end

    #Calculate the number of GRanges needed
    count::Int32
    count = floor((len(gr)/each))
    if (len(gr)%each) != 0
        count = count + ((len(gr)%each)/(len(gr)%each))
    end

    #Load Data
    repo = Array(GRanges, count)
    for i = 1:count
        index = (each * (i - 1))+1
        if  (len(gr) - index) >= each
            boundary = (index+(each-1))
        else
            boundary = len(gr)
        end
        repo[i] = GRanges([gr.granges[index:boundary]])
    end
    return GRangesList(repo)

end
