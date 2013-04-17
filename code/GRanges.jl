type IRanges
    start::Int
    finish::Int
    width::Int
    names::Int
end

type GRanges
    seqnames::Array{ASCIIString}
    ranges::Array{IRanges}
    strand::Array{Char}
    score::Array{Int}
    GC::Array{Float64}
end

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
        print(string("NA", "\t"))
    end
    print("\n")
end

function seqnames(gr::GRanges)
end

function ranges(gr::GRanges)
end

function strand(gr::GRanges)
end

function mcols(gr::GRanges)
end

function names(gr::GRanges)
end

function len(gr::GRanges)
end


gr = GRanges(["chr1", "chr2", "chr3"], [IRanges(0,0,0,0), IRanges(0,0,0,0), IRanges(0,0,0,0)] , ['+','-','*'] , [1,2,3],  [0.0,0.0,0.0])

