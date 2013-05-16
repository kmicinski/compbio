# Parse BAM files in Julia.
#
# This file allows a user to parse a BAM alignment file and reify the
# results into a Julia GRanges structure.  
# 

module BamParse

import Base.*
import IRanges

export BamFile,
       open_bam,
       close_bam,
       open_close

libbam_so = "libbam"
simplebam_so = "simplebam"

# Abstract type representing bam file
type BamFile
    file::Ptr
end

# Type representing BAM parsing parameters.  This is mirrored in
# Rsamtools 
# 
# To parse a BAM file we need simply a list of ranges to parse.
type BamParseParams
    regions :: Vector{IRanges}
end

function open_bam(filename :: String) 
    handle = ccall((:samopen,"libbam"),
                   Ptr{Void}, (Ptr{Uint8},Ptr{Uint8},Ptr), 
                   filename, "rb", C_NULL);

    # Wrap it up and actually put it in a BamFile structure.
    BamFile(handle)
end

function close_bam(file :: BamFile)
    ccall((:samclose,"libbam"), Void, (Ptr{Void},), file.file)
end

# Parse a BAM file, given a set of scan parameters.  Use the
# simplified wrapper to SAMtools to call into it and use its parsing
# utilities.
function parse_bam(file :: BamFile, params :: BamScanParams)
    ccall((:parse_bam,"simplebam"), Void, (Ptr{Void},), file.file, start(params), finish(params))
end

function open_close(fn :: String)
    handle = open_bam(fn)
    close_bam(handle)
end

end