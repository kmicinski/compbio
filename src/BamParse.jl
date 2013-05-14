# Parse BAM files in Julia.
#
# This file allows a user to parse a BAM alignment file and reify the
# results into a Julia GRanges structure.

module BamParse

export BamFile,
       open_bam,
       close_bam,
       open_close

libbam_so = "libbam"

# Abstract type representing bam file
type BamFile
    file::Ptr
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

function open_close(fn :: String)
    handle = open_bam(fn)
    close_bam(handle)
end

end