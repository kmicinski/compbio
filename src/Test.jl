using BioRanges
gr = GRanges([
GRange("chr1", IRange(1,7,"a"), '+', 0, 0.0),
GRange("chr2", IRange(2,8,"b"), '-', 1, 0.1),
GRange("chr2", IRange(3,9,"c"), '-', 1, 0.1),
GRange("chr2", IRange(4,10,"d"), '-', 1, 0.1),
GRange("chr1", IRange(5,11,"e"), '+', 0, 0.0),
GRange("chr1", IRange(6,12,"f"), '+', 0, 0.2),
GRange("chr3", IRange(7,13,"g"), '+', 0, 0.0),
GRange("chr3", IRange(8,14,"h"), '+', 0, 0.0),
GRange("chr3", IRange(9,15,"i"), '+', 0, 0.0),
GRange("chr3", IRange(10,16,"j"), '+', 0, 0.0)
])
