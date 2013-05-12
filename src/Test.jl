using BioRanges
gr = GRanges([GRange("chr1", IRange(1,7,"a"), '+', 0, 0.0),
GRange("chr1", IRange(1,7,"d"), '+', 0, 0.0),
GRange("chr2", IRange(8,20,"b"), '-', 1, 0.1),
GRange("chr2", IRange(8,20,"b"), '-', 1, 0.1),
GRange("chr2", IRange(8,20,"b"), '-', 1, 0.1),
GRange("chr1", IRange(1,7,"d"), '+', 0, 0.0),
GRange("chr1", IRange(10,42,"c"), '+', 0, 0.2),
GRange("chr3", IRange(1,7,"d"), '+', 0, 0.0),
GRange("chr3", IRange(1,7,"d"), '+', 0, 0.0),
GRange("chr3", IRange(1,7,"d"), '+', 0, 0.0),
GRange("chr3", IRange(1,7,"d"), '+', 0, 0.0)
])
