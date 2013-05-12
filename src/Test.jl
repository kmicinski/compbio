using BioRanges
gr = GRanges([
GRange("chr1", IRange(1,7,"a"), '-', 1, 1.0000000),
GRange("chr2", IRange(2,8,"b"), '+', 2, 0.8888889),
GRange("chr2", IRange(3,9,"c"), '+', 3, 0.7777778),
GRange("chr2", IRange(4,10,"d"), '*', 4, 0.6666667),
GRange("chr1", IRange(5,11,"e"), '*', 5, 0.5555556),
GRange("chr1", IRange(6,12,"f"), '+', 6, 0.4444444),
GRange("chr3", IRange(7,13,"g"), '+', 7, 0.3333333),
GRange("chr3", IRange(8,14,"h"), '+', 8, 0.2222222),
GRange("chr3", IRange(9,15,"i"), '-', 9, 0.1111111),
GRange("chr3", IRange(10,16,"j"), '-', 10, 0.0000000)
])
