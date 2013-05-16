n = 1000000
set.seed(63)
starts = sample( 1:1000, n, replace = T)
finishes = starts + sample( 190:210, n, replace = T)
intervals = IRanges( starts, finishes)

tree = IntervalTree(intervals)

options("digits.secs")
time = function(f, n) {
  times = vector(length=n)
  for (i in 1:n) {
    start = Sys.time()
    f()
    end = Sys.time()
    times[i] = end - start
  }
  print(times)
  print(mean(times))
  print(var(times))
}

tree = IntervalTree(intervals)
construct = function() {
  tree = IntervalTree(intervals)
}
query1 = function() {
  query = IRanges(c(0),c(1000))
  countOverlaps(query,tree)
}
query2 = function() {
  query = IRanges(c(100),c(300))
  countOverlaps(query,tree)
}
query3 = function() {
  query = IRanges(c(650),c(850))
  countOverlaps(query,tree)
}
time(construct,20)
time(query1,20)
time(query2,20)
time(query3,20)