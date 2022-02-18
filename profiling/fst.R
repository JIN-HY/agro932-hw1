geno <- read.table("largedata/geno.txt", header=FALSE)
names(geno)[1:4] <- c("chr", "pos", "ref", "alt")

geno[, 5:ncol(geno)] <- apply(geno[, 5:ncol(geno)], c(1,2), function(x) {sum(as.numeric(strsplit(x, split = "/")[[1]]))})

#geno$p <- apply(geno[, 5:ncol(geno)], 1, function(x) {prop.table(table(x))})


geno$p <- apply(geno[, 5:ncol(geno)], 1, function(x) {sum(is.na(x))})
geno$pmissing <- apply(geno[, 5:(ncol(geno)-1)], 1, function(x) {sum(is.na(x))})
geno$p <- geno$p/ (2*(ncol(geno)-6-geno$pmissing))

geno$p1 <- apply(geno[, 5:200], 1, function(x) {sum(is.na(x))})
geno$p1missing <- apply(geno[, 5:200], 1, function(x) {sum(is.na(x))})
geno$p1 <- geno$p1/ (2*196-geno$p1missing)

geno$p2 <- apply(geno[, 201:ncol(geno)], 1, function(x) {sum(is.na(x))})
geno$p2missing <- geno$pmissing - geno$p1missing
geno$p2 <- geno$p2/ (2*211 - geno$p2missing)


geno$fst <- with(geno, ((p1-p)^2 + (p2-p)^2)/(2*p*(1-p)) )

write.table(geno[,c(1:4,ncol(geno))], "cache/fst.tsv", sep="\t", row.names = FALSE, quote=FALSE)

fst <- read.csv("cache/fst.tsv", sep = "\t")

pdf("graphs/fst_persite.pdf")
plot(fst$pos, fst$fst, xlab="Physical position", ylab="Fst value", main="", xlim=c(2000000,3000000))
dev.off()

# theta calculation
pi <- function(n=10, p=0.1){
  return(n/(n-1)*(1-p^2-(1-p)^2))
}
geno$h1 <- 0
geno$h2 <- 0

for (row in 1:nrow(geno)) {
  geno$h1[row] <- pi(n=2*(ncol(geno)-8-geno$p1missing[row]), p=geno$p1[row])
  geno$h2[row] <- pi(n=2*(ncol(geno)-8-geno$p2missing[row]), p=geno$p2[row])
}
theta1 <- sum(geno$h1)/ 1e6
theta2 <- sum(geno$h2)/ 1e6



# plot by sliding window
win <- 10000
step <- 2000
start0 <- 2000000
end0 <- 3000000
start <- seq(start0, end0-win, by = step)
end <- start + win -1
fst_win <- data.frame(start = start, end = end, fst = 0, nsite = 0)

for (row in 1:nrow(fst)) {
  pos <- fst$pos[row]
  pos.fst <- fst$fst[row]
  for (r in 1:nrow(fst_win)) {
    if ((pos >= fst_win$start[r]) & (pos <= fst_win$end[r])) {
      fst_win$nsite[r] = fst_win$nsite + 1
      fst_win$fst[r] = fst_win$fst + pos.fst
    }
  }
}
fst_win$fst <- fst_win$fst/ fst_win$nsite

pdf("graphs/fst_bywindow.pdf")
plot(fst_win$start, fst_win$fst, type = "l", xlab="Physical position", ylab="Fst value", main="", xlim=c(2000000,3000000))
dev.off()
