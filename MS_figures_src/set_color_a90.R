
b.colors <- brewer.pal(7,'Dark2') # 'Dark2', 'Set2', 'Paired'

a.colors <- c(b.colors,b.colors)

s.hm <- s.hb <- s.hr <- 5 # diamond
s.mb <- 0 # square
s.mr <- 1 # circle
s.hc <- s.mc <- 3 # plus
# 4 -- 'x'

## symbols

## a.sym <- c(s.mr, s.hm, s.hr, s.hb, s.mb, s.hc, s.mc)
vrt.sym <- rep(3,5)
bct.sym <- rep(4,2)
yst.sym <- rep(1,2)
a.sym <- c(vrt.sym, bct.sym, yst.sym)

## 3=plus, 4=x, 1=o
mv.sym <- c(3, 4, 1)
vby.sym <- c(4, 3, 1)

tax.sym=c('rod'=0,'mam'=1, 'pln'=1, 'vrt'=5, 'yst'=6,'bct'=2, 'arch'=5)
