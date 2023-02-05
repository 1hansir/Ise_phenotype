nlab = 498
basemax = 4
Cmax = 25

source('./set_mdir.R')
dat = read.csv(paste0(mdir,
                      "/Real/T2D/T2D_lab", nlab[1],
                      "base",basemax,
                      "max",Cmax,
                      ".csv"))

# Create longitudinal follow up data
dat$I[is.na(dat$I)] = 0
dat$T[is.na(dat$T)] = max(dat$C)+1
dat$S[is.na(dat$S)] = max(dat$C)+1
dat$rY0 = ifelse(dat$I==0,
                 dat$C+1,
                 dat$T)
dat$rY1 = ifelse(dat$I==0,
                 0,
                 dat$C+1-dat$T)   # 是不是反了？

longdat = data.frame(
  ID =rep(dat$ID, dat$C+1), 
  T = unlist(sapply(dat$C, function(x) 0:x)),
  Y = unlist(mapply(function(x,y) rep(0:1,c(x,y)), 
                    x=dat$rY0,y=dat$rY1 ))
)

longdat = cbind(longdat,
                dat[rep(1:nrow(dat),dat$C+1),
                    5:(ncol(dat)-2)])

rownames(longdat) = 1:nrow(longdat)
row.lab = sum(dat$C[1:nlab]+1)

longlab = cbind(1:row.lab,longdat[1:row.lab,])
colnames(longlab)[1]='X'

longrest = cbind(1:as.numeric(length(longdat[,1])-row.lab),longdat[-1:-row.lab,])
colnames(longrest)[1]='X'

write.csv(longlab, file = paste0(mdir,
                                 "/Real/T2D/T2D_labeled.csv"))
write.csv(longrest, file = paste0(mdir,
                                 "/Real/T2D/T2D_unlabeled.csv"))
print("Creating Successfully")