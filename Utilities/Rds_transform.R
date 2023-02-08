nlab = 498
basemax = 4
Cmax = 25

source('./set_mdir.R')
dat = read.csv(paste0(mdir,
                      "/Real/T2D/T2D_lab", nlab[1],
                      "base",basemax,
                      "max",Cmax,
                      ".csv"))

dataname = paste0(mdir, "/Real/T2D/T2D.Rds")

dat$I[is.na(dat$I)] = 0
dat$T[is.na(dat$T)] = max(dat$C)+5
dat$S[is.na(dat$S)] = max(dat$C)+1
X = data.frame(cbind(dat$X.1,dat$X.2,dat$X.3,dat$X.4,dat$X.5,dat$X.6,dat$X.7,dat$X.8,dat$X.9,dat$X.10,dat$X.11,dat$X.12))

saveRDS(list(X=X,T=dat$T,C=dat$C,S=dat$S,I=dat$I,
               #P=P,
               MAP=dat$MAP),file = dataname)

dat_read = readRDS(dataname)
print(length(dat_read$C))