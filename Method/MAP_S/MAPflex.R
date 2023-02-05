#' Internal function to do inverse logit transformation
#' @noRd
g.logit = function(xx=NULL){exp(xx)/(exp(xx)+1)}

#' Internal function to do logit transformation
#' @noRd
logit = function(xx=NULL){log(xx/(1-xx))}

#' Internal function to get the probabilities that match the prevalence
#' @noRd
Prob.Scale = function(pp=NULL, prev=NULL){
    logit_pp = logit(pp)
    logit_pp[logit_pp == -Inf] = min(logit_pp[logit_pp > -Inf], na.rm=TRUE)
    logit_pp[logit_pp ==  Inf] = max(logit_pp[logit_pp < Inf], na.rm=TRUE)
    cc = uniroot(function(xx){mean(g.logit(logit_pp-xx), na.rm=TRUE)-prev},interval=c(-10,10))$root
    g.logit(logit_pp-cc)
}

#' Internal function to get the summary table
#' @noRd
IDtab = function(data=NULL){
    mat = cbind(data$mat, data$note)
    colnames(mat)=c(colnames(data$mat),"note")
    df = matrix(c("NO","YES")[as.matrix(!is.na(mat))+1],ncol=ncol(mat))
    df = as.data.frame(df)
    colnames(df) = colnames(mat)
    as.data.frame(table(df))
}


#' Internal function to do MAP algorithm
#' data: a list with three objects (ID, mat, note)
#' yes.con: concomitant not implemented yet, not sure what to include
#' family: has to be specified 
#' for each variable, clustering patients with non-zero values and pateints with 
#' zero value will be assigned to have zero probability.
#' For each variable, patients with zero/NA note or NA value in this variable will be 
#' assigned to be NA probabilities
#' Kmeans: use log-transformed data.
#' full results: probabilities are not re-scaled by the prevalence
#' prevalence estimated based on filter-positive patients

#' @noRd
MAP_flex_core = function(data = NULL,  yes.con = FALSE, full.output = TRUE){

    #vname.log = paste(vname,"_log",sep="")

    #prior.all = NULL
    #post.all = NULL
    #group.all = NULL
    IDtab = IDtab(data)
    
    vname = colnames(data$mat)
    vname.log = paste(vname,"_log",sep="")
    name.all = c(vname,vname.log)
    data$mat = cbind(data$mat, log(1+data$mat))
    colnames(data$mat) = name.all
    data$note = log(data$note+1)
    
    totalN = length(data$ID)
    post.all = Matrix(0, nrow=totalN, ncol=ncol(data$mat))
    prob.all = Matrix(0, nrow=totalN, ncol=ncol(data$mat))
    family = c( rep("poisson", length(vname)), rep("gaussian", length(vname.log)) )
    
    for(i in seq_along(name.all)){
        tmpfm = as.formula(paste(name.all[i],"~1"))
        if(yes.con){
            #ii = i%%(length(vname))
            #if(ii==0){ii = length(vname)}
            #vname.log.out = setdiff(vname.log, vname.log[ii])
            #tmpfm2 = as.formula(paste("~", paste(c("note",vname.log.out),collapse="+")))
            #tmpfm2 = FLXPmultinom(tmpfm2)
            #na.ind = apply(dat,1,function(x) { any(is.na(x)) } )
            #dat.tmp = dat[!na.ind,]
            tmpfm2 = FLXPconstant()
        }else{
            tmpfm2 = FLXPconstant()
            zero.ind = is.element(data$mat[,i], c(0))
            #zero.ind = rep(F,dim(data$mat)[1])
            na.ind = (is.element(data$mat[,i],c(NA)) | is.element(data$note[,1],c(NA)))
            exclude.ind = (zero.ind | na.ind)
            dat.tmp = data.frame(data$mat[!exclude.ind,i],data$note[!exclude.ind,1])
            colnames(dat.tmp) = c(name.all[i], "note")
        }
        totalN00 = nrow(dat.tmp)
        
        set.seed(1)
        n.clust = 1
        iter = 0
        while(n.clust < 2 & iter < 5){
            tmpfit = flexmix(tmpfm, k = 2,
                             model = FLXMRglmfix(fixed =~note,varFix=FALSE, family=family[i]),
                             concomitant=tmpfm2,control=list(tol = 1e-8), data=dat.tmp)
            n.clust = length(unique(tmpfit@cluster))
            iter = iter+1
        }

        ##if(name.all[i]=="nlp_log"){browser()}
        if(n.clust>1){
            avgcount = dat.tmp[,name.all[i]]
            tmpdiff = mean(avgcount[tmpfit@cluster==2]) - mean(avgcount[tmpfit@cluster==1])
            tmpind =  as.numeric((tmpdiff > 0) + 1)
            qprobtmp = qnorm(posterior(tmpfit))
            qprob = qprobtmp[,tmpind]
            qprob[is.infinite(qprob)] = -1*qprobtmp[is.infinite(qprob),3-tmpind]
            ### deal with some extreme clustering results ###
            if(sum(!is.infinite(qprob))>=2){
                qprob[qprob == Inf] = max(qprob[!is.infinite(qprob)])
                qprob[qprob == -Inf] = min(qprob[!is.infinite(qprob)])
            }else if(sum(!is.infinite(qprob))==1){
                if(qprob[!is.infinite(qprob)] >= 0){
                    qprob[qprob == Inf] = qprob[!is.infinite(qprob)]
                    qprob[qprob == -Inf] = qnorm(1/totalN00)
                }else{
                    qprob[qprob == Inf] = qnorm(1-1/totalN00)
                    qprob[qprob == -Inf] = qprob[!is.infinite(qprob)]
                }
            }else{
                qprob[qprob == Inf] = qnorm(1-1/totalN00)
                qprob[qprob == -Inf] = qnorm(1/totalN00)
            }
            #############
            # these details seem not useful, commented out for now #
            
            # qpost.tmp = data.frame("qprob" = qprob,
            #                    "group" = cbind(2-tmpfit@cluster,tmpfit@cluster-1)[,tmpind],
            #                    "prior" = tmpfit@prior[tmpind])
            # 
            # qpost.tmp = as.matrix(qpost.tmp)
            # qpost = Matrix(0, nrow=totalN, ncol=3)
            # colnames(qpost) = colnames(qpost.tmp)
            # qpost[!na.ind,] = qpost.tmp
            ##
            
            
        }else{
            qprob = rep( qnorm(1-1/totalN00), nrow(dat.tmp))
            # qpost = data.frame("qprob" = rep( qnorm(1-1/nrow(dat)), nrow(dat)),
            #                    "group" = rep(1,nrow(dat)),
            #                    "prior" = 1)
            # qpost = as.matrix(qpost)
            # qpost[na.ind, ] = NA
            warning(paste("The clustering step does not converge for variable ",
                          name.all[i], "!", sep="") )
        }
        
        post.all[!exclude.ind,i] = qprob
        if(sum(na.ind)>0){
          post.all[na.ind,i] = NA
        }
        prob.all[!exclude.ind,i] = pnorm(qprob)
        if(sum(na.ind)>0){
          prob.all[na.ind,i] = NA
        }
        
        
        #post.all = cbind(post.all,qpost[,"qprob"])
        #prior.all = c(prior.all,qpost[1,"prior"])
        #post.all = cbind(post.all,qpost[,"qprob"])
        #group.all = cbind(group.all, qpost[,"group"])
        
    }
    
    colnames(post.all) = name.all
    colnames(prob.all) = name.all 
    
    ## select eligible patients to estimate the prevalence 
    ## and re-scale predicted prob.
    prob.all00 = prob.all
    na.id = (rowSums(is.na(prob.all00))==ncol(prob.all00))
    prob.all00[is.na(prob.all00)] = 0
    mu00 = rowMeans(prob.all00)
    exclude.ind = (mu00==0)
    post.all00 = as.matrix(post.all[!exclude.ind,])
    mat00 = data$mat[!exclude.ind,]
    rm(post.all)
    
    final.score00 = (rowMeans(pnorm(post.all00),na.rm=TRUE) + 
                       pnorm(rowMeans(post.all00,na.rm=TRUE)))/2
    final.score00[is.na(final.score00)] = NA

    avgcount = rowMeans(mat00[,vname.log,drop=FALSE], na.rm=TRUE)
    avgpost = rowMeans(post.all00[,vname.log,drop=FALSE], na.rm=TRUE)
    avgcount = avgcount[!is.na(avgpost)]
    avgpost = avgpost[!is.na(avgpost)]

    cluster.ori = kmeans(cbind(avgcount,avgpost),centers=2)
    #cluster.ori = kmeans(cbind(dat[,vname.log],post.all[,vname.log]),centers=2) ### standardize?
    class.pos = 1 + 1*(cor(cluster.ori$cluster==2, avgcount)>0)
    prev.est0 = mean(cluster.ori$cluster==class.pos)
    prev.est = (mean(final.score00,na.rm=TRUE)+prev.est0)/2

    final.score00 = Prob.Scale(final.score00, prev.est)
    cut.MAP = as.numeric(quantile(final.score00,1-prev.est,na.rm = TRUE))
    final.score = Matrix(0,nrow=totalN,ncol=1)
    final.score[!exclude.ind,1] = final.score00
    if(sum(na.id)>0){
      final.score[na.id,1] = NA
    }
    
    
    
    IDtab = IDtab[IDtab$Freq>0,]
    cat("####################### \n")
    cat("MAP only considers pateints who have note count data and
        at least one nonmissing variable!\n")
    cat("####\nHere is a summary of the input data:\n")
    cat("Total number of unique IDs:", sum(IDtab$Freq), "\n")
    print(IDtab)
    cat("#### \n")
    
    
    if(full.output){
        list("scores" = final.score, "cut.MAP" = cut.MAP, "scores.all" = prob.all)
    }else{
        list("scores" = final.score, "cut.MAP" = cut.MAP)
    }

}


