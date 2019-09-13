### corresponding interal eRM objects
# # https://github.com/cran/eRm/blob/master/R/

#' @importFrom stats median model.matrix na.exclude pchisq
#' @import stats graphics methods lattice grDevices
#' @importFrom splines interpSpline
#' @importFrom MASS mvrnorm
#' @importFrom Matrix bdiag

`cmlprep` <-
  function(X01,mt_vek,mpoints,Groups,W,gmemb)
  {

    levs <- (gmemb-1)*max(Groups)+Groups              #merge Groups and gmemb vector into level vector

    if (length(Groups)==1) {                          #if no group contrast
      x_mt <- colSums(X01,na.rm=TRUE)                #item category raw scores as vector
      #eventuell x_mtlist auf NA gruppen aufbrechen
      x_mtlist <- list(x_mt)
      ngroups <- 1

    } else {                                            #if groups defined
      ngroups <- max(Groups)                            #number of groups
      x_mtlist <- by(X01,levs,colSums,na.rm=TRUE)       #item-category raw scores for each group (as list)
      x_mtlist.G <- by(X01,Groups,colSums,na.rm=TRUE)   #item-category raw scores for each group (as list)
      #FIXME!!! use x_mtlist??
      x_mt <- as.vector(unlist(x_mtlist.G))             #as vector: g1|g2|...
    }

    end1 <- length(mt_vek)*mpoints*ngroups
    mt_ind <- rep(1:end1,rep(mt_vek,mpoints*ngroups)) #category index vector (for converting x_mt into list)
    x_tmt <- split(x_mt,mt_ind)                       #list for likelihood: item-wise * ngroups
    rtot <- sum(mt_vek)*mpoints

    ics <-  rep(sequence(mt_vek),mpoints)                 #item category scores for each item as vector
    rv <- apply(X01,1,function(x) {                       #person raw scores of 0/1 matrix
      ics[!is.na(x)]%*%na.exclude(x)})

    #--- preparing index vector for item parameters ---
    if (ngroups > 1) {                                    #groups
      seglen <- sum(mt_vek)                               #length of beta vector (segment)
      gind <- rep(rep(1:ngroups,rep(seglen,ngroups)),mpoints) #parameter index vector for group extraction
    } else {
      gind <- rep(1,dim(W)[1])
    }

    #--- preparing lists for person splits ---
    rvlist <- split(rv,levs)                    #split person raw scores due to levels (NAgroup AND treatment)
    nrlist <- lapply(rvlist,function(rvel) {    #list with item raw score frequencies for each group (transposed)
      rvtab <- table(rvel)                            #raw score frequencies
      dnamevek <- as.numeric(unlist(dimnames(rvtab))) #different raw scores for 0 fill up
      nr <- rep (0,rtot+1)                            #setting 0 raw score frequencies
      nr[dnamevek+1] <- rvtab #vector with person raw scores from 1:rtot (with 0 fill up)
      nr <- nr[-1]
      return(nr)
    })


    if ((ngroups > 1) && (length(unique(gmemb)) > 1)) {          #NA groups AND Groups
      gg <- table(Groups,gmemb)
      #gg[gg > 0] <- 1
      g_NA <- as.vector(rowSums(gg))                             #How many NA-sub groups in each Group
      #grgm <- cbind(Groups, gmemb)
      #grgmst <- apply(grgm,1,function(x) {                       #merge indexes to characters
      #            paste(x[1],x[2]) })
      #GGind <- rank(unique(grgmst))
      #levtab <- table(levs)                                      #frequencies of levels
      #FIXME!!! following line wrong index
      #gby <- rep(GGind,levtab)                                   #ordering by NAgroups nested in Group

      #this probably does the job
      gby <- levs
    } else {
      g_NA <- 1
      gby <- gmemb
    }

    NAstruc <- by(!is.na(X01),gby,function(x) {                  #list of unique NA structures for each Group
      x.u <- unique(x)
      as.numeric(as.matrix(x.u))}) #NA's are coded with 0

    NAcheck <- sapply(NAstruc,sum)                         #if for certain NAgroups only 1 item was presented

    list(x_mt=x_mt,mt_ind=mt_ind,x_tmt=x_tmt,rtot=rtot,nrlist=nrlist,gind=gind,x_mtlist=x_mtlist,
         NAstruc=NAstruc,g_NA=g_NA,gby=gby)
  }

#### datcheck.LRtest
datcheck.LRtest <- function(x, X, model)
{
  #sanity checks for LRtest (internal function of LRtest.R)
  #x...submatrix (splitted with "splitcr" and called within Xlist)
  #X...original data matrix (from model fit)

  exclude <- NULL                                             #vector with items to be excluded

  #----check full/0 responses------
  n.NA <- colSums(apply(X,2,is.na))                                   #number of NA's per column
  maxri <- (dim(X)[1]*(apply(X,2,max,na.rm=TRUE)))-n.NA               #maximum item raw scores with NA
  ri <- apply(x,2,sum,na.rm=TRUE)                              #item raw scores
  exclude <- c(exclude,which((ri==maxri) | (ri==0)))

  #----check full(-1) NA's---------
  allna.vec <- apply(x,2,function(y) {
    naTF <- is.na(y)
    (sum(naTF) >= length(y-1))
  })
  exclude <- c(exclude,which(allna.vec))

  #----minimum category = 0--------
  ri.min <- apply(x,2,min,na.rm=TRUE)                                 #if no 0 responses
  exclude <- c(exclude,which(ri.min!=0))

  #----RSM-checks for same number of categories--------
  if ((model == "RSM") || (model == "LRSM")) {
    highcat <- max(X, na.rm=TRUE)                    #highest category in original data
    highcat.sub <- apply(x,2,max,na.rm=TRUE)             #RSM check for equal number of categories
    exclude <- c(exclude,which(highcat.sub != highcat))
  }

  #---PCM checks for all categories responses---------
  if ((model=="PCM") || (model=="LPCM")) {                         #check if there are missing categories for PCM (for RSM doesn't matter)
    cat.data <- apply(X,2,function(y) list(unique(na.exclude(y)))) #categories of orginal data
    cat.sub <- apply(x,2,function(y) list(unique(na.exclude(y))))  #categories of subgroup data
    catcomp <- mapply(function(y.s,y.d) {
      (length(y.s[[1]]) == (length(y.d[[1]])))
    },cat.sub,cat.data)
    exclude <- c(exclude,which(!catcomp))
  }

  return(unique(exclude))             #return vector with items to be eliminated
}

### datprep
`datprep_LLTM` <-
  function(X,W,mpoints,Groups,sum0)
  {
    # Design matrix see Fischer & Molenaar, p. 159

    #TFrow <- (rowSums(X)==0 | rowSums(X)==(dim(X)[2]))  #el. persons with 0/K rawscore
    #X <- X[!TFrow,]

    ngroups <- max(Groups)
    X01 <- X
    N <- dim(X)[1]                                  #number of persons
    K <- dim(X)[2]/mpoints                            #number of items
    mt_vek <- rep(1,K)

    #automatized generation of the design matrix W
    if (length(W)==1) {
      W11diag <- diag(1,(sum(mt_vek)-1))                #build up design matrix
      if (sum0) {
        w110 <- rep(-1,(sum(mt_vek)-1))                 #sum0 restriction
      } else {
        w110 <- rep(0,(sum(mt_vek)-1))                  #first item category parameter set to 0
      }
      W11 <- rbind(w110,W11diag)                        #RM design matrix
      ZW <- dim(W11)[1]

      W1 <- NULL
      for (i in 1:(mpoints*ngroups)) W1 <- rbind(W1,W11)    #first part with virtual items

      if (mpoints > 1) {                                    #more than 1 measurement points
        if (ngroups > 1) {                                  #more than 1 group/more mpoints
          t_mp1 <- rep(1:mpoints,rep(ZW*ngroups,mpoints))
          t_mp <- factor(t_mp1)
          g_ng1 <- rep(rep(1:ngroups,rep(ZW,ngroups)),mpoints)
          g_ng <- factor(g_ng1)
          W2 <- model.matrix(~t_mp+g_ng)[,-1]               #main effects g and mp
          W2[1:(ZW*ngroups),] <- 0                          #remove main effects for the first test occasion
        } else {                                            #1 group/more mpoints
          t_mp <- gl(mpoints,ZW)                            #factor for measurement points
          W2 <- model.matrix(~t_mp)[,-1] }
      } else if (ngroups > 1) {                             #1 mpoint/more groups
        g_ng <- gl(ngroups,ZW)
        W2 <- model.matrix(~g_ng)[,-1]
        warning("Group contrasts without repeated measures can not be estimated!")
      } else if (ngroups == 1) W2 <- NULL                   #1 mpoint/1 group

      W <- cbind(W1,W2)
      colnames(W) <- NULL
      rownames(W) <- NULL
    }

    list(X=X,X01=X01,mt_vek=mt_vek,W=W)
    #Output: X01      ... 0/1 response matrix of dimension N*rtot
    #        mt_vek   ... vector of length K with number of categories - 1 (for each item)
    #        W        ... design matrix of dimension (K*T)*((K-1)*(T-1)+1)
  }

`datprep_PCM` <-
  function(X,W,sum0)
  {
    #... X: data matrix with response categories to be converted into 0/1 matrix

    #TFrow <- (rowSums(X)==0)  #el. persons with 0/K rawscore
    #X <- X[!TFrow,]

    #converting into 0/1 matrix
    N <- dim(X)[1]                                  #number of persons
    mt_vek <- apply(X,2,max,na.rm=TRUE)             #number of categories - 1 for each item
    mt_vek_0 <- mt_vek+1                            #number of categories for each item
    X01_0 <- matrix(rep(0,(N*sum(mt_vek_0))),nrow=N)#empty 0/1 matrix
    K <- length(mt_vek)                             #number of items
    cummt0 <- c(0,cumsum(mt_vek_0)[1:(K-1)])+1      #index vector for 0th category
    indmatp <- apply(X,1,function(xi) {xi+cummt0})  #preparing index matrix for 1 responses
    imp1 <- as.vector(indmatp)
    imp2 <- rep(1:N,rep(K,N))
    indmat <- cbind(imp2,imp1)                      #final index matrix for 1 responses
    X01_0[indmat] <- 1                              #0/1 matrix with 0th category

    NAindmat <- rbind(imp2,rep(1:K,N),c(t(X)))         #impose NA structure
    rownames(NAindmat) <- NULL
    NAind <- t(NAindmat[1:2,is.na(NAindmat[3,])])      #index matrix for NA's in X

    if (length(NAind) > 0) {
      NAindlist <- apply(NAind,1,function(x){
        co <- seq(cummt0[x[2]],cummt0[x[2]]+mt_vek[x[2]])
        NAind01 <- cbind(rep(x[1],length(co)),co)
        data.frame(NAind01,row.names=NULL)                                               #list with NA indices
      })
      indmatNA <- matrix(unlist(lapply(NAindlist, function(x) {t(as.matrix(x))})),ncol=2,byrow=TRUE)   #matrix with NA indices
      X01_0[indmatNA] <- NA
    }

    X01 <- X01_0[,-cummt0]                          #delete 0-category answers --> final 0/1 pattern matrix (dim N*sum(mt_vek))


    #automatized generation of the design matrix W
    if (length(W)==1) {
      W1 <- diag(1,(sum(mt_vek)-1))                   #build up design matrix
      if (sum0) {
        w1 <- rep(-1,(sum(mt_vek)-1))                         #sum0 restriction
      } else {
        w1 <- rep(0,(sum(mt_vek)-1))                          #first item parameter set to 0
      }
      W <- rbind(w1,W1)                               #PCM design matrix
      colnames(W) <- NULL
      rownames(W) <- NULL
    }

    list(X=X,X01=X01,mt_vek=mt_vek,W=W)
    #Output: X01      ... 0/1 response matrix of dimension N*rtot
    #        mt_vek   ... vector of length K with number of categories - 1 (for each item)
    #        W        ... design matrix of dimension sum(mt_vek)*sum(mt_vek)
  }


`datprep_RM` <-
  function(X,W,sum0)                       #prepares data matrix for Rasch model
  {
    X01 <- X                                        #X is already X(0,1)

    mt_vek <- rep(1,dim(X01)[2])                    #number of categories for each item
    K <- length(mt_vek)

    #automatized generation of the design matrix W
    if (length(W)==1) {
      W1 <- diag(1,(K-1))                           #build up design matrix
      if (sum0) {
        w1 <- rep(-1,(K-1))                         #sum0 restriction
      } else {
        w1 <- rep(0,(K-1))                          #first item parameter set to 0
      }
      W <- rbind(w1,W1)                             #RM design matrix
      colnames(W) <- NULL
      rownames(W) <- NULL
    }
    list(X=X,X01=X01,mt_vek=mt_vek,W=W)
    #Output: X01      ... 0/1 response matrix of dimension N*rtot
    #        mt_vek   ... 1-vector of length K
    #        W        ... design matrix of dimension K*K
  }


`datprep_RSM` <-
  function(X,W,sum0)
  {
    #... X: data matrix with response categories to be converted into 0/1 matrix

    max.it <- apply(X,2,max,na.rm=TRUE)             #RSM check for equal number of categories
    if (length(table(max.it)) > 1) stop("RSM can not be computed since number of categories are not the same for each item!\n")

    N <- dim(X)[1]                                  #number of persons
    K <- dim(X)[2]                                  #number of items
    hmax <- max(X,na.rm=TRUE)                       #highest category
    mt_vek <- rep(hmax,K)                           #vector with number of categories - 1 for each item

    mt_vek_0 <- mt_vek+1                            #number of categories for each item
    X01_0 <- matrix(rep(0,(N*sum(mt_vek_0))),nrow=N) #empty 0/1 matrix
    K <- length(mt_vek)
    cummt0 <- c(0,cumsum(mt_vek_0)[1:(K-1)])+1      #index vector for 0th category
    indmatp <- apply(X,1,function(xi) {xi+cummt0})  #preparing index matrix for 1 responses
    imp1 <- as.vector(indmatp)
    imp2 <- rep(1:N,rep(K,N))
    indmat <- cbind(imp2,imp1)                      #final index matrix for 1 responses
    X01_0[indmat] <- 1                              #0/1 matrix with 0th category

    NAindmat <- rbind(imp2,rep(1:K,N),c(t(X)))         #impose NA structure
    rownames(NAindmat) <- NULL
    NAind <- t(NAindmat[1:2,is.na(NAindmat[3,])])      #index matrix for NA's in X

    if (length(NAind) > 0) {
      NAindlist <- apply(NAind,1,function(x){
        co <- seq(cummt0[x[2]],cummt0[x[2]]+mt_vek[x[2]])
        NAind01 <- cbind(rep(x[1],length(co)),co)
        data.frame(NAind01,row.names=NULL)                                               #list with NA indices
      })
      indmatNA <- matrix(unlist(lapply(NAindlist, function(x) {t(as.matrix(x))})),ncol=2,byrow=TRUE)   #matrix with NA indices
      X01_0[indmatNA] <- NA
    }

    X01 <- X01_0[,-cummt0]                          #delete 0-category answers --> final 0/1 pattern matrix (dim N*sum(mt_vek))

    #automatized generation of the design matrix W
    if (length(W)==1) {
      e_it <- gl(K,hmax)                              #factor for item parameters
      e_cat <- gl(hmax,1,K*hmax)                      #factor for category par

      if (sum0) {
        Xm <- model.matrix(~e_it+e_cat)[,-1]          #dummy coding
        Xm[1:hmax,1:(K-1)] <- -1                      #first item to be sum0 normalized
      } else {
        Xm <- model.matrix(~e_it+e_cat)[,-1]          #design matrix with 0/1 contrasts (without intercept)
      }

      catvek <- 1:hmax                                #preparing the item design vectors
      e_itnew <- catvek*Xm[,1:(K-1)]
      Xm[,1:(K-1)] <- e_itnew
      W <- Xm                                         #final design matrix
      colnames(W) <- NULL
      rownames(W) <- NULL
    }

    list(X=X,X01=X01,mt_vek=mt_vek,W=W)
    #Output: X01      ... 0/1 response matrix of dimension N*rtot
    #        mt_vek   ... vector of length K with number of categories - 1 (for each item)
    #        W        ... design matrix of dimension sum(mt_vek)*((K-1)+(hmax-1))
  }


# https://github.com/cran/eRm/blob/master/R/zzz.R


setClass("prediction",
         representation(predictions = "list",
                        labels      = "list",
                        cutoffs     = "list",
                        fp          = "list",
                        tp          = "list",
                        tn          = "list",
                        fn          = "list",
                        n.pos       = "list",
                        n.neg       = "list",
                        n.pos.pred  = "list",
                        n.neg.pred  = "list"))

setClass("performance",
         representation(x.name       = "character",
                        y.name       = "character",
                        alpha.name   = "character",
                        x.values     = "list",
                        y.values     = "list",
                        alpha.values = "list" ))

#setMethod("plot",signature(x="performance",y="missing"),
#          function(x,y,...) {
#              .plot.performance(x,...)
#          })

prettyPaste <- function(...){
  paste(strwrap(paste0(..., collapse = ""), width = getOption("width")), sep="\n", collapse="\n")
}
