######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# eRm_cml
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine related to computation of
# score funtion and Hessian matrix.
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2019, Last Modified 15/08/2019
######################################################################
#' @importFrom numDeriv hessian jacobian

eRm_cml <- function(X, eta, W, model = "RM") {
  # X = observed data matrix
  # eta - item easiness parameter
  # W = design matrix
  # model = RM, PCM, RSM, LLTM

  call<-match.call()

  # Default values
  mpoints <- 1    # LPCM, LRSM, LLTM
  groupvec <- 1   # LPCM, LRSM, LLTM
  sum0 = FALSE
  se = FALSE

  if (missing(eta)) stop("eta Vector missing!")
  else eta <- as.vector(eta)

  if (missing(W)) W <- NA
  else W <- as.matrix(W)

  ###################################################################
  # cml() sub function abstracted from fitcml.R, part of eRm package
  # https://github.com/cran/eRm/blob/master/R/fitcml.R
  ###################################################################
  #cml function for call in nlm
  cml <- function(eta){

    beta <- as.vector(W %*% eta)
    #FIXME!!! gby??
    beta.list <- split(beta,gind)  #gind index for treatment groups
    beta.list1 <- beta.list

    #beta and NAstructure (over Groups): 1st line parameter values, 2nd line which item NA
    betaNA <- mapply(function(x,y) {rbind(x,y)},beta.list1,NAstruc,SIMPLIFY=FALSE)

    #likelihood term based on gamma functions for each Group x NAgroup combination
    Lg <- lapply(betaNA, function(betaNAmat) {
      beta.vec <- betaNAmat[1,]                #get parameter vector beta

      #gamma functions for each NAgroup within Groups
      Lg.NA <- apply(matrix(betaNAmat[-1,],ncol=length(beta.vec)),1, function(NAvec) {

        #list of virtual item-category parameters per item
        beta_list <- as.list(split(beta.vec[NAvec==1],mt_ind[1:(length(beta.vec[NAvec==1]))]))
        parlist <- lapply(beta_list,exp)                                #initial epsilon as list

        #------------------gamma functions----------------------
        g_iter <- NULL                                                  #computation of the gamma functions
        K <- length(parlist)
        for (t in 1:(K-1)) {                                            #building up J1,...,Jt,...,Js

          if (t==1) {                                                   #first iteration step
            gterm <- c(1,parlist[[t]])                                  #0th element included
          }else
          {
            gterm <- g_iter                                   #gamma previous iteration with 0th el
            g_iter <- NULL
          }

          parvek <- c(1,parlist[[t+1]])                      #eps vector in current iteration with 0th el
          h <- length(parvek)                                #dimensions for matrix
          mt <- length(gterm)
          rtot1 <- h+mt-1                                    #number of possible raw scores (0 included)

          gtermvek <- rep(c(gterm,rep(0,h)),h)                          #building up matrix for gamma term
          gtermvek <- gtermvek[-((length(gtermvek)-h+1):length(gtermvek))]      #eliminating last h 0's
          gmat <- matrix(gtermvek,nrow=rtot1,ncol=h)
          emat <- matrix(rep(parvek,rep(rtot1,h)),ncol=h,nrow=rtot1)    #building up matrix for eps term
          gmat_new <- gmat*emat                                                 #merge matrices
          g_iter <- rowSums(gmat_new)                     #gamma functions in current iteration are rowsums
        }
        #----------------- end gamma functions ------------------

        Lg.NA <- as.vector(g_iter[2:(rtot+1)])     #final gamma vector stored in gamma (without gamma0)
        return(Lg.NA)
      })
    })
    #----------------- compute likelihood components -----------------------
    L1 <- sum(mapply(function(x,z) {
      x[!is.na(z)]%*%na.exclude(z)
    },nrlist,lapply(Lg,log)))        #sum up L1-terms (group-wise)

    L2 <- sum(mapply("%*%",x_mtlist,beta.list1))        #sum up L2-terms (group-wise)
    L1-L2                                               #final likelihood value
  }
  #----------------- end likelihood -----------------------

  #########################################################################
  # R code abstracted and adapted from likeLR.R, part of eRm package
  # URL: https://github.com/cran/eRm/blob/master/R/likLR.R
  #########################################################################
  Groups <- groupvec  # function(,... Groups,...)

  if ( any( is.na(X) )) {
    dichX <- ifelse( is.na(X), 1,0)
    strdata <- apply( dichX, 1, function(x) {paste(x,collapse="")})
    gmemb <- as.vector( data.matrix( data.frame( strdata ) ))
  } else {
    gmemb <- rep( 1, dim(X)[1] )
  } # end if

  #data preparation, design matrix generation for various models
  if (model=="RM") { Xprep <- datprep_RM( X, W, sum0)
  } else if (model=="LLTM") {Xprep <- datprep_LLTM( X, W, mpoints, Groups, sum0)
  } else if (model=="RSM") {Xprep <- datprep_RSM( X, W, sum0)
  } else if (model=="PCM") {Xprep <- datprep_PCM( X, W, sum0)
#  } else if (model=="LRSM") {Xprep <- datprep_LRSM( X, W, mpoints, Groups, sum0)
#  } else if (model=="LPCM") {Xprep <- datprep_LPCM(X, W, mpoints, Groups, sum0)
  }

  Lprep <- cmlprep( X01 =    Xprep$X01,
                          mt_vek = Xprep$mt_vek,
                          mpoints = mpoints,
                          Groups = groupvec,
                          W = Xprep$W,
                          gmemb = gmemb)

  #################################################################
  # likeLR.R
  #################################################################

  # arguments used in cml() function stored in eRm_cml environment
  mt_ind = Lprep$mt_ind
  nrlist = Lprep$nrlist
  x_mt =   Lprep$x_mt
  rtot =   Lprep$rtot
  W = Xprep$W
  ngroups = max(Groups)
  gind = Lprep$gind
  x_mtlist = Lprep$x_mtlist
  NAstruc = Lprep$NAstruc
  g_NA = Lprep$g_NA
  st.err = se
  gby = Lprep$gby

  # etapar <- -etapar          # RM() output difficulty for RM, RSM, PCM --> input = easiness! AK
  if (model=="RM" || model=="RSM" || model=="PCM") { eta <- -eta}

  loglik <- - cml(eta)  # conditional log likelihood evaluated at eta

#  e <- globalenv()   # global environment
#  final.gamma <- e$Lg[[1]]   # get final gamma vector from global environment

  ######################################################################
  # calculate score function and Hessian matrix using numDeriv package
  #####################################################################
  scorefun <- as.vector(numDeriv::jacobian(cml, eta))    # Jacobian = scorefun as vector
  r.hessian <- numDeriv::hessian(cml, eta) # hessian matrix

 # rm("Lg", envir = e)  # remove Lg from global environment

  par.list <-      list( "eta" = eta,
                         "mt_ind" = mt_ind,
                         "nrlist" = nrlist,
                         "x_mt" = x_mt,
                         "rtot" = rtot,
                         "W" = W,
                         "ngroups" = max(Groups),
                         "gind" = gind,
                         "x_mtlist" = x_mtlist,
                         "NAstruc" = NAstruc,
                         "g_NA" = g_NA,
                         "st.err" = se,
                         "gby" = gby,
                         "X01" = Xprep$X01,
                         "call" = call)

  return( list( "loglik" = loglik,
            #    "gamma"= final.gamma,
                "scorefun" = scorefun,
                "hessian" = r.hessian,
                "par.list" = par.list ) )
}
