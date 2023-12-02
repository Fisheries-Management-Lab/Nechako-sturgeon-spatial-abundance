# SMR_TMB.R
# Brett van Poorten
# April 2, 2023
# simulates data for the spatial mark recapture model on Nechako sturgeon
# to allow simulation-evaluation; also calls real data to do estimation

require(TMB)
library(tmbstan)
require(ggplot2)
require(reshape2)
require(tidyr)
require(dplyr)
library(cowplot)
library(grid)
library(gridExtra)
library(sf)
library(raster)
library(spData)
library(spDataLarge)
library(tmap)
library(leaflet)
library(bcmaps)

###########################################################################
########################### available functions ###########################
###########################################################################
#
# Simulator() - workhorse simulation function
# rand_par()  - randomizes parameters
# invgrav() - take proportional distribution and calculate Markov and gravity parameters
# IBM() - individually move fish among areas
# tagprocessing() - individual tagging history (marks and subsequent recaptures)
# catprocessing() - first captures put in a matrix (regardless of where you find them again)
# nocapprocessing() - summarize tagging data into numbers moving between areas
# allcapprocessing() - individual tagging history of every capture (including last observation)
# buildrealdata() - transform loaded datasets into appropriate objects for TMB import
# optimize() - simulation-estimation comparisons using different movement models
# fit_data() - fit the model to data (frequentist and Bayesian)
# PostPred() - take posterior estimates and recreate the population
# TransMat() - take the fitted model and provide estimated parameters and Pdist
### PLOTS
# bias() - do simulation-evaluation and THEN plot it
# plotTrend() - plot population trajectory and estimated parameter output across space models
# compareM() - look at posterior updates across different M priors
# Figure1Map() - show study area

sim     <- TRUE  # if TRUE, runs simulation-evaluation
est     <- FALSE # if TRUE, estimates parameters from real data
path    <- getwd()

###########################################################################
################################# indices #################################
###########################################################################

ny      <- 21    # number of years
na      <- 4     # number of areas
ntel    <- 1496   # number of telemetry observations

###########################################################################
################################ parameters ###############################
###########################################################################

meanM   <- 0.04           # mean instantaneous mortality
N1      <- 1000           # initial abundance
TL      <- 0 #0.00204        # instantaneous tag loss rate
Pdist   <- rep(1/na,na)   # proportion in each area
gravsd  <- 1              # standard deviation in gravity across areas
visc    <- 2              # weight on maintaining residency
fullmov <- TRUE
muP     <- 0.05           # mean of capture probability
cvP     <- 0.1            # cv of capture probability
Ptarg   <- 1#0.4          # target P ( ?hyperstability parameter? )
DR      <- 1            # capture rate adjustment for years after first

###########################################################################
################################### data ##################################
###########################################################################

real_tagdat <- read.csv(file=paste0(path,"/Data/Nechako tagging.csv"))
real_tel    <- read.csv(file=paste0(path,"/Data/Nechako telemetry.csv"))

###########################################################################
################################ functions ################################
###########################################################################

"Simulator" <- function(ny,na,M,N1,TL,fullmov,Pdist,gravsd,visc,ntel,muP,cvP,
                        Ptarg,DR){
  mean_pars <- list(
    M=M, 
    N1=N1, 
    TL=TL,
    fullmov=fullmov, 
    Pdist=Pdist, 
    gravsd=gravsd, 
    visc=visc
  )
  ipars   <- rand_par(na, mean_pars)
  N1  <- unlist( ipars$N1 )
  TL  <- unlist( ipars$TL )
  M   <- unlist( ipars$M )
  mov <- unlist( ipars$mov )
  G   <- unlist( ipars$G )
  
  # N is total number of fish
  # TT is number of tag captures
  N  <- TT <- array( dim=c( N1, ny, na ))

  # random detection probability in an area and year for binomial sampling
  PP      <- array( data=rlnorm(ny*na, log(muP), cvP), dim=c(ny,na))   
  Ct      <- array( data=NA, dim=c(ny,na))
  
  # set equilibrium distribution (called d in paper)
  initdist <- apply( mov,2,sum )/na
  for( j in 1:400 )
    initdist <- apply( initdist*mov, 2, sum )
  
  # randomly distribute initial population based on proportional distribution
  N[,1,]    <- t( rmultinom( N1, 1, initdist ))
  # calculate hyperstable sampling probability 
  totalfrac <- sum( initdist*PP[1,])
  newfrac   <- sum( initdist*PP[1,]^Ptarg )
  Prat      <- totalfrac/newfrac
  
  # TRUE/FALSE if a fish was in an area and captured
  SAMP      <- array(runif(na*N1)<rep(Prat*PP[1,]^Ptarg,each=N1),c(N1,na))  
  # calculate number of fish captured in first time period
  TT[,1,]   <- N[,1,] * SAMP
  
  # simulate for remaining time steps
  for( y in 2:ny){
    # mortality
    N[,y,]     <- N[,y-1,] * rbinom( n=N1, size=1, prob=exp(-M) ) 
    # redistribute among areas
    N[,y,]     <- IBM( N[,y,], mov )                              
    # distribution among areas
    Ndist      <- apply(N[,y,],2,sum) / sum( N[,y,] )             
    # proportion potentially sampled
    totalfrac  <- sum( Ndist*PP[y,] )                             
    # corrected with Ptarg to account for hyperstability
    newfrac    <- sum( Ndist*PP[1,]^Ptarg )                       
    # proportional sampling rate
    Prat       <- totalfrac/newfrac                               
    # binomial sampling
    SAMP       <- array( runif(na*N1)<rep(Prat*PP[y,]^Ptarg,each=N1), c(N1,na))
    # binomial tag loss
    TagLoss    <- array( runif(na*N1)<rep(exp(-TL), each=N1*na), c(N1,na))
    # reduce sampling probability using DR in a second binomial sampling
    TT[,y,]    <- N[,y,]*SAMP*(runif(N1)<DR)#*TagLoss
  }
  
  # summarize data
  # Create capture history (if needed)
  CH <- apply( TT, 1:2, sum )
  CH <- CH[apply(CH, 1, sum) > 0]
  
  # individual tagging history (for fish that were recaptured)
  tagdat <- allcapprocessing( TT )
  
  # year, area, and number of last captures
#  nocap  <- nocapprocessing( TT ) - this is now in allcapprocessing()
  
  # year, area, and number of first captures
  Ct     <- catprocessing( TT )
  
  # total number of captures by year and area
  Ctot   <- apply( TT, 2:3, sum )
  
  # telemetry transitions - all areas over all years; equal area distribution
  # This just distributes the ntel observations across area combinations (a to k)
  teldat <- array( rmultinom( 1, ntel, mov ), dim=c(na,na) )
  
  # put data objects in the global environment

  assign("tagdat",tagdat,envir=.GlobalEnv)  # tag-recapture data
  assign("Ct",Ct,envir=.GlobalEnv)          # initial captures
  assign("Ctot",Ctot,envir=.GlobalEnv)      # total captures per area/time
  assign("teldat",teldat,envir=.GlobalEnv)  # telemetry movements
#  assign("nocap",nocap,envir=.GlobalEnv)
  assign("N",N,envir=.GlobalEnv)
  assign("TT",TT,envir=.GlobalEnv)
  assign("mov",mov,envir=.GlobalEnv)
  assign("PP",PP,envir=.GlobalEnv)
  assign("idist",initdist,envir=.GlobalEnv)

  out <- list(
    M     = M,
    N1    = N1,
    TL    = TL,
    mov   = mov,
    N     = N,
    G     = G,
    idist = initdist
  )
  return( out )
} # simulator()

"rand_par" <- function( na, pars ){  # randomize parameters
  M       <- unlist( pars[1] )
  N1      <- unlist( pars[2] ) 
  TL      <- unlist( pars[3] )
  fullmov <- unlist( pars[4] )
  Pdist   <- unlist( pars[5] )
  gravsd  <- unlist( pars[6] )
  visc    <- unlist( pars[7] )
  
  iN1    <- rpois( 1, N1 )
  iTL    <- rlnorm( 1, log(TL), 0.1 )
  iM     <- rlnorm( 1, log(M), 0.1 )
  visc   <- rlnorm(1, log(visc), 0.2 )
  # define transition matrix
  if( fullmov ){
    # full Markov matrix
    # take Pdist, calculate G matrix and add variation - 
    gravs <- c( rep(0,na),rep(invgrav(Pdist)[2:na],each=na)+
                  rnorm((na-1)*na,0,gravsd))  
    mov <- array(gravs,dim=c(na,na)) # put into G matrix
    ind<-cbind(1:na,1:na)        # create indices relating to a==k
    mov[ind]<-mov[ind]+visc      # add on residency
    mov<-exp(mov)                # exponentiate
    mov<-mov/rep(apply(mov,1,sum),na) # calculate transition matrix
    movest <- mov
    
    # now calculate initial G matrix
    den <- vector()
    gmat <- array(NA,dim=c(na,na-1))
    for(i in 1:na){                # move down columns
      den[i] <- 1/movest[i,1]
      for( j in 1:(na-1) ){            # move across rows
        gmat[i,j] <- log(movest[i,j+1]*den[i])
      }
    }
    G <- as.vector(gmat)
    
  } else {
    # gravity matrix
    # take Pdist, calculate vector and add variation 
    gravs <- c(0, invgrav(Pdist)[2:na]+rnorm(na-1,0,sd=gravsd))                        
    mov <- array(rep(gravs,each=na),dim=c(na,na)) # put into G matrix
    ind <- cbind(1:na,1:na)      # create indices relating to a==k
    mov[ind] <- mov[ind] + visc  # add residency
    mov <- exp(mov)              # exponentiate
    mov <- mov/rep(apply(mov,1,sum),na) # calculate transition matrix
    movest <- mov

    # now calculate g vector and v 
    G <- vector()
    den <- vector()
    expG <- vector()
    for(i in 2:na){
      den[i] <- 1/movest[i,1]
      if(i < na){
        G[i] <- log(movest[i,i+1]*den[i])
      } else {
        G[1] <- log(movest[i,2]*den[i])
      }
    }
    expG <- exp( c( 0, G ) )
    G[na] <- log( den[na] - sum(expG[1:(na-1)]) ) - G[na-1]   #viscosity parameter
    G[na] <- log(G[na])  # in TMB, viscosity parameter is read in log space

  }
  
  out   <- list(
    N1 = iN1,
    TL = iTL,
    M = iM,
    mov = mov,
    G = G
  )
  
  return( out )
} # rand_par()

# take Markov gravity from a row of the transition matrix and calculate g parameters
"invgrav" <- function(x){  
  g              <- rep( 0, length(x) )
  g[2:length(x)] <- log(x[2:length(x)]/x[1])
  return(g)
} # invgrav()

# redistribute individuals among areas
"IBM" <- function( nn, mov ){ 
  # nn is number of fish across areas
  # mov is transition matrix
  na          <- nrow(mov)    # number of areas
  N1          <- dim(nn)[1]   # number of fish (alive or not)
  loc.matrix  <- array( rep( 1:na, each=N1 ), dim=c(N1,na)) # locations
  aa          <- apply( loc.matrix*nn, 1, sum ) # what area each fish is currently in
  nn2         <- array( 0, dim(nn))
  stillalive  <- which(aa > 0 )  # which haven't died already
  for( j in stillalive )
    nn2[j,]    <- rmultinom(1,1,mov[aa[j],]) # randomly redistribute
  return(nn2)
} # IBM()

# individual tagging history (marks and subsequent recaptures)
"tagprocessing" <- function( tt ){ 
  N1      <- dim( tt )[1]
  ny      <- dim( tt )[2]
  na      <- dim( tt )[3]
  
  tags    <- NULL
  yind    <- rep( 1:ny, each=na )
  aind    <- rep( 1:na, ny )
  
  for( i in 1:N1 ){
    yrs   <- yind[ as.vector( t(tt[i,,]))==1 ] # years when fish was seen
    areas <- aind[as.vector( t(tt[i,,]))==1 ]  # areas where fish was seen
    nobs  <- length(yrs)
    if( nobs>1 ){
      for( k in 1:(nobs-1) ){
        # add a line with area from, year1; area to, year2
        tags <- rbind(tags, c(areas[k],yrs[k],yrs[k+1],areas[k+1]))  
      } # i
    } # if(nobs>1)
  } # for i
  
  # summarize tagging data to provide numbers moving between all combination of areas
  tagdat <- aggregate(rep(1,nrow(tags)), 
                      by=list(tags[,1],tags[,2],tags[,3],tags[,4]),sum ) 
  names(tagdat)   <- c("Arel","Yrel","Ycap","Acap","N")
  return( tagdat )
} # tagprocessing()

# first captures put in a matrix (regardless of where you find them again)
"catprocessing" <- function( tt ){  
  N1      <- dim( tt )[1]
  ny      <- dim( tt )[2]
  na      <- dim( tt )[3]
  
  cats    <- NULL
  caught  <- apply( tt, 1, sum )  # number of times each fish was captured
  yind    <- rep( 1:ny, each=na )
  aind    <- rep( 1:na, ny )
  
  for( k in (1:N1)[caught>0] ){
    # times and areas when this fish was first caught
    ind   <- match(1, as.vector( t( tt[k,,] ))) 
    cats  <- rbind(cats, c( yind[ind], aind[ind] ))
  } # k
  
  # summarize tagging data to provide times and areas when each fish was captured 
  catdat_r  <- aggregate(rep(1,nrow(cats)), by=list( cats[,1], cats[,2]), sum )
  names( catdat_r ) <- c("Year","Area","Catch")
  catdat    <- array( 0, dim=c(ny,na) )
  catdat[as.matrix(catdat_r[,1:2])] <- catdat_r[,3]
  return( catdat )
} # catprocessing()

# summarize tagging data to provide numbers moving between all combination of areas
"nocapprocessing" <- function( tt ){
  N1      <- dim( tt )[1]
  ny      <- dim( tt )[2]
  na      <- dim( tt )[3]
  
  nocaps    <- NULL
  caught  <- apply( tt, 1, sum )  # number of times each fish was captured
  yind    <- rep( 1:ny, each=na )
  aind    <- rep( 1:na, ny )
  
  for( k in (1:N1)[caught > 0]){
    TN      <- t( tt[k, ny:1,])     # captures by area and time (indexed backwards)
    ind     <- match( 1, as.vector(TN) ) # index of times and areas last captured
    nocaps  <- rbind( nocaps, c(yind[ind],aind[ind]) )
  } # k
  
  nocap_r  <- aggregate( rep(1,nrow(nocaps)), 
                         by=list( nocaps[,1],nocaps[,2]), sum )
  names( nocap_r ) <- c("Year","Area","Catch")
  nocap    <- array( 0, dim=c(ny,na) )
  nocap[as.matrix(nocap_r[,1:2])] <- nocap_r[,3]
  return( nocap )
} # nocapprocessing()

# individual tagging history of every capture (including last observation)
"allcapprocessing" <- function( tt ){ 
  N1      <- dim( tt )[1]
  ny      <- dim( tt )[2]
  na      <- dim( tt )[3]
  
  tags    <- NULL
  yind    <- rep( 1:ny, each=na )
  aind    <- rep( 1:na, ny )
  
  for( i in 1:N1 ){
    yrs   <- yind[ as.vector( t(tt[i,,]))==1 ] # years when fish was seen
    areas <- aind[as.vector( t(tt[i,,]))==1 ]  # areas where fish was seen
    id    <- (yrs-1)*na+areas
    nobs  <- length(yrs)
    if( nobs>1 ){
      for( k in 1:(nobs-1) ){
        # add a line with area from, year1; area to, year2
        tags <- rbind(tags, c(id[k],areas[k],yrs[k],id[k+1],yrs[k+1],areas[k+1]))  
      } # i
    }
    if( nobs >= 1 ){
      # add on last observation
      tags <- rbind(tags, c(id[nobs],areas[nobs],yrs[nobs],na*ny+1,ny+1,na+1) ) 
    } # if(nobs>=1)
  } # for i
  
  # summarize tagging data to provide numbers moving between all combination of areas
  tagdat <- aggregate(rep(1,nrow(tags)), 
                      by=list(tags[,1],tags[,2],tags[,3],
                              tags[,4],tags[,5],tags[,6]),sum ) 
  names(tagdat)   <- c("idrel","Arel","Yrel","idcap","Ycap","Acap","N")
  tagdat <- tagdat[ with( tagdat, order( idrel, idcap )), ]
  return( tagdat )
} # allcapprocessing()

"buildrealdata" <- function(maxyr=2023){
  
  # take tagging data, include years>2000; find only first observation in a 
  #   year, sort by ID and year
  rtagdat        <- real_tagdat %>% 
    filter( Year > 2000 & Year < maxyr ) %>%  # only obs after 2000
    mutate( ID = match( PIT, unique(PIT ))) %>% # create ID based on tag #
    mutate( YR = Year-min(Year)+1 ) %>% # create YR index
    distinct( ID, YR, .keep_all = TRUE ) %>% # keep only first observation each year
    arrange( ID, YR, AREA ) # sort by ID, Year, Area
  
  # find indices of area, year and tags
  na             <- max( unique( rtagdat$AREA ))
  ny             <- max( unique( rtagdat$YR ))
  ntag           <- length( unique( rtagdat$ID ))

  # populate tagging capture history across areas
  tt             <- array( data=0, dim=c(ntag,ny,na) )
  indices        <- with( rtagdat, cbind( ID, YR, AREA ))
  
  tt[indices]    <- 1

  # individual tagging history (for fish that were recaptured)
  tagdat <- allcapprocessing( tt )
  
  # year, area, and number of first captures
  Ct     <- catprocessing( tt )
  
  # total number of captures by year and area
  Ctot   <- apply( tt, 2:3, sum )
  
  # take telemetry dataset and create transition matrix
  rteldat         <- real_tel %>%
    filter( Year > 2000 & Year < maxyr ) %>%  # keep only obs after 2000
    mutate( ID = match( PIT, unique(PIT ))) %>% # create ID based on pit code
    mutate( YR = Year-min(Year)+1 ) %>% # create year index
    distinct( ID, YR, .keep_all = TRUE ) %>% # keep only first obs each year
    arrange( ID, YR, AREA ) # sort by ID, Year, Area
  
  tagID   <- sort(unique(rteldat$PIT)) # unique tag codes
  tel     <- array( data=0, dim=c(na,na))  # array to fill
  for( i in 1:length(tagID) ){
    inddat    <- rteldat %>% filter( PIT==tagID[i] ) # isolate obs by fish
    nobs      <- dim( inddat )[1] # number of obs for that fish
    if( nobs > 1 ){
      for( j in 2:dim(inddat)[1] ){
        tel[inddat$AREA[j-1],inddat$AREA[j]] <- 
          tel[inddat$AREA[j-1],inddat$AREA[j]] + 1
      } # j
    } # if(nobs)
  } # i
  
  # put data objects in the global environment
  
  assign("tagdat",tagdat,envir=.GlobalEnv)  # tag-recapture data
  assign("Ct",Ct,envir=.GlobalEnv)          # initial captures
  assign("Ctot",Ctot,envir=.GlobalEnv)      # total captures per area/time
  assign("teldat",tel,envir=.GlobalEnv)     # telemetry movements
  cat("Data loaded and put into appropriate objects\n")
}

###########################################################################
############################### main section ##############################
###########################################################################

# simulation-estimation comparisons using different movement models
"optimize" <- function(nrep=1, simMove=c("Markov","Grav"), estMove=c("Markov","Grav")){
  sim_N <- vector()
  sim_decline <- vector()
  est_N <- vector()
  est_decline <- vector()
  # TRUE/FALSE flag indicating if simulation/estimation is Markov or Gravity
  Sim_MovMod <- simMove == "Markov"   
  Est_MovMod <- estMove == "Markov"

  for( r in 1:nrep ){
    # generate random data
    sim_dat <- Simulator(ny,na,meanM,N1,TL,fullmov=Sim_MovMod,
                         Pdist,gravsd,visc,ntel,muP,cvP,Ptarg,DR)
    
    movest <- sim_dat$mov
    G <- vector()
    if(estMove == "Grav"){  # calculate G parameters for gravity model
      den <- vector()
      expG <- vector()
      for(i in 2:na){
        den[i] <- 1/movest[i,1]
        if(i < na){
          G[i] <- log(movest[i,i+1]*den[i])
        } else {
          G[1] <- log(movest[i,2]*den[i])
        }
      }
      expG <- exp( c( 0, G ) )
      G[na] <- log( den[na] - sum(expG[1:(na-1)]) ) - G[na-1]   #viscosity parameter
      G[na] <- log(G[na])  # in TMB, viscosity parameter is read in log space
    } else {                  # calculate G parameters for Markov model
      den <- vector()
      gmat <- array(NA,dim=c(na,na-1))
      for(i in 1:na){                # move down columns
        den[i] <- 1/movest[i,1]
        for( j in 1:(na-1) ){            # move across rows
          gmat[i,j] <- log(movest[i,j+1]*den[i])
        }
      }
      G <- as.vector(gmat)
    }
    if( is.na(sum(G)) )
      next
    
    # get TMB ready
    years  <-array( 1:ny, dim=c(ny,na))
    areas  <-array(rep(1:na,each=ny),dim=c(ny,na))
    Find   <- cbind(years[Ctot>0],areas[Ctot>0])
    Fmat   <- -log( 1- Ctot/(apply(sim_dat$N,c(2,3),sum)+1) )
    Fvec   <- Fmat[Find]
    
    ## Change data names to capture and recapture 
    DATA <- list(
      est       = 0,
      fullmov   = sum(Est_MovMod),
      CV        = 0.05,
      Mmu       = log(meanM),
      Msd       = 0.4,
      tagind    = as.matrix(tagdat[,1:6])-1, # tag-recpature indices
      tagdat    = as.vector(tagdat[,7]),     # tag-recaptures
#      nocap     = nocap,
      cat       = Ct,                        # initial captures
      cattot    = Ctot,                      # total captures per area/time
      tel       = teldat,
      Find      = Find-1                     # F's when captures happened
    )
  
    PARAMETERS <- list(
      lnN0      = log( sim_dat$N1 ),
      lnM       = log( sim_dat$M ),
      lnF       = log( Fvec ),
      movest    = G
    )
    
    MAP        <- list(
      #lnN0     = factor(NA),
      #lnM      = factor(NA),
      lnF      = factor(rep(NA,dim(Find)[1]))
    )
  
    compile("SMRtmb.cpp")
    dyn.load(dynlib("SMRtmb"))
    obj <- MakeADFun(
                      data = DATA, 
                      parameters = PARAMETERS,
                      hessian = TRUE,
                      #map = MAP,
                      DLL = "SMRtmb"
                    )
    opt <- nlminb(objective = obj$fn, 
               start = obj$par, 
               gradient = obj$gr, 
               control=list(eval.max=1000,
                            iter.max=1000),
               lower=c(-20,-20,rep(-10,dim(Find)[1]),rep(-6,na*na)),
               upper=c(20,20,rep(0.4,dim(Find)[1]),rep(6,na*na))
            )
    
    sim_N[r]         <- sum( sim_dat$N[,ny,])
    sim_decline[r]   <- sum( sim_dat$N[,ny,]) / sum( sim_dat$N[,1,] )
    est_N[r]         <- sum( obj$report()$N[ny,] )
    est_decline[r]   <- sum( obj$report()$N[ny,] ) / sum( obj$report()$N[1,] ) 
  } # r
  
  out <- list(
    sim_dat     = sim_dat,
    DATA        = DATA,
    PARAMETERS  = PARAMETERS,
    MAP         = MAP,
    obj         = obj,
    opt         = opt, 
    sim_N       = sim_N,
    est_N       = est_N,
    sim_decline = sim_decline,
    est_decline = est_decline
  )
  
  return( out )
}

# fit the model to data (frequentist and Bayesian)
"fit_data" <- function(maxyr = 2023, 
                       bayes = TRUE, 
                       meanM = 0.04, sdM = 0.4,
                       estMove=c("Markov","Grav")){
  buildrealdata(maxyr=maxyr)
  
  ### NOTE: NEED TO BUILD TAG LOSS INTO ESTIMATION

  # TRUE/FALSE flag indicating if simulation/estimation is Markov or Gravity
  Est_MovMod <- estMove == "Markov"
  # indices
  ny     <- dim(Ctot)[1]
  na     <- dim(Ctot)[2]

  # get TMB ready
  years  <-array( 1:ny, dim=c(ny,na))
  areas  <-array(rep(1:na,each=ny),dim=c(ny,na))
  Find   <- cbind(years[Ctot>0],areas[Ctot>0])
  N0     <- 500
  M      <- meanM
  Fvec   <- rep( 0.1, dim(Find)[1] )
  if( Est_MovMod ){
    G <- rnorm(na*(na-1))
  } else {
    G <- rnorm(na)
  } # if(Est_MovMod)
  

  ## Change data names to capture and recapture 
  DATA <- list(
    est       = 1,
    fullmov   = sum(Est_MovMod),
    CV        = 0.05,
    Mmu       = log(meanM),
    Msd       = sdM,
    tagind    = as.matrix(tagdat[,1:6])-1, # tag-recpature indices
    tagdat    = as.vector(tagdat[,7]),     # tag-recaptures
    #      nocap     = nocap,
    cat       = Ct,                        # initial captures
    cattot    = Ctot,                      # total captures per area/time
    tel       = teldat,
    Find      = Find-1                     # F's when captures happened
  )
  
  PARAMETERS <- list(
    lnN0      = log( N0 ),
    lnM       = log( M ),
    lnF       = log( Fvec ),
    movest    = G
  )
  
  MAP        <- list(
    #lnN0     = factor(NA),
    #lnM      = factor(NA),
    lnF      = factor(rep(NA,dim(Find)[1]))
  )
  
  compile("SMRtmb.cpp")
  dyn.load(dynlib("SMRtmb"))
  obj <- MakeADFun(
    data = DATA, 
    parameters = PARAMETERS,
    hessian = TRUE,
    #map = MAP,
    DLL = "SMRtmb"
  )
  opt <- nlminb(objective = obj$fn, 
                start = obj$par, 
                gradient = obj$gr, 
                control=list(eval.max=1000,
                             iter.max=1000),
                lower=c(-20,-20,rep(-10,dim(Find)[1]),rep(-6,na*na)),
                upper=c(20,20,rep(0.4,dim(Find)[1]),rep(6,na*na))
  )
  
  if( bayes ){
    cores <- 3
    parnames <- names(opt$par)
    plnN0  <- opt$par[which(parnames=="lnN0")]
    plnM   <- opt$par[which(parnames=="lnM")]
    plnF   <- opt$par[which(parnames=="lnF")]
    movest <- opt$par[which(parnames=="movest")]
    options(mc.cores = cores)
    init.fn <- function(){
      list(
        lnN0      = rnorm(1,plnN0,0.1), 
        lnM       = rnorm(1,plnM,0.01), 
        lnF       = rnorm(length(Fvec),plnF,0.001), 
        movest    = rnorm(length(G),movest,0.1)
      )
    } # init.fn()
    fit <- tmbstan(obj, chains=cores, open_progress=FALSE, init=init.fn ,
                   iter=10000, control=list(max_treedepth=12) )
    #plot( fit )
    #bayesplot::mcmc_trace( fit )
    #bayesplot::mcmc_pairs( fit, pch=19 )
    
    save( fit, file=paste0("Posterior ",estMove,meanM,sdM,".RData"))
  } # if(bayes)
  
  out <- list(
    DATA       = DATA,
    PARAMETERS = PARAMETERS, 
    obj        = obj,
    opt        = opt
  )
  
  return(out)
} # fit_data()

# take posterior estimates and recreate the population
"PostPred" <- function(estMove=c("Markov","Grav")){
  load( paste0( "Posterior ",estMove,"0.040.4.RData" ))
  buildrealdata(maxyr=2023)

  # indices
  ny     <- dim(Ctot)[1]
  na     <- dim(Ctot)[2]
  
  years  <- array( 1:ny, dim=c(ny,na))
  areas  <- array(rep(1:na,each=ny),dim=c(ny,na))
  Find   <- cbind(years[Ctot>0],areas[Ctot>0])
  nMP    <- ifelse( estMove=="Markov", na*(na-1), na )
  nF     <- dim(Find)[1]
  
  post   <- as.data.frame( fit )
  npost  <- dim(post)[1]
  N0     <- exp( post$lnN0 )
  M      <- exp( post$lnM )
  Fc     <- array( dim=c(npost,nF) )
  for( i in 1:nF ){
    Fc[,i] <- exp( post[,2+i])
  }
  movest <- array( dim=c(npost,nMP) )
  for( i in 1:nMP ){
    movest[,i] <- post[,2+nF+i]
  }
  
  # gravity matrix
  mov    <- array( data=0, dim=c(npost,na,na) )
  if( estMove == "Markov" ){
    # fully described Markov model
    ii=1
    for( a2 in 2:na ){
      for( a in 1:na ){
        mov[,a,a2] <- movest[,ii]       # probability of moving to cell is proportional to its estimated G
        ii <- ii + 1
      } # a
    } # a2
  } else {
    # Gravity model
    for( a in 1:na ){
      mov[,a,2:na] <- movest[,1:(na-1)]
      mov[,a,a] <- mov[,a,a] + exp(movest[,na]) # viscosity to staying in the same area
    } # a
  } # if(Markov)
  mov     <- exp( mov )
  mov_den <- apply( mov, c(1,2), sum )
  for( i in 1:npost ){
    mov[i,,]  <- mov[i,,] / mov_den[i,]
  }
  
  initdist <- apply( mov,c(1,3),sum )/na  # initialize
  twoDprod <- function(i){                # dotproduct function for each posterior
    return( initdist[i,] %*% mov[i,,] )
  }
  for( j in 1:400 ){
    initdist <- t(mapply( twoDprod, 1:npost ))
  }
  
  N        <- array( dim=c(npost,ny,na) )
  N[,1,]   <- N0 * initdist
  nextyr   <- function( i, Ntmp ){
    return( ( Ntmp[i,] %*% mov[i,,] )*exp(-M[i]) )
  }
  for( y in 2:ny ){
    N[,y,] <- t(mapply( nextyr, 1:npost, MoreArgs = list(Ntmp=N[,y-1,]) ))
  }
  
  out          <- list()
  out$N        <- N
  out$M        <- M
  out$mov      <- mov
  out$initdist <- initdist
  
  return(out)
} # PostPred()

"TransMat" <- function(estMove=c("Grav","Markov")){
  est <- fit_data(bayes=FALSE, estMove=estMove)
  movest <- est$opt$par[which(names(est$opt$par)=='movest')]
  
  # gravity matrix
  mov    <- array( data=0, dim=c(na,na) )
  if( estMove == "Markov" ){
    # fully described Markov model
    ii=1
    for( a2 in 2:na ){
      for( a in 1:na ){
        mov[a,a2] <- movest[ii]       # probability of moving to cell is proportional to its estimated G
        ii <- ii + 1
      } # a
    } # a2
  } else {
    # Gravity model
    for( a in 1:na ){
      mov[a,2:na] <- movest[1:(na-1)]
      mov[a,a] <- mov[a,a] + exp(movest[na]) # viscosity to staying in the same area
    } # a
  } # if(Markov)
  mov     <- exp( mov )
  mov_den <- rowSums(mov)
  mov     <- mov/mov_den

  initdist <- apply( mov,2,sum )/na
  for( j in 1:400 )
    initdist <- apply( initdist*mov, 2, sum )
  
  out <- list(
    mov=mov,
    initdist=initdist
  )
  
  return(out)
  
} # TransMat

###########################################################################
############################### plot section ##############################
###########################################################################

"bias" <- function(load=TRUE, tagloss=0.00204, tenyear=TRUE,nrep=500){
  
  bias_N <- bias_dec <- matrix(nrow=4,ncol=6) # 90th, 50th percentile
  sim_MovMod        <- c( "Grav", "Grav", "Markov", "Markov" )
  est_MovMod        <- c( "Grav", "Markov", "Grav", "Markov" )
  trial             <- list()
  N_bias_comp       <- list()
  decline_bias_comp <- list()
  TL                <<- tagloss
  if(tenyear){
    ny <<- 10
    ntel <<- 300
    if(load){
      load(file=paste0(path,"/SimEval10years.RData"))
    } else {
      for( i in 1:4 ){
        trial[[i]] <- optimize( nrep, simMove=sim_MovMod[i], estMove=est_MovMod[i] )
        if( !is.finite( sum( trial[[i]]$obj$he() ) ) ) next
        N_bias_comp[[i]] <- (trial[[i]]$est_N-trial[[i]]$sim_N)/trial[[i]]$sim_N
        decline_bias_comp[[i]] <- (trial[[i]]$est_decline-trial[[i]]$sim_decline)/
          trial[[i]]$sim_decline
        
        bias_N[i,c(1:3,5:6)] <- quantile(N_bias_comp[[i]],probs=c(0.1,0.25,0.5,0.75,0.9),
                                         na.rm=TRUE)
        bias_N[i,4] <- mean(N_bias_comp[[i]], na.rm=TRUE)
        bias_dec[i,c(1:3,5:6)] <- quantile(decline_bias_comp[[i]],
                                           probs=c(0.1,0.25,0.5,0.75,0.9), na.rm=TRUE)
        bias_dec[i,4] <- mean(decline_bias_comp[[i]], na.rm=TRUE)
      }
      
      outcomes <- c(rep("(a) Population decline",4),
                    rep("(b) Current population numbers",4))
      treatments <- rep(c("Sim Gravity - \n Est Gravity","Sim Gravity - \n Est Markov",
                          "Sim Markov - \n Est Gravity","Sim Markov - \n Est Markov"),2)
      df1 <- data.frame( outcome = outcomes,
                         treatment = treatments,
                         upper90 = c(bias_dec[,1],bias_N[,1]),
                         upper50 = c(bias_dec[,2],bias_N[,2]),
                         median = c(bias_dec[,3],bias_N[,3]),
                         mean = c(bias_dec[,4],bias_N[,4]),
                         lower50 = c(bias_dec[,5],bias_N[,5]),
                         lower90 = c(bias_dec[,6],bias_N[,6]))
      save( trial, df1, file=paste0(path,"/SimEval10years.RData"))
    }
    ggplot(data=df1,aes(x=treatment,y=median,ymin=lower90,ymax=upper90)) +
      geom_hline(yintercept=seq(-0.2,0.2,0.05),colour="grey80",lty=3) +
      geom_pointrange(fatten=6) +
      geom_errorbar(aes(x=treatment,y=median,ymin=lower50,ymax=upper50),width=0,lwd=2,
                    show.legend = FALSE) +
      theme_classic(base_size=18) +
      geom_hline(yintercept=0,lty=2) +
      facet_grid( outcome~., switch = "x") +
      theme(strip.background = element_blank(), strip.text = element_blank()) +
      geom_text(aes(label=outcome), x=.6, y=Inf, hjust=0, vjust=1, size=5) +
      ylab("Bias: (model predicted - simulated)/simulated") + 
      xlab("Simulation type - model type")
  } # tenyear
  else {
    if(load){
      load(file=paste0(path,"/SimEval.RData"))
    } else {
      # build the dataset with no tag loss
      TL <<- 0
      for( i in 1:4 ){
        trial[[i]] <- optimize( nrep, simMove=sim_MovMod[i], estMove=est_MovMod[i] )
        if( !is.finite( sum( trial[[i]]$obj$he() ) ) ) next
        N_bias_comp[[i]] <- (trial[[i]]$est_N-trial[[i]]$sim_N)/trial[[i]]$sim_N
        decline_bias_comp[[i]] <- (trial[[i]]$est_decline-trial[[i]]$sim_decline)/
          trial[[i]]$sim_decline
  
        bias_N[i,c(1:3,5:6)] <- quantile(N_bias_comp[[i]],probs=c(0.1,0.25,0.5,0.75,0.9),
                                         na.rm=TRUE)
        bias_N[i,4] <- mean(N_bias_comp[[i]], na.rm=TRUE)
        bias_dec[i,c(1:3,5:6)] <- quantile(decline_bias_comp[[i]],
                                           probs=c(0.1,0.25,0.5,0.75,0.9), na.rm=TRUE)
        bias_dec[i,4] <- mean(decline_bias_comp[[i]], na.rm=TRUE)
      }
      outcomes <- c(rep("(a) Population decline",4),
                    rep("(b) Current population numbers",4))
      treatments <- rep(c("Sim Gravity - \n Est Gravity","Sim Gravity - \n Est Markov",
                          "Sim Markov - \n Est Gravity","Sim Markov - \n Est Markov"),2)
      df1 <- data.frame( tagloss = rep("No tag loss",8),
                         outcome = outcomes,
                         treatment = treatments,
                         upper90 = c(bias_dec[,1],bias_N[,1]),
                         upper50 = c(bias_dec[,2],bias_N[,2]),
                         median = c(bias_dec[,3],bias_N[,3]),
                         mean = c(bias_dec[,4],bias_N[,4]),
                         lower50 = c(bias_dec[,5],bias_N[,5]),
                         lower90 = c(bias_dec[,6],bias_N[,6]))
      save( trial, df1, file=paste0(path,"/SimEval.RData"))
    } # load==FALSE
    # then: build the dataset WITH tag loss
    TL <<- 0.00204
    for( i in 1:4 ){
      trial[[i]] <- optimize( nrep, simMove=sim_MovMod[i], estMove=est_MovMod[i] )
      if( !is.finite( sum( trial[[i]]$obj$he() ) ) ) next
      N_bias_comp[[i]] <- (trial[[i]]$est_N-trial[[i]]$sim_N)/trial[[i]]$sim_N
      decline_bias_comp[[i]] <- (trial[[i]]$est_decline-trial[[i]]$sim_decline)/
        trial[[i]]$sim_decline
      
      bias_N[i,c(1:3,5:6)] <- quantile(N_bias_comp[[i]],probs=c(0.1,0.25,0.5,0.75,0.9),
                                       na.rm=TRUE)
      bias_N[i,4] <- mean(N_bias_comp[[i]], na.rm=TRUE)
      bias_dec[i,c(1:3,5:6)] <- quantile(decline_bias_comp[[i]],
                                         probs=c(0.1,0.25,0.5,0.75,0.9), na.rm=TRUE)
      bias_dec[i,4] <- mean(decline_bias_comp[[i]], na.rm=TRUE)
    }
    outcomes <- c(rep("(a) Population decline",4),
                  rep("(b) Current population numbers",4))
    treatments <- rep(c("Sim Gravity - \n Est Gravity","Sim Gravity - \n Est Markov",
                        "Sim Markov - \n Est Gravity","Sim Markov - \n Est Markov"),2)
    df2 <- data.frame( tagloss = rep("Tag loss",8),
                       outcome = outcomes,
                       treatment = treatments,
                       upper90 = c(bias_dec[,1],bias_N[,1]),
                       upper50 = c(bias_dec[,2],bias_N[,2]),
                       median = c(bias_dec[,3],bias_N[,3]),
                       mean = c(bias_dec[,4],bias_N[,4]),
                       lower50 = c(bias_dec[,5],bias_N[,5]),
                       lower90 = c(bias_dec[,6],bias_N[,6]))
  }
  
  df3 <- rbind( df1, df2 )
  
  ggplot(data=df3,aes(x=treatment,y=median,ymin=lower90,ymax=upper90,colour=tagloss)) +
    geom_hline(yintercept=seq(-0.25,0.25,0.05),colour="grey80",lty=3) +
    geom_pointrange(fatten=6,position=position_dodge(0.3)) +
    geom_errorbar(aes(x=treatment,y=median,ymin=lower50,ymax=upper50),width=0,lwd=2,
                  position=position_dodge(0.3),show.legend = FALSE) +
    geom_errorbar(aes(x=treatment,y=median,ymin=lower90,ymax=upper90),width=0,lwd=1,
                  position=position_dodge(0.3),show.legend = FALSE) +
    theme_classic(base_size=18) +
    geom_hline(yintercept=0,lty=2) +
    facet_grid( outcome~., switch = "x") +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    geom_text(aes(label=outcome), colour='grey40', x=.6, y=Inf,  hjust=0, vjust=1, size=5,show.legend = FALSE) +
    ylab("Bias: (model predicted - simulated)/simulated") + 
    xlab("Simulation type - model type")
}

"plotTrend" <- function(estMove, plottrend=c( "trend", "compare")){
  
  if(plottrend == "trend"){
    posterior <- PostPred(estMove)
    N        <- posterior$N
    Npost    <- data.frame()
    empty    <- vector()
    for( y in 1:ny ){
      Nsum       <- rowSums( N[,y,] )/1000
      Npost      <- rbind( Npost,
                           c( y+2000, 0, quantile( Nsum, probs=c(0.05,0.1,0.5,0.9,0.95))))
      for( a in 1:na ){
        empty[1]   <- y+2000
        empty[2]   <- a
        empty[3:7] <- quantile( N[,y,a], probs=c(0.05,0.1,0.5,0.9,0.95) )/1000
        Npost      <- rbind( Npost, empty )
      } # a
    } # y
    Npost$label    <- rep( c( "Total",
                         "Area 1: Lower Nechako and Stuart River",
                         "Area 2: Mid-Nechako",
                         "Area 3: Upper core area",
                         "Area 4: Upper Nechako"),
                         ny )
    names( Npost ) <- c( "Year", "Area", 
                         "lower5", "lower10", "median", "upper90", "upper95",
                         "label" )
    bounds            <- matrix( nrow=na+1, ncol=3 )
    bounds[,1]        <- 0:na
    bounds[1,2]       <- median( rowSums( N[,1,] ))/1000
    bounds[1+1:na,2]  <- apply( N[,1,]/1000, 2, median )
    bounds[1,3]       <- median( rowSums( N[,ny,] ))/1000
    bounds[1+1:na,3]  <- apply( N[,ny,]/1000, 2, median )
    bounds            <- as.data.frame( bounds )
    names(bounds)     <- c('Area','Initial','Final')
    
    trend <- ggplot( Npost, aes( x=Year, y=median ) ) +
      facet_grid( Area~. ) +
      geom_ribbon(aes(ymin = lower10, ymax = upper90), alpha = .5,
                  fill = "darkseagreen3", color = "transparent") +
      geom_line(color = "aquamarine4", lwd = .7) +
      theme_classic(base_size=18) +
      theme(strip.background = element_blank(), 
            strip.text = element_blank(),
            panel.border = element_rect(colour="black",fill=NA)) +
      geom_text(aes(label=label), x=ny+2000, y=Inf, hjust=1, vjust=2.4, size=4.5) +
      scale_y_continuous( breaks = seq(0,1,0.5) ) +
      ylab("Numbers (thousands)") + 
      xlab("Year") + 
      geom_hline(data=bounds, aes( yintercept=Initial ),
                 linetype=2, col="red") +
      geom_hline(data=bounds, aes( yintercept=Final ),
                 linetype=2, col="red")
    return(trend)
  } #if(plottrend)
  
  if( plottrend== "compare" ){
    mark     <- PostPred(estMove="Markov")
    grav     <- PostPred(estMove="Grav")
    N.mark   <- mark$N
    N.grav   <- grav$N
    npost    <- dim(N.grav)[1]
    Mpost    <- array(dim=c(npost,2))
    Mpost[,1]<- grav$M
    Mpost[,2]<- mark$M
    Nsum <- dec <- ycrit <- array( dim=c(npost,2))
    Nsum[,1] <- rowSums( N.grav[,ny,] )
    Nsum[,2] <- rowSums( N.mark[,ny,] )
    dec[,1]  <- rowSums( N.grav[,ny,] ) / rowSums( N.grav[,1,] )
    dec[,2]  <- rowSums( N.mark[,ny,] ) / rowSums( N.mark[,1,] )
    ycrit[,1]<- -log(200/rowSums(N.grav[,ny,]) )/Mpost[,1]
    ycrit[,2]<- -log(200/rowSums(N.mark[,ny,]) )/Mpost[,2]
    var <- rep( c( "N2023","decline","M","ycrit" ),2 )
    mod <- rep(c("Gravity","Markov"), each = 4)
    postdf   <- data.frame()
    for( i in 1:2 ){
      postdf  <- rbind( postdf, 
                        c( quantile( Nsum[,i], probs=c(0.1,0.25,0.5,0.75,0.90)), mean( Nsum[,i] )),
                        c( quantile( dec[,i], probs=c(0.1,0.25,0.5,0.75,0.90)), mean( dec[,i] )),
                        c( quantile( Mpost[,i], probs=c(0.1,0.25,0.5,0.75,0.90)), mean( Mpost[,i] )),
                        c( quantile( ycrit[,i], probs=c(0.1,0.25,0.5,0.75,0.90)), mean( ycrit[,i] )))
    } # i
    postdf$var <- var
    postdf$mod <- mod
    names(postdf) <- c( "lower10", "lower25", "median", "upper75", "upper90",
                        "mean","Variable","Model" )
    
    df1 <- postdf %>% filter(Variable=="N2023")
    p1 <- ggplot(data=df1,aes(x=Model,y=median,ymin=lower10,ymax=upper90)) +
      geom_pointrange(fatten=6) +
      geom_errorbar(aes(x=Model,y=median,ymin=lower25,ymax=upper75),width=0,lwd=2) +
      theme_classic(base_size=18) +
      scale_y_continuous(breaks=seq(450,625,25)) +
      geom_hline(yintercept=seq(450,625,25), lty=3, col="grey"  ) +
      ylab(expression(paste("N"[2022]))) + xlab("")
    
    df2 <- postdf %>% filter(Variable=="decline")
    p2 <- ggplot(data=df2,aes(x=Model,y=median,ymin=lower10,ymax=upper90)) +
      geom_pointrange(fatten=6) +
      geom_errorbar(aes(x=Model,y=median,ymin=lower25,ymax=upper75),width=0,lwd=2) +
      theme_classic(base_size=18) +
      geom_hline(yintercept=seq(0.40,0.6,0.05), lty=3, col="grey" ) +
      ylab(expression(paste("Decline (N"[2022]," / N"[2001],")"))) + xlab("")
    
    df3 <- postdf %>% filter(Variable=="M")
    p3 <- ggplot(data=df3,aes(x=Model,y=median,ymin=lower10,ymax=upper90)) +
      geom_pointrange(fatten=6) +
      geom_errorbar(aes(x=Model,y=median,ymin=lower25,ymax=upper75),width=0,lwd=2) +
      theme_classic(base_size=18) +
      geom_hline(yintercept=seq(0.025,0.045,0.005), lty=3, col="grey"  ) +
      ylab("Natural mortality rate") + xlab("")
    
    df4 <- postdf %>% filter(Variable=="ycrit")
    p4 <- ggplot(data=df4,aes(x=Model,y=median,ymin=lower10,ymax=upper90)) +
      geom_pointrange(fatten=6) +
      geom_errorbar(aes(x=Model,y=median,ymin=lower25,ymax=upper75),width=0,lwd=2) +
      theme_classic(base_size=18) +
      geom_hline(yintercept=seq(20,45,5), lty=3, col="grey"  ) +
      ylab("Years until N < 200") + xlab("")
    
    x.grob   <- textGrob( "Model", gp=gpar( fontsize=15 ), just="bottom" )
    compare <- grid.arrange(p1, p2, p3, p4, nrow=2,ncol=2, bottom=x.grob)
    return(compare)

  } # if( compare )
  
} #plotTrend

"compareM" <- function(){
  movemod   <- rep(c("Grav","Markov"),4)
  Mmean     <- c( 0.04, 0.04, 0.02, 0.02, 0.08, 0.08, 0.04, 0.04 )
  Msd       <- c( rep( 0.4, 6 ), rep( 0.8, 2 ))
  postmat   <- matrix(nrow=15000,ncol=8)
  priormat  <- matrix(nrow=15000,ncol=4)
  for( i in 1:8 ){
    load(paste0("Posterior ",movemod[i],Mmean[i],Msd[i],".RData"))
    postmat[,i]    <- exp( as.matrix(fit)[,'lnM'] )
  } #i
  for( i in 1:4 ){
    priormat[,i] <- exp(rnorm(15000, log(Mmean[(i)*2]), Msd[(i)*2] ))
  }
  
  postdf         <- as.data.frame( postmat )
  names(postdf)  <- c("Grav0.040.4","Markov0.040.4",
                     "Grav0.020.4","Markov0.020.4",
                     "Grav0.080.4","Markov0.080.4",
                     "Grav0.041","Markov0.041")
  priordf        <- as.data.frame( priormat )
  
  names(priordf) <- c("Prior0.040.4",
                     "Prior0.020.4",
                     "Prior0.080.4",
                     "Prior0.041")
  
#  fulldf         <- cbind(priordf,postdf)
  
  label          <- c( "(a) Default",
                       "(b) Low M",
                       "(c) High M",
                       "(d) Vague M")
  
  priordf$label  <- label

  df1         <- postdf[,1:2]
  names(df1)  <- c("Gravity","Markov")
  df1         <- melt( df1 )
  p1 <- ggplot(data=priordf, aes(x=Prior0.040.4)) + 
    geom_density( adjust=2, colour="lightblue", fill="lightblue") +
    theme_classic(base_size=18) +
    ylab("") + xlab("") + xlim(c(0,0.1)) + ylim(c(0,80)) +
    geom_text(aes(label=label[1]), x=0, y=Inf, hjust=0, vjust=1, size=4.5) +
    theme(legend.position = "none", 
          strip.background = element_blank(), 
          strip.text = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0.4,0.4,-0.3,-0.3),'cm'),
          #panel.border = element_rect(colour="black",fill=NA),
          panel.grid.minor = element_line( size=0.25, linetype="solid", 
                                           colour="grey"),
          panel.grid.major = element_line( size=0.25, linetype="solid", 
                                           colour="grey")) +
    geom_density(data=df1,aes( x=value, linetype=variable))

  df2         <- postdf[,3:4]
  names(df2)  <- c("Gravity","Markov")
  df2         <- melt( df2 )
  p2 <- ggplot(data=priordf, aes(x=Prior0.020.4)) + 
    geom_density( adjust=2, colour="lightblue", fill="lightblue" ) +
    theme_classic(base_size=18) +
    ylab("") + xlab("") + xlim(c(0,0.1)) + ylim(c(0,80)) +
    geom_text(aes(label=label[2]), x=0, y=Inf, hjust=0, vjust=1, size=4.5) +
    theme(legend.position = c(0.95,.95), legend.justification =c("right","top"),
          legend.title=element_blank(),
          strip.background = element_blank(), 
          strip.text = element_blank(),
          axis.text = element_blank(),
          plot.margin = unit(c(0.4,0.4,-0.3,-0.3),'cm'),
          #panel.border = element_rect(colour="black",fill=NA),
          panel.grid.minor = element_line( size=0.25, linetype="solid", 
                                           colour="grey"),
          panel.grid.major = element_line( size=0.25, linetype="solid", 
                                           colour="grey")) +
    geom_density(data=df2,aes( x=value, linetype=variable))
    
  df3         <- postdf[,5:6]
  names(df3)  <- c("Gravity","Markov")
  df3         <- melt( df3 )
  p3 <- ggplot(data=priordf, aes(x=Prior0.080.4)) + 
    geom_density( adjust=2, colour="lightblue", fill="lightblue" ) +
    theme_classic(base_size=18) +
    ylab("") + xlab("") + xlim(c(0,0.1)) + ylim(c(0,80)) +
    geom_text(aes(label=label[3]), x=0, y=Inf, hjust=0, vjust=1, size=4.5) +
    theme(legend.position = "none", 
          strip.background = element_blank(), 
          strip.text = element_blank(), 
          plot.margin = unit(c(0.4,0.4,-0.3,-0.3),'cm'),
          #panel.border = element_rect(colour="black",fill=NA),
          panel.grid.minor = element_line( size=0.25, linetype="solid", 
                                           colour="grey"),
          panel.grid.major = element_line( size=0.25, linetype="solid", 
                                           colour="grey")) +
    geom_density(data=df3,aes( x=value, linetype=variable))
  
  df4         <- postdf[,7:8]
  names(df4)  <- c("Gravity","Markov")
  df4         <- melt( df4 )
  p4 <- ggplot(data=priordf, aes(x=Prior0.041)) + 
    geom_density( adjust=2, colour="lightblue", fill="lightblue" ) +
    theme_classic(base_size=18) +
    ylab("") + xlab("") + xlim(c(0,0.1)) + ylim(c(0,80)) +
    geom_text(aes(label=label[4]), x=0, y=Inf, hjust=0, vjust=1, size=4.5) +
    theme(legend.position = "none", 
          strip.background = element_blank(), 
          strip.text = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0.4,0.4,-0.3,-0.3),'cm'),
          #panel.border = element_rect(colour="black",fill=NA),
          panel.grid.minor = element_line( size=0.25, linetype="solid", 
                                           colour="grey"),
          panel.grid.major = element_line( size=0.25, linetype="solid", 
                                           colour="grey")) +
    geom_density(data=df4,aes( x=value, linetype=variable))
    
  x.grob   <- textGrob( "Natural mortality rate, M", 
                        gp=gpar( fontsize=15 ), just="bottom",vjust=0,hjust=0.25 )
  y.grob   <- textGrob( "Probability density", rot=90, hjust=0.25, vjust=1,
                        gp=gpar( fontsize=15 ), just="bottom" )
  compare.M <- grid.arrange(p1, p2, p3, p4, nrow=2,ncol=2, bottom=x.grob, left=y.grob,
               widths=c(3.4,3.1), heights=c(3.15,3.4)) 
  return(compare.M)
    
} #compareM()
