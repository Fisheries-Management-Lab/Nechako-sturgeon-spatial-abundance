// spatial mark-recapture model
#include <TMB.hpp>

template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}

template<class Type>
Type objective_function<Type>::operator() () {
  // data:
  DATA_INTEGER(est);                // flag for whether you are estimating real data or simulated (real needs prior on F)
  DATA_INTEGER(fullmov);            // flag for whether you are estimating gravity or Markov movement
  DATA_SCALAR(CV);                  // CV in catch
  DATA_SCALAR(Mmu);                 // mean of M prior
  DATA_SCALAR(Msd);                 // sd of M prior
  DATA_IMATRIX(tagind);             // indices of area and year for tagdat
  DATA_VECTOR(tagdat);              // tag data file
//  DATA_MATRIX(nocap);               // number of tags that were never recaptured in a particular release year and area
  DATA_MATRIX(cat);                 // number of initial tagging events (first capture)
  DATA_MATRIX(cattot);              // all catch - tagged and untagged
  DATA_MATRIX(tel);                 // telemetry data - movement
  DATA_IMATRIX(Find);               // F-index - allows F by area and time as a vector

  // parameters:
  PARAMETER(lnN0);                  // initial abundance
  PARAMETER(lnM);                   // instantaneous natural mortality
  PARAMETER_VECTOR(lnF);            // instantaneous capture probability
  PARAMETER_VECTOR(movest);         // estimated movement parameters
   
  // procedures: (transformed parameters)
  Type N0 = exp(lnN0);              // initial abundance
  Type M = exp(lnM);                // natural mortality rate

  Type tiny = 0.00001;              // small number
  Type pen = 0.0;                   // penalty for negative numbers 
  int ny = cat.rows();              // number of years
  int na = cat.cols();              // number of ages
  int ntag = tagdat.rows();         // number of tags
  int nMP = movest.size();          // number of estimated movement parameters
  int nF = Find.rows();             // number of recaptures
  int tagdim = na*ny+1;             // dimensions of tag matrices
//  int nnocap = nocap.size();        // number of no recaptures
  Type one=1.;
  Type zero=0.;
  vector<Type> MR(na);              // mark rate by area
  vector<Type> idist(na);           // initial distribution of fish based on transition matrix
  vector<Type> Frow(na);            // annual sampling rate
  vector<Type> Nrow(na);            // annual total abundance
  vector<Type> NMrow(na);           // annual unmarked abundance before sampling
  vector<Type> NTrow(na);           // annual unmarked abundance after sampling
  vector<Type> tempS(na);           // temporary survival placeholder
  matrix<Type> movTrans(na,na);     // transposed movement matrix
  matrix<Type> F(ny,na);            // sampling rate
  matrix<Type> N(ny,na);            // number of fish by area and year
  matrix<Type> NT(ny,na);           // number of (UN)tagged fish by area and year (before sampling) ? (BRETT IS CHANGING THIS TO MARKED FISH)
  matrix<Type> NM(ny,na);           // number of tagged fish by area and year (after sampling) ? (BRETT IS CHANGING THIS TO UNTAGGED FISH)
  matrix<Type> catpred(ny,na);      // predicted recaptures
  matrix<Type> nocappred(ny,na);    // predicted new captures
  matrix<Type> Pnocap(ny,na);       // probability of never seeing a fish again
  matrix<Type> mov(na,na);          // gravity matrix
  matrix<Type> emov(na,na);         // exponentiated gravity matrix
  matrix<Type> obsCat(tagdim,tagdim);// observed number of captures
  matrix<Type> pCat(tagdim,tagdim);  // predicted probability of capture
  array<Type> Surv(na,ny,ny,na);    // survival probability array
  array<Type> Pcap(na,ny,ny,na);    // capture probability array
  array<Type> cappred(na,ny,ny,na); // predicted total captures
  
  Type nll = 0.;                    // initialize negative log likelihood
  F.setZero();
  N.setZero();
  Pcap.setZero();
  Surv.setZero();
  obsCat.setZero();

  // INITIALIZE VARIABLES
  for( int i=0; i<nF; i++ ){
    int y=Find(i,0);
    int a=Find(i,1);
    F(y,a) = exp(lnF(i));           // place estimated Fs into F matrix
  } // i
  
  // set movement matrix
  if( fullmov == 0 ){
    // Gravity model calculation
    for( int a=0; a<na; a++ ){
      mov(a,0) = 0.;                // first row is set to zero
    } // a
    for( int a=0; a<na; a++ ){
      for( int a2=1; a2<na; a2++ ){
        mov(a,a2) = movest(a2-1);   // gravity weight for the receiving area
      } // a2
      mov(a,a) += exp(movest(nMP-1)); // viscosity to staying in the same area
    } // a
  } // if(!fulmov)

  if( fullmov > 0 ){
    // fully described Markov model
    for( int a=0; a<na; a++ ){
      mov(a,0) = 0.;                // first row is set to zero
    } // a
    int ii=0;
    
    for( int a2=1; a2<na; a2++ ){
      for( int a=0; a<na; a++ ){
        mov(a,a2)=movest(ii);       // probability of moving to cell is proportional to its estimated G
        ii++;
      } // a
    } // a2
  } // if(fullmov)

  for( int a=0; a<na; a++ ){
    for( int a2=0; a2<na; a2++){
      mov(a,a2) = exp(mov(a,a2));
    }
    vector<Type> mov_a = mov.row(a);
    mov.row(a) = mov_a/sum(mov_a);  // make rows of movement matrix proper probabilities
  } // a
    //std::cout<<"mov\n"<<mov<<"\n";
     
  // initialize initial distribution to uniform across areas
  for( int a=0; a<na; a++ ){
    idist(a) = one/na;
  } // a
  
  // numerically solve for initial distribution
  movTrans = mov.transpose();
  for( int i=0; i<400; i++ ){
    idist = movTrans * idist;              // numerically solve equilibrium initial distribution
  } // i
    //std::cout<<"idist\n"<<idist<<"\n";

  // POPULATION DYNAMICS
  N.row(0) = N0 * idist;                // distribute initial abundance (marked and unmarked) across areas
  NM.row(0) = N0 * idist;               // unmarked fish prior to tagging
  Frow = F.row(0);
  NT.row(0) = N0 * idist * exp(-Frow);  // unmarked fish immediately after tagging (removes newly tagged fish)
  
  // redistribute fish through years and areas and remove fish due to mortality
  for( int y=1; y<ny; y++){
    // set abundance across areas to that from last year
    Nrow = N.row(y-1);
    NTrow = NT.row(y-1);      // fish about to be marked is based on those that weren't captured previously 
    // BRETT this is different than tpl code - check to make sure it makes sense
    // survive them
    Nrow = movTrans * Nrow;    // total pop: move first then survive
    N.row(y) = Nrow * exp(-M); 
    NMrow = movTrans * NTrow;  // unmarked pop: move unmarked first then survive
    NM.row(y) = NMrow * exp(-M);
    Frow = F.row(y);
    NTrow = NTrow * exp(-M-Frow); // unmarked after capture: survive first then move
    NT.row(y) = movTrans * NTrow;
  } // y
    //std::cout<<"N\n"<<N<<"\n";

    //std::cout<<"F\n"<<F<<"\n";
  // SURVIVAL AND CAPTURE PROBABILITY
  // probability of surviving and being captured in state a
    //std::cout<<"movTrans\n"<<movTrans<<"\n";
  for( int a=0; a<na; a++ ){
    for( int y=0; y<ny; y++ ){
      // probability of capture and survival in  year of marking
      Surv(a,y,y,a) = one; // exp(-M/Type(2.));
      Pcap(a,y,y,a) = zero; // exp(-M+F(y,a)/2)*(one-exp(-F(y,a)/Type(2.)));
    } // y
  } // a
  
  for( int a=0; a<na; a++ ){
    for( int y=0; y<ny; y++ ){
      for( int y2=y+1; y2<ny; y2++ ){
        for( int a2=0; a2<na; a2++ ){
          tempS(a2) = Surv(a,y,y2-1,a2)*exp(-M);
          //std::cout<<"Surv\t"<<Surv(a,y,y2-1,a2)<<"\t"<<F(y2-1,a2)<<"\t"<<exp(-M-F(y2-1,a2))<<"\n";
        }
        //Frow = F.row(y2);
          //std::cout<<"a1 "<<a<<"\ty1 "<<y<<"\ty2 "<<y2<<"1st tempS\n"<<tempS<<"\n";
        tempS = movTrans * tempS;   // matrix multiply
        for( int a2=0; a2<na; a2++ ){
//          Surv(a,y,y2,a2) = tempS(a2);
          Pcap(a,y,y2,a2) = tempS(a2) * (one-exp(-F(y2,a2)));
          Surv(a,y,y2,a2) = tempS(a2) * exp(-F(y2,a2));
          //std::cout<<"Pcap\t"<<Surv(a,y,y2,a2)<<"\t"<<Pcap(a,y,y2,a2)<<"\n";
        }
        //Pcap.col(a).col(y).col(y2) = Surv.col(a).col(y).col(y2) * (one-exp(-Frow));  // element-wise multiply
      }
    }
  }
    //std::cout<<"Pcap\n"<<Pcap<<"\n";
  
  // CALCULATE LIKELIHOODS
  Type thetaC = 0.;
  Type thetaT = 0.;
  Type thetaY = 0.;
  Type prior = 0.;
  
  // Theta(C) - initial captures
  for( int i=0; i<nF; i++ ){
    int y = Find( i, 0 );
    int a = Find( i, 1 );
    
    catpred(y,a) = N(y,a) * ( one - exp(-F(y,a)) );
    //thetaC -= dbinom( cat(y,a), NM(y,a), one-exp(-F(y,a)) );
    thetaC -= dnorm( log(cattot(y,a)), log(catpred(y,a)), CV, true );
  } // i
    //std::cout<<"cat\n"<<cat<<"\ncatpred\n"<<catpred<<"\nthetaC\t"<<thetaC<<"\n";
  
  // Theta(T) - all tags - recaptured and never seen again
  // set pcapture to anything so blank cells don't get NA when logged
  for( int i=0; i<tagdim; i++ ){
    for( int ii=0; ii<tagdim; ii++ ){
      pCat(i,ii) = Type(1e-10); 
    } // ii
  } // i
  
  for( int i=0; i<ntag; i++ ){
    int tagi = tagind(i,0);
    int a = tagind(i,1);
    int y = tagind(i,2);
    int reci = tagind(i,3);
    int y2 = tagind(i,4);
    int a2 = tagind(i,5);
    
    obsCat( reci, tagi ) += tagdat(i);
    if( reci < tagdim-1 ){  // if this WAS recaptured
      pCat(reci,tagi) = Pcap(a,y,y2,a2);
    } else {  // if never seen again
      pCat(reci,tagi) = one - sum( Pcap.rotate(2).col(y).col(a) );
    } // if(reci<tagdim-1)
  } // i
  vector<Type> obsCatvec(tagdim);
  vector<Type> pCatvec(tagdim);
  for( int i=0; i<tagdim; i++ ){
    obsCatvec = obsCat.col(i);
    pCatvec = pCat.col(i);
     //std::cout<<"tagid "<<i<<"obsCatvec\n"<<obsCatvec<<"\npCatvec\n"<<pCatvec<<"\n";
    thetaT -= dmultinom( obsCatvec, pCatvec, true );
  } // i

  // Theta(Y) - telemetry observations
  vector<Type> movvec(na*na);  // create whole vector of move estimates so multinomial scales properly
  vector<Type> telvec(na*na);  // create whole vector of tel estimates so multinomial scales properly
  int ii=0;
  for( int a=0; a<na; a++ ){
    for( int a2=0; a2<na; a2++ ){
      movvec(ii) = mov(a,a2);
      telvec(ii) = tel(a,a2);
      ii += 1;
    } // a2
  } // a
  thetaY -= dmultinom( telvec, movvec, true );
    //std::cout<<"tel\n"<<tel<<"\nmov\n"<<mov<<"\nthetaY\t"<<thetaY<<"\n";
  
  // Prior on M
  prior -= dnorm( lnM, Mmu, Msd, true );
  prior -= dnorm( movest, 0, 5, true ).sum();
  if(est==1)
    prior -= dnorm( lnF, -5, 5, true ).sum();
  //std::cout<<"prior\n"<<prior<<"\n";
  
  //std::cout<<"nll is turned off\n";
  nll = thetaC + thetaT + thetaY + prior;

  ADREPORT(N0);                   // Initial abundance
  ADREPORT(M);                    // Instantaneous natural mortality
  ADREPORT(F);                    // Instantaneous sampling rate
  ADREPORT(mov);                  // movement rates
  ADREPORT(N);
  
  REPORT(N);
  REPORT(NM);
  REPORT(NT);

  return nll;
}
