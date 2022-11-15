comparebcJ <-function (marker1, marker2, D, alpha, plots){

  if ((length(marker1) != length(D)) | (length(marker2) != length(D))) {
    stop("ERROR: The length of the 'marker' and 'D' inputs must be equal.")
  } else if (min(D) != 0 | max(D) != 1) {
    stop("ERROR: Controls must be assigned a value of 0; cases must be assigned a value of 1. Both controls and cases should be included in the dataset.")
  } else if (alpha <= 0 | alpha >= 1) {
    stop("ERROR: The level of significance, alpha, should be set between 0 and 1. A common choice is 0.05.")
  } else if (sum(is.na(marker1)) > 0 | sum(is.na(marker2)) > 0 | sum(is.na(D)) > 0) {
    stop("ERROR: Please remove all missing data before running this function.")
  } else {

    erf <- function (x) 2 * pnorm(x * sqrt(2)) - 1
    erfc <- function (x) 2 * pnorm(x * sqrt(2), lower.tail = FALSE)
    erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
    erfcinv <- function (x) qnorm(x/2, lower.tail = FALSE)/sqrt(2)

    pmalpha=qnorm(1-alpha/2);

    if (plots!="on"){plots="off"}


    W1a=marker1[D==0]
    W1b=marker1[D==1]
    W2a=marker2[D==0]
    W2b=marker2[D==1]




    Wa=cbind(W1a,W2a);
    Wb=cbind(W1b,W2b);
    na=length(W1a)
    nb=length(W1b)

    likbox2D<-function(W1a,W1b,W2a,W2b,lam){
      na=length(W1a);
      nb=length(W1b);
      out=c();
      #for (i in 1:length(h)){
      W1alam= (W1a^lam[1]-1)/lam[1];
      W2alam= (W2a^lam[2]-1)/lam[2];
      W1blam= (W1b^lam[1]-1)/lam[1];
      W2blam= (W2b^lam[2]-1)/lam[2];
      Walam=cbind(W1alam,W2alam)
      Wblam=cbind(W1blam,W2blam)
      #cova=1/na*sum((W1alam-mean(W1alam))*(W2alam-mean(W2alam)))
      #covb=1/nb*sum((W1blam-mean(W1blam))*(W2blam-mean(W2blam)))

      cova= cov(Walam)*(na-1)/na;
      covb= cov(Wblam)*(nb-1)/nb;
      mualam=c(mean(W1alam), mean(W2alam))
      mublam=c(mean(W1blam), mean(W2blam))

      #dmvn(c(0,0), mu, mcov, log = F)

      #   out=          (-sum(log(dmvn(Walam,mualam,cova)))-sum(log(dmvn(Wblam,mublam,covb))) -(lam[1]-1)*sum(log(W1a))-(lam[2]-1)*sum(log(W2a))-(lam[1]-1)*sum(log(W1b))-(lam[2]-1)*sum(log(W2b)));
      out=          -sum(log(dmvnorm(Walam,mualam,cova)))-sum(log(dmvnorm(Wblam,mublam,covb))) -(lam[1]-1)*sum(log(W1a))-(lam[2]-1)*sum(log(W2a))-(lam[1]-1)*sum(log(W1b))-(lam[2]-1)*sum(log(W2b));

      #out1=-sum(log(dmvnorm(Walam,mualam,cova)))
      #out2=-sum(log(dmvnorm(Wblam,mublam,covb)))
      #out3=-(lam[1]-1)*sum(log(W1a))
      #out4=-(lam[2]-1)*sum(log(W2a))
      #out5=-(lam[1]-1)*sum(log(W1b))
      #out6=-(lam[2]-1)*sum(log(W2b))

      #out=c(out,out1,out2,out3,out4,out5,out6, cova, covb)
      #out=          (-sum(log(mvnpdf2([W1alam(g(1)) W2alam(g(2))],[mean(W1alam(g(1))) mean(W2alam(g(2)))],[std(W1alam(g(1)),1).^2   cova;cova  (std(W2alam(g(2)),1)).^2])))+...
      #               -sum(log(mvnpdf2([W1blam(g(1)) W2blam(g(2))],[mean(W1blam(g(1))) mean(W2blam(g(2)))],[std(W1blam(g(1)),1).^2   covb;covb  (std(W2blam(g(2)),1)).^2])))+...
      #               -(g(1)-1).*sum(log(W1a))-(g(2)-1).*sum(log(W2a))-(g(1)-1).*sum(log(W1b))-(g(2)-1).*sum(log(W2b)));

      return(out)
    }

    lam=c(2,2)
    likbox2D(W1a,W1b,W2a,W2b,lam)

    logL<-function(h){
      likbox2D(W1a,W1b,W2a,W2b,h)
    }

    # init=c(0.8,0.8)
    # lam=fminsearch(logL,init)
    # lam=c(lam$optbase$xopt)
    lam<- optim(c(1,1),logL,gr=NULL,method="BFGS", control=list(maxit=10000))
    lam=c(lam$par)
    lam

    #init=c(-10,10)
    #lam=optimize(logL,init)
    #lam

    #f <- function (x, a) (x - a)^2
    #xmin <- optimize(f, c(0, 1), tol = 0.0001, a = 1/3)
    #xmin

    out <- optim(c(1,1), logL, method = "Nelder-Mead")
    lam = out$par


    ###################################
    W1alam= (W1a^lam[1]-1)/lam[1];
    W2alam= (W2a^lam[2]-1)/lam[2];
    W1blam= (W1b^lam[1]-1)/lam[1];
    W2blam= (W2b^lam[2]-1)/lam[2];

    Walam=cbind(W1alam,W2alam)
    Wblam=cbind(W1blam,W2blam)

    m1ahat=mean(W1alam)
    m1bhat=mean(W1blam)
    m2ahat=mean(W2alam)
    m2bhat=mean(W2blam)

    s1ahat=sqrt( var(W1alam)*(na-1)/na )
    s1bhat=sqrt( var(W1blam)*(nb-1)/nb )
    s2ahat=sqrt( var(W2alam)*(na-1)/na )
    s2bhat=sqrt( var(W2blam)*(nb-1)/nb )

    cova= cov(Walam)*(na-1)/na;cova=cova[1,2];
    covb= cov(Wblam)*(nb-1)/nb;covb=covb[1,2];

    lam1=lam[1]
    lam2=lam[2]

    ##################################

    h1_1=sum(-2/(2*s1ahat^2 - (2*cova^2)/s2ahat^2))
    h1_2=sum(-(2*s1ahat*s2ahat^2*(cova*lam1 - lam2*s2ahat^2 - W2a^lam2*cova*lam1 + W1a^lam1*lam2*s2ahat^2 + cova*lam1*lam2*m2ahat - lam1*lam2*m1ahat*s2ahat^2))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h1_3=0
    h1_4=0;
    h1_5=sum(cova/(- cova^2 + s1ahat^2*s2ahat^2));
    h1_6=sum((s2ahat*(2*cova^2*(lam2 - W1a^lam1*lam2 + lam1*lam2*m1ahat) - 2*cova*s1ahat^2*(lam1 - W2a^lam2*lam1 + lam1*lam2*m2ahat)))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h1_7=0
    h1_8=0
    h1_9=sum((cova^2*(lam2*m2ahat - W2a^lam2 + 1) + s2ahat^2*((lam2*m2ahat - W2a^lam2 + 1)*s1ahat^2 - 2*cova*lam2*m1ahat))/(lam2*(s1ahat^2*s2ahat^2 - cova^2)^2) - (cova*s2ahat^2*(2*lam2 - 2*W1a^lam1*lam2))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h1_10=0
    h1_11=sum((2*s2ahat^2*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1))/(lam1^2*(- 2*cova^2 + 2*s1ahat^2*s2ahat^2)))
    h1_12=sum(-(cova*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1))/(lam2^2*(- cova^2 + s1ahat^2*s2ahat^2)))


    h2_1=sum(-(2*s1ahat*s2ahat^2*(cova*lam1 - lam2*s2ahat^2 - W2a^lam2*cova*lam1 + W1a^lam1*lam2*s2ahat^2 + cova*lam1*lam2*m2ahat - lam1*lam2*m1ahat*s2ahat^2))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h2_2=sum(- (3*s1ahat^2*s2ahat^6 + cova^4*(lam1^2*m2ahat^2 + lam1^2*s2ahat^2) - cova^3*(2*lam1*m2ahat*s2ahat^2 - 2*W1a^lam1*lam1*m2ahat*s2ahat^2 + 2*lam1^2*m1ahat*m2ahat*s2ahat^2) - cova*(6*lam1*m2ahat*s1ahat^2*s2ahat^4 - 6*W1a^lam1*lam1*m2ahat*s1ahat^2*s2ahat^4 + 6*lam1^2*m1ahat*m2ahat*s1ahat^2*s2ahat^4) + cova^2*(W1a^(2*lam1)*s2ahat^4 - 2*W1a^lam1*s2ahat^4 + s2ahat^4 + 2*lam1*m1ahat*s2ahat^4 + lam1^2*m1ahat^2*s2ahat^4 - 2*W1a^lam1*lam1*m1ahat*s2ahat^4 + 3*lam1^2*m2ahat^2*s1ahat^2*s2ahat^2) - 6*W1a^lam1*s1ahat^2*s2ahat^6 + 3*W1a^(2*lam1)*s1ahat^2*s2ahat^6 - lam1^2*s1ahat^4*s2ahat^6 + 3*lam1^2*m1ahat^2*s1ahat^2*s2ahat^6 + 6*lam1*m1ahat*s1ahat^2*s2ahat^6 - 6*W1a^lam1*lam1*m1ahat*s1ahat^2*s2ahat^6)/(lam1^2*(s1ahat^2*s2ahat^2 - cova^2)^3) - (cova^4*(W2a^(2*lam2)*lam1^2 - 2*W2a^lam2*lam1^2 + lam1^2) + cova^2*(3*lam1^2*s1ahat^2*s2ahat^2 + 3*W2a^(2*lam2)*lam1^2*s1ahat^2*s2ahat^2 - 6*W2a^lam2*lam1^2*s1ahat^2*s2ahat^2) - lam2*((2*W2a^lam2*lam1^2*m2ahat - 2*lam1^2*m2ahat)*cova^4 + (2*lam1*s2ahat^2 + 2*lam1^2*m1ahat*s2ahat^2 - 2*W1a^lam1*lam1*s2ahat^2 - 2*W2a^lam2*lam1*s2ahat^2 - 2*W2a^lam2*lam1^2*m1ahat*s2ahat^2 + 2*W1a^lam1*W2a^lam2*lam1*s2ahat^2)*cova^3 + (6*W2a^lam2*lam1^2*m2ahat*s1ahat^2*s2ahat^2 - 6*lam1^2*m2ahat*s1ahat^2*s2ahat^2)*cova^2 + (6*lam1*s1ahat^2*s2ahat^4 + 6*lam1^2*m1ahat*s1ahat^2*s2ahat^4 - 6*W1a^lam1*lam1*s1ahat^2*s2ahat^4 - 6*W2a^lam2*lam1*s1ahat^2*s2ahat^4 - 6*W2a^lam2*lam1^2*m1ahat*s1ahat^2*s2ahat^4 + 6*W1a^lam1*W2a^lam2*lam1*s1ahat^2*s2ahat^4)*cova))/(lam1^2*lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3))
    h2_3=sum(0)
    h2_4=sum(0)
    h2_5=sum((s1ahat*(2*cova^2*(lam1 - W2a^lam2*lam1 + lam1*lam2*m2ahat) - 2*cova*s2ahat^2*(lam2 - W1a^lam1*lam2 + lam1*lam2*m1ahat)))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h2_6=sum((lam2*((2*cova*s2ahat^3*(2*lam1 + 2*lam1^2*m1ahat - 2*W1a^lam1*lam1 - 2*W2a^lam2*lam1 + 2*W1a^lam1*W2a^lam2*lam1 - 2*W2a^lam2*lam1^2*m1ahat) - 2*cova*s2ahat*(4*cova*lam1^2*m2ahat - 4*W2a^lam2*cova*lam1^2*m2ahat))*s1ahat^3 + 2*cova*s2ahat*(2*cova^2*lam1 + 2*cova^2*lam1^2*m1ahat - 2*W1a^lam1*cova^2*lam1 - 2*W2a^lam2*cova^2*lam1 - 2*W2a^lam2*cova^2*lam1^2*m1ahat + 2*W1a^lam1*W2a^lam2*cova^2*lam1)*s1ahat) - 2*cova*s1ahat^3*s2ahat*(2*cova*lam1^2 + 2*W2a^(2*lam2)*cova*lam1^2 - 4*W2a^lam2*cova*lam1^2))/(lam1^2*lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3) - ((4*cova^2*lam1^2*m2ahat^2*s2ahat - 2*cova*s2ahat^3*(2*lam1*m2ahat + cova*lam1^2 + 2*lam1^2*m1ahat*m2ahat - 2*W1a^lam1*lam1*m2ahat))*s1ahat^3 + (2*cova*(2*cova + 2*W1a^(2*lam1)*cova - 4*W1a^lam1*cova + 2*cova*lam1^2*m1ahat^2 + 4*cova*lam1*m1ahat - 4*W1a^lam1*cova*lam1*m1ahat)*s2ahat^3 + 2*cova*(cova^3*lam1^2 - 2*cova^2*lam1*m2ahat + 2*W1a^lam1*cova^2*lam1*m2ahat - 2*cova^2*lam1^2*m1ahat*m2ahat)*s2ahat)*s1ahat)/(lam1^2*(s1ahat^2*s2ahat^2 - cova^2)^3))
    h2_7=sum(0)
    h2_8=sum(0)
    h2_9=sum(- (cova*(s2ahat^4*(- 4*lam2^2*m1ahat^2*s1ahat + 2*lam2^2*s1ahat^3) - s1ahat^3*s2ahat^2*(4*lam2*m2ahat - 4*W2a^lam2 + 2*W2a^(2*lam2) + 2*lam2^2*m2ahat^2 - 4*W2a^lam2*lam2*m2ahat + 2)) - cova^3*(s1ahat*(4*lam2*m2ahat - 4*W2a^lam2 + 2*W2a^(2*lam2) + 2*lam2^2*m2ahat^2 - 4*W2a^lam2*lam2*m2ahat + 2) + 2*lam2^2*s1ahat*s2ahat^2) + s1ahat^3*s2ahat^4*(2*lam2*m1ahat + 2*lam2^2*m1ahat*m2ahat - 2*W2a^lam2*lam2*m1ahat) + cova^2*s1ahat*s2ahat^2*(6*lam2*m1ahat + 6*lam2^2*m1ahat*m2ahat - 6*W2a^lam2*lam2*m1ahat))/(lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3) - (lam1*((6*lam2 + 6*lam2^2*m2ahat - 6*W1a^lam1*lam2 - 6*W2a^lam2*lam2 + 6*W1a^lam1*W2a^lam2*lam2 - 6*W1a^lam1*lam2^2*m2ahat)*cova^2*s1ahat*s2ahat^2 + (8*W1a^lam1*lam2^2*m1ahat - 8*lam2^2*m1ahat)*cova*s1ahat*s2ahat^4 + (2*lam2 + 2*lam2^2*m2ahat - 2*W1a^lam1*lam2 - 2*W2a^lam2*lam2 + 2*W1a^lam1*W2a^lam2*lam2 - 2*W1a^lam1*lam2^2*m2ahat)*s1ahat^3*s2ahat^4) - cova*s1ahat*s2ahat^4*(4*W1a^(2*lam1)*lam2^2 - 8*W1a^lam1*lam2^2 + 4*lam2^2))/(lam1^2*lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3))
    h2_10=sum(0)
    h2_11=sum((2*s1ahat*s2ahat^2*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(cova*lam1 - lam2*s2ahat^2 - W2a^lam2*cova*lam1 + W1a^lam1*lam2*s2ahat^2 + cova*lam1*lam2*m2ahat - lam1*lam2*m1ahat*s2ahat^2))/(lam1^3*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h2_12=sum(-(2*cova*s1ahat*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1)*(cova*lam1 - lam2*s2ahat^2 - W2a^lam2*cova*lam1 + W1a^lam1*lam2*s2ahat^2 + cova*lam1*lam2*m2ahat - lam1*lam2*m1ahat*s2ahat^2))/(lam1*lam2^3*(s1ahat^2*s2ahat^2 - cova^2)^2))


    h3_1=sum(0)
    h3_2=sum(0)
    h3_3=sum(-2/(2*s1bhat^2 - (2*covb^2)/s2bhat^2))
    h3_4=sum(-(2*s1bhat*s2bhat^2*(covb*lam1 - lam2*s2bhat^2 - W2b^lam2*covb*lam1 + W1b^lam1*lam2*s2bhat^2 + covb*lam1*lam2*m2bhat - lam1*lam2*m1bhat*s2bhat^2))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h3_5=sum(0)
    h3_6=sum(0)
    h3_7=sum(covb/(- covb^2 + s1bhat^2*s2bhat^2))
    h3_8=sum((s2bhat*(2*covb^2*(lam2 - W1b^lam1*lam2 + lam1*lam2*m1bhat) - 2*covb*s1bhat^2*(lam1 - W2b^lam2*lam1 + lam1*lam2*m2bhat)))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h3_9=sum(0)
    h3_10=sum((covb^2*(lam2*m2bhat - W2b^lam2 + 1) + s2bhat^2*((lam2*m2bhat - W2b^lam2 + 1)*s1bhat^2 - 2*covb*lam2*m1bhat))/(lam2*(s1bhat^2*s2bhat^2 - covb^2)^2) - (covb*s2bhat^2*(2*lam2 - 2*W1b^lam1*lam2))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h3_11=sum((2*s2bhat^2*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1))/(lam1^2*(- 2*covb^2 + 2*s1bhat^2*s2bhat^2)))
    h3_12=sum(-(covb*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1))/(lam2^2*(- covb^2 + s1bhat^2*s2bhat^2)))



    h4_1=sum(0)
    h4_2=sum(0)
    h4_3=sum(-(2*s1bhat*s2bhat^2*(covb*lam1 - lam2*s2bhat^2 - W2b^lam2*covb*lam1 + W1b^lam1*lam2*s2bhat^2 + covb*lam1*lam2*m2bhat - lam1*lam2*m1bhat*s2bhat^2))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h4_4=sum(- (3*s1bhat^2*s2bhat^6 + covb^4*(lam1^2*m2bhat^2 + lam1^2*s2bhat^2) - covb^3*(2*lam1*m2bhat*s2bhat^2 - 2*W1b^lam1*lam1*m2bhat*s2bhat^2 + 2*lam1^2*m1bhat*m2bhat*s2bhat^2) - covb*(6*lam1*m2bhat*s1bhat^2*s2bhat^4 - 6*W1b^lam1*lam1*m2bhat*s1bhat^2*s2bhat^4 + 6*lam1^2*m1bhat*m2bhat*s1bhat^2*s2bhat^4) + covb^2*(W1b^(2*lam1)*s2bhat^4 - 2*W1b^lam1*s2bhat^4 + s2bhat^4 + 2*lam1*m1bhat*s2bhat^4 + lam1^2*m1bhat^2*s2bhat^4 - 2*W1b^lam1*lam1*m1bhat*s2bhat^4 + 3*lam1^2*m2bhat^2*s1bhat^2*s2bhat^2) - 6*W1b^lam1*s1bhat^2*s2bhat^6 + 3*W1b^(2*lam1)*s1bhat^2*s2bhat^6 - lam1^2*s1bhat^4*s2bhat^6 + 3*lam1^2*m1bhat^2*s1bhat^2*s2bhat^6 + 6*lam1*m1bhat*s1bhat^2*s2bhat^6 - 6*W1b^lam1*lam1*m1bhat*s1bhat^2*s2bhat^6)/(lam1^2*(s1bhat^2*s2bhat^2 - covb^2)^3) - (covb^4*(W2b^(2*lam2)*lam1^2 - 2*W2b^lam2*lam1^2 + lam1^2) + covb^2*(3*lam1^2*s1bhat^2*s2bhat^2 + 3*W2b^(2*lam2)*lam1^2*s1bhat^2*s2bhat^2 - 6*W2b^lam2*lam1^2*s1bhat^2*s2bhat^2) - lam2*((2*W2b^lam2*lam1^2*m2bhat - 2*lam1^2*m2bhat)*covb^4 + (2*lam1*s2bhat^2 + 2*lam1^2*m1bhat*s2bhat^2 - 2*W1b^lam1*lam1*s2bhat^2 - 2*W2b^lam2*lam1*s2bhat^2 - 2*W2b^lam2*lam1^2*m1bhat*s2bhat^2 + 2*W1b^lam1*W2b^lam2*lam1*s2bhat^2)*covb^3 + (6*W2b^lam2*lam1^2*m2bhat*s1bhat^2*s2bhat^2 - 6*lam1^2*m2bhat*s1bhat^2*s2bhat^2)*covb^2 + (6*lam1*s1bhat^2*s2bhat^4 + 6*lam1^2*m1bhat*s1bhat^2*s2bhat^4 - 6*W1b^lam1*lam1*s1bhat^2*s2bhat^4 - 6*W2b^lam2*lam1*s1bhat^2*s2bhat^4 - 6*W2b^lam2*lam1^2*m1bhat*s1bhat^2*s2bhat^4 + 6*W1b^lam1*W2b^lam2*lam1*s1bhat^2*s2bhat^4)*covb))/(lam1^2*lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3))
    h4_5=sum(0)
    h4_6=sum(0)
    h4_7=sum((s1bhat*(2*covb^2*(lam1 - W2b^lam2*lam1 + lam1*lam2*m2bhat) - 2*covb*s2bhat^2*(lam2 - W1b^lam1*lam2 + lam1*lam2*m1bhat)))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h4_8=sum((lam2*((2*covb*s2bhat^3*(2*lam1 + 2*lam1^2*m1bhat - 2*W1b^lam1*lam1 - 2*W2b^lam2*lam1 + 2*W1b^lam1*W2b^lam2*lam1 - 2*W2b^lam2*lam1^2*m1bhat) - 2*covb*s2bhat*(4*covb*lam1^2*m2bhat - 4*W2b^lam2*covb*lam1^2*m2bhat))*s1bhat^3 + 2*covb*s2bhat*(2*covb^2*lam1 + 2*covb^2*lam1^2*m1bhat - 2*W1b^lam1*covb^2*lam1 - 2*W2b^lam2*covb^2*lam1 - 2*W2b^lam2*covb^2*lam1^2*m1bhat + 2*W1b^lam1*W2b^lam2*covb^2*lam1)*s1bhat) - 2*covb*s1bhat^3*s2bhat*(2*covb*lam1^2 + 2*W2b^(2*lam2)*covb*lam1^2 - 4*W2b^lam2*covb*lam1^2))/(lam1^2*lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3) - ((4*covb^2*lam1^2*m2bhat^2*s2bhat - 2*covb*s2bhat^3*(2*lam1*m2bhat + covb*lam1^2 + 2*lam1^2*m1bhat*m2bhat - 2*W1b^lam1*lam1*m2bhat))*s1bhat^3 + (2*covb*(2*covb + 2*W1b^(2*lam1)*covb - 4*W1b^lam1*covb + 2*covb*lam1^2*m1bhat^2 + 4*covb*lam1*m1bhat - 4*W1b^lam1*covb*lam1*m1bhat)*s2bhat^3 + 2*covb*(covb^3*lam1^2 - 2*covb^2*lam1*m2bhat + 2*W1b^lam1*covb^2*lam1*m2bhat - 2*covb^2*lam1^2*m1bhat*m2bhat)*s2bhat)*s1bhat)/(lam1^2*(s1bhat^2*s2bhat^2 - covb^2)^3))
    h4_9=sum(0)
    h4_10=sum(- (covb*(s2bhat^4*(- 4*lam2^2*m1bhat^2*s1bhat + 2*lam2^2*s1bhat^3) - s1bhat^3*s2bhat^2*(4*lam2*m2bhat - 4*W2b^lam2 + 2*W2b^(2*lam2) + 2*lam2^2*m2bhat^2 - 4*W2b^lam2*lam2*m2bhat + 2)) - covb^3*(s1bhat*(4*lam2*m2bhat - 4*W2b^lam2 + 2*W2b^(2*lam2) + 2*lam2^2*m2bhat^2 - 4*W2b^lam2*lam2*m2bhat + 2) + 2*lam2^2*s1bhat*s2bhat^2) + s1bhat^3*s2bhat^4*(2*lam2*m1bhat + 2*lam2^2*m1bhat*m2bhat - 2*W2b^lam2*lam2*m1bhat) + covb^2*s1bhat*s2bhat^2*(6*lam2*m1bhat + 6*lam2^2*m1bhat*m2bhat - 6*W2b^lam2*lam2*m1bhat))/(lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3) - (lam1*((6*lam2 + 6*lam2^2*m2bhat - 6*W1b^lam1*lam2 - 6*W2b^lam2*lam2 + 6*W1b^lam1*W2b^lam2*lam2 - 6*W1b^lam1*lam2^2*m2bhat)*covb^2*s1bhat*s2bhat^2 + (8*W1b^lam1*lam2^2*m1bhat - 8*lam2^2*m1bhat)*covb*s1bhat*s2bhat^4 + (2*lam2 + 2*lam2^2*m2bhat - 2*W1b^lam1*lam2 - 2*W2b^lam2*lam2 + 2*W1b^lam1*W2b^lam2*lam2 - 2*W1b^lam1*lam2^2*m2bhat)*s1bhat^3*s2bhat^4) - covb*s1bhat*s2bhat^4*(4*W1b^(2*lam1)*lam2^2 - 8*W1b^lam1*lam2^2 + 4*lam2^2))/(lam1^2*lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3))
    h4_11=sum((2*s1bhat*s2bhat^2*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(covb*lam1 - lam2*s2bhat^2 - W2b^lam2*covb*lam1 + W1b^lam1*lam2*s2bhat^2 + covb*lam1*lam2*m2bhat - lam1*lam2*m1bhat*s2bhat^2))/(lam1^3*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h4_12=sum(-(2*covb*s1bhat*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1)*(covb*lam1 - lam2*s2bhat^2 - W2b^lam2*covb*lam1 + W1b^lam1*lam2*s2bhat^2 + covb*lam1*lam2*m2bhat - lam1*lam2*m1bhat*s2bhat^2))/(lam1*lam2^3*(s1bhat^2*s2bhat^2 - covb^2)^2))


    h5_1=sum(cova/(- cova^2 + s1ahat^2*s2ahat^2))
    h5_2=sum((s1ahat*(2*cova^2*(lam1 - W2a^lam2*lam1 + lam1*lam2*m2ahat) - 2*cova*s2ahat^2*(lam2 - W1a^lam1*lam2 + lam1*lam2*m1ahat)))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h5_3=0
    h5_4=0
    h5_5=sum(-2/(2*s2ahat^2 - (2*cova^2)/s1ahat^2))
    h5_6=sum(-(2*s1ahat^2*s2ahat*(cova*lam2 - lam1*s1ahat^2 - W1a^lam1*cova*lam2 + W2a^lam2*lam1*s1ahat^2 + cova*lam1*lam2*m1ahat - lam1*lam2*m2ahat*s1ahat^2))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h5_7=0
    h5_8=0
    h5_9=sum((cova^2*(lam1*m1ahat - W1a^lam1 + 1) + s1ahat^2*((lam1*m1ahat - W1a^lam1 + 1)*s2ahat^2 - 2*cova*lam1*m2ahat))/(lam1*(s1ahat^2*s2ahat^2 - cova^2)^2) - (cova*s1ahat^2*(2*lam1 - 2*W2a^lam2*lam1))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h5_10=0
    h5_11=sum(-(cova*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1))/(lam1^2*(- cova^2 + s1ahat^2*s2ahat^2)))
    h5_12=sum((2*s1ahat^2*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1))/(lam2^2*(- 2*cova^2 + 2*s1ahat^2*s2ahat^2)))



    h6_1=sum((s2ahat*(2*cova^2*(lam2 - W1a^lam1*lam2 + lam1*lam2*m1ahat) - 2*cova*s1ahat^2*(lam1 - W2a^lam2*lam1 + lam1*lam2*m2ahat)))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h6_2=sum((lam2*((2*cova*s2ahat^3*(2*lam1 + 2*lam1^2*m1ahat - 2*W1a^lam1*lam1 - 2*W2a^lam2*lam1 + 2*W1a^lam1*W2a^lam2*lam1 - 2*W2a^lam2*lam1^2*m1ahat) - 2*cova*s2ahat*(4*cova*lam1^2*m2ahat - 4*W2a^lam2*cova*lam1^2*m2ahat))*s1ahat^3 + 2*cova*s2ahat*(2*cova^2*lam1 + 2*cova^2*lam1^2*m1ahat - 2*W1a^lam1*cova^2*lam1 - 2*W2a^lam2*cova^2*lam1 - 2*W2a^lam2*cova^2*lam1^2*m1ahat + 2*W1a^lam1*W2a^lam2*cova^2*lam1)*s1ahat) - 2*cova*s1ahat^3*s2ahat*(2*cova*lam1^2 + 2*W2a^(2*lam2)*cova*lam1^2 - 4*W2a^lam2*cova*lam1^2))/(lam1^2*lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3) - ((4*cova^2*lam1^2*m2ahat^2*s2ahat - 2*cova*s2ahat^3*(2*lam1*m2ahat + cova*lam1^2 + 2*lam1^2*m1ahat*m2ahat - 2*W1a^lam1*lam1*m2ahat))*s1ahat^3 + (2*cova*(2*cova + 2*W1a^(2*lam1)*cova - 4*W1a^lam1*cova + 2*cova*lam1^2*m1ahat^2 + 4*cova*lam1*m1ahat - 4*W1a^lam1*cova*lam1*m1ahat)*s2ahat^3 + 2*cova*(cova^3*lam1^2 - 2*cova^2*lam1*m2ahat + 2*W1a^lam1*cova^2*lam1*m2ahat - 2*cova^2*lam1^2*m1ahat*m2ahat)*s2ahat)*s1ahat)/(lam1^2*(s1ahat^2*s2ahat^2 - cova^2)^3))
    h6_3=0
    h6_4=0
    h6_5=sum(-(2*s1ahat^2*s2ahat*(cova*lam2 - lam1*s1ahat^2 - W1a^lam1*cova*lam2 + W2a^lam2*lam1*s1ahat^2 + cova*lam1*lam2*m1ahat - lam1*lam2*m2ahat*s1ahat^2))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h6_6=sum(- (3*s1ahat^6*s2ahat^2 + cova^4*(lam2^2*m1ahat^2 + lam2^2*s1ahat^2) - cova^3*(2*lam2*m1ahat*s1ahat^2 - 2*W2a^lam2*lam2*m1ahat*s1ahat^2 + 2*lam2^2*m1ahat*m2ahat*s1ahat^2) - cova*(6*lam2*m1ahat*s1ahat^4*s2ahat^2 - 6*W2a^lam2*lam2*m1ahat*s1ahat^4*s2ahat^2 + 6*lam2^2*m1ahat*m2ahat*s1ahat^4*s2ahat^2) + cova^2*(W2a^(2*lam2)*s1ahat^4 - 2*W2a^lam2*s1ahat^4 + s1ahat^4 + 2*lam2*m2ahat*s1ahat^4 + lam2^2*m2ahat^2*s1ahat^4 - 2*W2a^lam2*lam2*m2ahat*s1ahat^4 + 3*lam2^2*m1ahat^2*s1ahat^2*s2ahat^2) - 6*W2a^lam2*s1ahat^6*s2ahat^2 + 3*W2a^(2*lam2)*s1ahat^6*s2ahat^2 - lam2^2*s1ahat^6*s2ahat^4 + 3*lam2^2*m2ahat^2*s1ahat^6*s2ahat^2 + 6*lam2*m2ahat*s1ahat^6*s2ahat^2 - 6*W2a^lam2*lam2*m2ahat*s1ahat^6*s2ahat^2)/(lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3) - (cova^4*(W1a^(2*lam1)*lam2^2 - 2*W1a^lam1*lam2^2 + lam2^2) + cova^2*(3*lam2^2*s1ahat^2*s2ahat^2 + 3*W1a^(2*lam1)*lam2^2*s1ahat^2*s2ahat^2 - 6*W1a^lam1*lam2^2*s1ahat^2*s2ahat^2) - lam1*((2*W1a^lam1*lam2^2*m1ahat - 2*lam2^2*m1ahat)*cova^4 + (2*lam2*s1ahat^2 + 2*lam2^2*m2ahat*s1ahat^2 - 2*W1a^lam1*lam2*s1ahat^2 - 2*W2a^lam2*lam2*s1ahat^2 - 2*W1a^lam1*lam2^2*m2ahat*s1ahat^2 + 2*W1a^lam1*W2a^lam2*lam2*s1ahat^2)*cova^3 + (6*W1a^lam1*lam2^2*m1ahat*s1ahat^2*s2ahat^2 - 6*lam2^2*m1ahat*s1ahat^2*s2ahat^2)*cova^2 + (6*lam2*s1ahat^4*s2ahat^2 + 6*lam2^2*m2ahat*s1ahat^4*s2ahat^2 - 6*W1a^lam1*lam2*s1ahat^4*s2ahat^2 - 6*W2a^lam2*lam2*s1ahat^4*s2ahat^2 - 6*W1a^lam1*lam2^2*m2ahat*s1ahat^4*s2ahat^2 + 6*W1a^lam1*W2a^lam2*lam2*s1ahat^4*s2ahat^2)*cova))/(lam1^2*lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3))
    h6_7=0
    h6_8=0
    h6_9=sum(- (cova*(s1ahat^4*(- 4*lam1^2*m2ahat^2*s2ahat + 2*lam1^2*s2ahat^3) - s1ahat^2*s2ahat^3*(4*lam1*m1ahat - 4*W1a^lam1 + 2*W1a^(2*lam1) + 2*lam1^2*m1ahat^2 - 4*W1a^lam1*lam1*m1ahat + 2)) - cova^3*(s2ahat*(4*lam1*m1ahat - 4*W1a^lam1 + 2*W1a^(2*lam1) + 2*lam1^2*m1ahat^2 - 4*W1a^lam1*lam1*m1ahat + 2) + 2*lam1^2*s1ahat^2*s2ahat) + s1ahat^4*s2ahat^3*(2*lam1*m2ahat + 2*lam1^2*m1ahat*m2ahat - 2*W1a^lam1*lam1*m2ahat) + cova^2*s1ahat^2*s2ahat*(6*lam1*m2ahat + 6*lam1^2*m1ahat*m2ahat - 6*W1a^lam1*lam1*m2ahat))/(lam1^2*(s1ahat^2*s2ahat^2 - cova^2)^3) - (lam2*((6*lam1 + 6*lam1^2*m1ahat - 6*W1a^lam1*lam1 - 6*W2a^lam2*lam1 + 6*W1a^lam1*W2a^lam2*lam1 - 6*W2a^lam2*lam1^2*m1ahat)*cova^2*s1ahat^2*s2ahat + (8*W2a^lam2*lam1^2*m2ahat - 8*lam1^2*m2ahat)*cova*s1ahat^4*s2ahat + (2*lam1 + 2*lam1^2*m1ahat - 2*W1a^lam1*lam1 - 2*W2a^lam2*lam1 + 2*W1a^lam1*W2a^lam2*lam1 - 2*W2a^lam2*lam1^2*m1ahat)*s1ahat^4*s2ahat^3) - cova*s1ahat^4*s2ahat*(4*W2a^(2*lam2)*lam1^2 - 8*W2a^lam2*lam1^2 + 4*lam1^2))/(lam1^2*lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3))
    h6_10=0
    h6_11=sum(-(2*cova*s2ahat*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(cova*lam2 - lam1*s1ahat^2 - W1a^lam1*cova*lam2 + W2a^lam2*lam1*s1ahat^2 + cova*lam1*lam2*m1ahat - lam1*lam2*m2ahat*s1ahat^2))/(lam1^3*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h6_12=sum((2*s1ahat^2*s2ahat*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1)*(cova*lam2 - lam1*s1ahat^2 - W1a^lam1*cova*lam2 + W2a^lam2*lam1*s1ahat^2 + cova*lam1*lam2*m1ahat - lam1*lam2*m2ahat*s1ahat^2))/(lam1*lam2^3*(s1ahat^2*s2ahat^2 - cova^2)^2))



    h7_1=0
    h7_2=0
    h7_3=sum(covb/(- covb^2 + s1bhat^2*s2bhat^2))
    h7_4=sum((s1bhat*(2*covb^2*(lam1 - W2b^lam2*lam1 + lam1*lam2*m2bhat) - 2*covb*s2bhat^2*(lam2 - W1b^lam1*lam2 + lam1*lam2*m1bhat)))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h7_5=0
    h7_6=0
    h7_7=sum(-2/(2*s2bhat^2 - (2*covb^2)/s1bhat^2))
    h7_8=sum(-(2*s1bhat^2*s2bhat*(covb*lam2 - lam1*s1bhat^2 - W1b^lam1*covb*lam2 + W2b^lam2*lam1*s1bhat^2 + covb*lam1*lam2*m1bhat - lam1*lam2*m2bhat*s1bhat^2))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h7_9=0
    h7_10=sum((covb^2*(lam1*m1bhat - W1b^lam1 + 1) + s1bhat^2*((lam1*m1bhat - W1b^lam1 + 1)*s2bhat^2 - 2*covb*lam1*m2bhat))/(lam1*(s1bhat^2*s2bhat^2 - covb^2)^2) - (covb*s1bhat^2*(2*lam1 - 2*W2b^lam2*lam1))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h7_11=sum(-(covb*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1))/(lam1^2*(- covb^2 + s1bhat^2*s2bhat^2)))
    h7_12=sum((2*s1bhat^2*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1))/(lam2^2*(- 2*covb^2 + 2*s1bhat^2*s2bhat^2)))


    h8_1=0
    h8_2=0
    h8_3=sum((s2bhat*(2*covb^2*(lam2 - W1b^lam1*lam2 + lam1*lam2*m1bhat) - 2*covb*s1bhat^2*(lam1 - W2b^lam2*lam1 + lam1*lam2*m2bhat)))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h8_4=sum((lam2*((2*covb*s2bhat^3*(2*lam1 + 2*lam1^2*m1bhat - 2*W1b^lam1*lam1 - 2*W2b^lam2*lam1 + 2*W1b^lam1*W2b^lam2*lam1 - 2*W2b^lam2*lam1^2*m1bhat) - 2*covb*s2bhat*(4*covb*lam1^2*m2bhat - 4*W2b^lam2*covb*lam1^2*m2bhat))*s1bhat^3 + 2*covb*s2bhat*(2*covb^2*lam1 + 2*covb^2*lam1^2*m1bhat - 2*W1b^lam1*covb^2*lam1 - 2*W2b^lam2*covb^2*lam1 - 2*W2b^lam2*covb^2*lam1^2*m1bhat + 2*W1b^lam1*W2b^lam2*covb^2*lam1)*s1bhat) - 2*covb*s1bhat^3*s2bhat*(2*covb*lam1^2 + 2*W2b^(2*lam2)*covb*lam1^2 - 4*W2b^lam2*covb*lam1^2))/(lam1^2*lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3) - ((4*covb^2*lam1^2*m2bhat^2*s2bhat - 2*covb*s2bhat^3*(2*lam1*m2bhat + covb*lam1^2 + 2*lam1^2*m1bhat*m2bhat - 2*W1b^lam1*lam1*m2bhat))*s1bhat^3 + (2*covb*(2*covb + 2*W1b^(2*lam1)*covb - 4*W1b^lam1*covb + 2*covb*lam1^2*m1bhat^2 + 4*covb*lam1*m1bhat - 4*W1b^lam1*covb*lam1*m1bhat)*s2bhat^3 + 2*covb*(covb^3*lam1^2 - 2*covb^2*lam1*m2bhat + 2*W1b^lam1*covb^2*lam1*m2bhat - 2*covb^2*lam1^2*m1bhat*m2bhat)*s2bhat)*s1bhat)/(lam1^2*(s1bhat^2*s2bhat^2 - covb^2)^3))
    h8_5=0
    h8_6=0
    h8_7=sum(-(2*s1bhat^2*s2bhat*(covb*lam2 - lam1*s1bhat^2 - W1b^lam1*covb*lam2 + W2b^lam2*lam1*s1bhat^2 + covb*lam1*lam2*m1bhat - lam1*lam2*m2bhat*s1bhat^2))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h8_8=sum(- (3*s1bhat^6*s2bhat^2 + covb^4*(lam2^2*m1bhat^2 + lam2^2*s1bhat^2) - covb^3*(2*lam2*m1bhat*s1bhat^2 - 2*W2b^lam2*lam2*m1bhat*s1bhat^2 + 2*lam2^2*m1bhat*m2bhat*s1bhat^2) - covb*(6*lam2*m1bhat*s1bhat^4*s2bhat^2 - 6*W2b^lam2*lam2*m1bhat*s1bhat^4*s2bhat^2 + 6*lam2^2*m1bhat*m2bhat*s1bhat^4*s2bhat^2) + covb^2*(W2b^(2*lam2)*s1bhat^4 - 2*W2b^lam2*s1bhat^4 + s1bhat^4 + 2*lam2*m2bhat*s1bhat^4 + lam2^2*m2bhat^2*s1bhat^4 - 2*W2b^lam2*lam2*m2bhat*s1bhat^4 + 3*lam2^2*m1bhat^2*s1bhat^2*s2bhat^2) - 6*W2b^lam2*s1bhat^6*s2bhat^2 + 3*W2b^(2*lam2)*s1bhat^6*s2bhat^2 - lam2^2*s1bhat^6*s2bhat^4 + 3*lam2^2*m2bhat^2*s1bhat^6*s2bhat^2 + 6*lam2*m2bhat*s1bhat^6*s2bhat^2 - 6*W2b^lam2*lam2*m2bhat*s1bhat^6*s2bhat^2)/(lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3) - (covb^4*(W1b^(2*lam1)*lam2^2 - 2*W1b^lam1*lam2^2 + lam2^2) + covb^2*(3*lam2^2*s1bhat^2*s2bhat^2 + 3*W1b^(2*lam1)*lam2^2*s1bhat^2*s2bhat^2 - 6*W1b^lam1*lam2^2*s1bhat^2*s2bhat^2) - lam1*((2*W1b^lam1*lam2^2*m1bhat - 2*lam2^2*m1bhat)*covb^4 + (2*lam2*s1bhat^2 + 2*lam2^2*m2bhat*s1bhat^2 - 2*W1b^lam1*lam2*s1bhat^2 - 2*W2b^lam2*lam2*s1bhat^2 - 2*W1b^lam1*lam2^2*m2bhat*s1bhat^2 + 2*W1b^lam1*W2b^lam2*lam2*s1bhat^2)*covb^3 + (6*W1b^lam1*lam2^2*m1bhat*s1bhat^2*s2bhat^2 - 6*lam2^2*m1bhat*s1bhat^2*s2bhat^2)*covb^2 + (6*lam2*s1bhat^4*s2bhat^2 + 6*lam2^2*m2bhat*s1bhat^4*s2bhat^2 - 6*W1b^lam1*lam2*s1bhat^4*s2bhat^2 - 6*W2b^lam2*lam2*s1bhat^4*s2bhat^2 - 6*W1b^lam1*lam2^2*m2bhat*s1bhat^4*s2bhat^2 + 6*W1b^lam1*W2b^lam2*lam2*s1bhat^4*s2bhat^2)*covb))/(lam1^2*lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3))
    h8_9=0
    h8_10=sum(- (covb*(s1bhat^4*(- 4*lam1^2*m2bhat^2*s2bhat + 2*lam1^2*s2bhat^3) - s1bhat^2*s2bhat^3*(4*lam1*m1bhat - 4*W1b^lam1 + 2*W1b^(2*lam1) + 2*lam1^2*m1bhat^2 - 4*W1b^lam1*lam1*m1bhat + 2)) - covb^3*(s2bhat*(4*lam1*m1bhat - 4*W1b^lam1 + 2*W1b^(2*lam1) + 2*lam1^2*m1bhat^2 - 4*W1b^lam1*lam1*m1bhat + 2) + 2*lam1^2*s1bhat^2*s2bhat) + s1bhat^4*s2bhat^3*(2*lam1*m2bhat + 2*lam1^2*m1bhat*m2bhat - 2*W1b^lam1*lam1*m2bhat) + covb^2*s1bhat^2*s2bhat*(6*lam1*m2bhat + 6*lam1^2*m1bhat*m2bhat - 6*W1b^lam1*lam1*m2bhat))/(lam1^2*(s1bhat^2*s2bhat^2 - covb^2)^3) - (lam2*((6*lam1 + 6*lam1^2*m1bhat - 6*W1b^lam1*lam1 - 6*W2b^lam2*lam1 + 6*W1b^lam1*W2b^lam2*lam1 - 6*W2b^lam2*lam1^2*m1bhat)*covb^2*s1bhat^2*s2bhat + (8*W2b^lam2*lam1^2*m2bhat - 8*lam1^2*m2bhat)*covb*s1bhat^4*s2bhat + (2*lam1 + 2*lam1^2*m1bhat - 2*W1b^lam1*lam1 - 2*W2b^lam2*lam1 + 2*W1b^lam1*W2b^lam2*lam1 - 2*W2b^lam2*lam1^2*m1bhat)*s1bhat^4*s2bhat^3) - covb*s1bhat^4*s2bhat*(4*W2b^(2*lam2)*lam1^2 - 8*W2b^lam2*lam1^2 + 4*lam1^2))/(lam1^2*lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3))
    h8_11=sum(-(2*covb*s2bhat*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(covb*lam2 - lam1*s1bhat^2 - W1b^lam1*covb*lam2 + W2b^lam2*lam1*s1bhat^2 + covb*lam1*lam2*m1bhat - lam1*lam2*m2bhat*s1bhat^2))/(lam1^3*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h8_12=sum((2*s1bhat^2*s2bhat*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1)*(covb*lam2 - lam1*s1bhat^2 - W1b^lam1*covb*lam2 + W2b^lam2*lam1*s1bhat^2 + covb*lam1*lam2*m1bhat - lam1*lam2*m2bhat*s1bhat^2))/(lam1*lam2^3*(s1bhat^2*s2bhat^2 - covb^2)^2))




    h9_1=sum((cova^2*(lam2*m2ahat - W2a^lam2 + 1) + s2ahat^2*((lam2*m2ahat - W2a^lam2 + 1)*s1ahat^2 - 2*cova*lam2*m1ahat))/(lam2*(s1ahat^2*s2ahat^2 - cova^2)^2) - (cova*s2ahat^2*(2*lam2 - 2*W1a^lam1*lam2))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h9_2=sum(- (cova*(s2ahat^4*(- 4*lam2^2*m1ahat^2*s1ahat + 2*lam2^2*s1ahat^3) - s1ahat^3*s2ahat^2*(4*lam2*m2ahat - 4*W2a^lam2 + 2*W2a^(2*lam2) + 2*lam2^2*m2ahat^2 - 4*W2a^lam2*lam2*m2ahat + 2)) - cova^3*(s1ahat*(4*lam2*m2ahat - 4*W2a^lam2 + 2*W2a^(2*lam2) + 2*lam2^2*m2ahat^2 - 4*W2a^lam2*lam2*m2ahat + 2) + 2*lam2^2*s1ahat*s2ahat^2) + s1ahat^3*s2ahat^4*(2*lam2*m1ahat + 2*lam2^2*m1ahat*m2ahat - 2*W2a^lam2*lam2*m1ahat) + cova^2*s1ahat*s2ahat^2*(6*lam2*m1ahat + 6*lam2^2*m1ahat*m2ahat - 6*W2a^lam2*lam2*m1ahat))/(lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3) - (lam1*((6*lam2 + 6*lam2^2*m2ahat - 6*W1a^lam1*lam2 - 6*W2a^lam2*lam2 + 6*W1a^lam1*W2a^lam2*lam2 - 6*W1a^lam1*lam2^2*m2ahat)*cova^2*s1ahat*s2ahat^2 + (8*W1a^lam1*lam2^2*m1ahat - 8*lam2^2*m1ahat)*cova*s1ahat*s2ahat^4 + (2*lam2 + 2*lam2^2*m2ahat - 2*W1a^lam1*lam2 - 2*W2a^lam2*lam2 + 2*W1a^lam1*W2a^lam2*lam2 - 2*W1a^lam1*lam2^2*m2ahat)*s1ahat^3*s2ahat^4) - cova*s1ahat*s2ahat^4*(4*W1a^(2*lam1)*lam2^2 - 8*W1a^lam1*lam2^2 + 4*lam2^2))/(lam1^2*lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3))
    h9_3=0
    h9_4=0
    h9_5=sum((cova^2*(lam1*m1ahat - W1a^lam1 + 1) + s1ahat^2*((lam1*m1ahat - W1a^lam1 + 1)*s2ahat^2 - 2*cova*lam1*m2ahat))/(lam1*(s1ahat^2*s2ahat^2 - cova^2)^2) - (cova*s1ahat^2*(2*lam1 - 2*W2a^lam2*lam1))/(lam1*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h9_6=sum(- (cova*(s1ahat^4*(- 4*lam1^2*m2ahat^2*s2ahat + 2*lam1^2*s2ahat^3) - s1ahat^2*s2ahat^3*(4*lam1*m1ahat - 4*W1a^lam1 + 2*W1a^(2*lam1) + 2*lam1^2*m1ahat^2 - 4*W1a^lam1*lam1*m1ahat + 2)) - cova^3*(s2ahat*(4*lam1*m1ahat - 4*W1a^lam1 + 2*W1a^(2*lam1) + 2*lam1^2*m1ahat^2 - 4*W1a^lam1*lam1*m1ahat + 2) + 2*lam1^2*s1ahat^2*s2ahat) + s1ahat^4*s2ahat^3*(2*lam1*m2ahat + 2*lam1^2*m1ahat*m2ahat - 2*W1a^lam1*lam1*m2ahat) + cova^2*s1ahat^2*s2ahat*(6*lam1*m2ahat + 6*lam1^2*m1ahat*m2ahat - 6*W1a^lam1*lam1*m2ahat))/(lam1^2*(s1ahat^2*s2ahat^2 - cova^2)^3) - (lam2*((6*lam1 + 6*lam1^2*m1ahat - 6*W1a^lam1*lam1 - 6*W2a^lam2*lam1 + 6*W1a^lam1*W2a^lam2*lam1 - 6*W2a^lam2*lam1^2*m1ahat)*cova^2*s1ahat^2*s2ahat + (8*W2a^lam2*lam1^2*m2ahat - 8*lam1^2*m2ahat)*cova*s1ahat^4*s2ahat + (2*lam1 + 2*lam1^2*m1ahat - 2*W1a^lam1*lam1 - 2*W2a^lam2*lam1 + 2*W1a^lam1*W2a^lam2*lam1 - 2*W2a^lam2*lam1^2*m1ahat)*s1ahat^4*s2ahat^3) - cova*s1ahat^4*s2ahat*(4*W2a^(2*lam2)*lam1^2 - 8*W2a^lam2*lam1^2 + 4*lam1^2))/(lam1^2*lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3))
    h9_7=0
    h9_8=0
    h9_9=sum(- (s1ahat^2*(cova^2*(6*lam2*m2ahat - 6*W2a^lam2 + 3*W2a^(2*lam2) + 3*lam2^2*m2ahat^2 - 6*W2a^lam2*lam2*m2ahat + 3) - cova*(6*lam2*m1ahat*s2ahat^2 - 6*W2a^lam2*lam2*m1ahat*s2ahat^2 + 6*lam2^2*m1ahat*m2ahat*s2ahat^2) + lam2^2*m1ahat^2*s2ahat^4) - cova^3*(2*lam2*m1ahat + 2*lam2^2*m1ahat*m2ahat - 2*W2a^lam2*lam2*m1ahat) + s1ahat^4*(W2a^(2*lam2)*s2ahat^2 - 2*W2a^lam2*s2ahat^2 + s2ahat^2 - lam2^2*s2ahat^4 + 2*lam2*m2ahat*s2ahat^2 + lam2^2*m2ahat^2*s2ahat^2 - 2*W2a^lam2*lam2*m2ahat*s2ahat^2) + cova^4*lam2^2 + 3*cova^2*lam2^2*m1ahat^2*s2ahat^2)/(lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3) - (cova^2*(3*lam2^2*s2ahat^2 - 6*W1a^lam1*lam2^2*s2ahat^2 + 3*W1a^(2*lam1)*lam2^2*s2ahat^2) + s1ahat^2*(lam2^2*s2ahat^4 - 2*W1a^lam1*lam2^2*s2ahat^4 + W1a^(2*lam1)*lam2^2*s2ahat^4) - lam1*(cova^3*(2*lam2 + 2*lam2^2*m2ahat - 2*W1a^lam1*lam2 - 2*W2a^lam2*lam2 + 2*W1a^lam1*W2a^lam2*lam2 - 2*W1a^lam1*lam2^2*m2ahat) - cova^2*(6*lam2^2*m1ahat*s2ahat^2 - 6*W1a^lam1*lam2^2*m1ahat*s2ahat^2) + s1ahat^2*(cova*(6*lam2*s2ahat^2 + 6*lam2^2*m2ahat*s2ahat^2 - 6*W1a^lam1*lam2*s2ahat^2 - 6*W2a^lam2*lam2*s2ahat^2 - 6*W1a^lam1*lam2^2*m2ahat*s2ahat^2 + 6*W1a^lam1*W2a^lam2*lam2*s2ahat^2) - 2*lam2^2*m1ahat*s2ahat^4 + 2*W1a^lam1*lam2^2*m1ahat*s2ahat^4)))/(lam1^2*lam2^2*(s1ahat^2*s2ahat^2 - cova^2)^3))
    h9_10=0
    h9_11=sum(-((W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(cova^2*lam1 - 2*cova*lam2*s2ahat^2 + lam1*s1ahat^2*s2ahat^2 - W2a^lam2*cova^2*lam1 + 2*W1a^lam1*cova*lam2*s2ahat^2 + cova^2*lam1*lam2*m2ahat - W2a^lam2*lam1*s1ahat^2*s2ahat^2 + lam1*lam2*m2ahat*s1ahat^2*s2ahat^2 - 2*cova*lam1*lam2*m1ahat*s2ahat^2))/(lam1^3*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h9_12=sum(-((W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1)*(cova^2*lam2 - 2*cova*lam1*s1ahat^2 + lam2*s1ahat^2*s2ahat^2 - W1a^lam1*cova^2*lam2 + 2*W2a^lam2*cova*lam1*s1ahat^2 + cova^2*lam1*lam2*m1ahat - W1a^lam1*lam2*s1ahat^2*s2ahat^2 + lam1*lam2*m1ahat*s1ahat^2*s2ahat^2 - 2*cova*lam1*lam2*m2ahat*s1ahat^2))/(lam1*lam2^3*(s1ahat^2*s2ahat^2 - cova^2)^2))







    h10_1=0
    h10_2=0
    h10_3=sum((covb^2*(lam2*m2bhat - W2b^lam2 + 1) + s2bhat^2*((lam2*m2bhat - W2b^lam2 + 1)*s1bhat^2 - 2*covb*lam2*m1bhat))/(lam2*(s1bhat^2*s2bhat^2 - covb^2)^2) - (covb*s2bhat^2*(2*lam2 - 2*W1b^lam1*lam2))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h10_4=sum(- (covb*(s2bhat^4*(- 4*lam2^2*m1bhat^2*s1bhat + 2*lam2^2*s1bhat^3) - s1bhat^3*s2bhat^2*(4*lam2*m2bhat - 4*W2b^lam2 + 2*W2b^(2*lam2) + 2*lam2^2*m2bhat^2 - 4*W2b^lam2*lam2*m2bhat + 2)) - covb^3*(s1bhat*(4*lam2*m2bhat - 4*W2b^lam2 + 2*W2b^(2*lam2) + 2*lam2^2*m2bhat^2 - 4*W2b^lam2*lam2*m2bhat + 2) + 2*lam2^2*s1bhat*s2bhat^2) + s1bhat^3*s2bhat^4*(2*lam2*m1bhat + 2*lam2^2*m1bhat*m2bhat - 2*W2b^lam2*lam2*m1bhat) + covb^2*s1bhat*s2bhat^2*(6*lam2*m1bhat + 6*lam2^2*m1bhat*m2bhat - 6*W2b^lam2*lam2*m1bhat))/(lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3) - (lam1*((6*lam2 + 6*lam2^2*m2bhat - 6*W1b^lam1*lam2 - 6*W2b^lam2*lam2 + 6*W1b^lam1*W2b^lam2*lam2 - 6*W1b^lam1*lam2^2*m2bhat)*covb^2*s1bhat*s2bhat^2 + (8*W1b^lam1*lam2^2*m1bhat - 8*lam2^2*m1bhat)*covb*s1bhat*s2bhat^4 + (2*lam2 + 2*lam2^2*m2bhat - 2*W1b^lam1*lam2 - 2*W2b^lam2*lam2 + 2*W1b^lam1*W2b^lam2*lam2 - 2*W1b^lam1*lam2^2*m2bhat)*s1bhat^3*s2bhat^4) - covb*s1bhat*s2bhat^4*(4*W1b^(2*lam1)*lam2^2 - 8*W1b^lam1*lam2^2 + 4*lam2^2))/(lam1^2*lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3))
    h10_5=0
    h10_6=0
    h10_7=sum((covb^2*(lam1*m1bhat - W1b^lam1 + 1) + s1bhat^2*((lam1*m1bhat - W1b^lam1 + 1)*s2bhat^2 - 2*covb*lam1*m2bhat))/(lam1*(s1bhat^2*s2bhat^2 - covb^2)^2) - (covb*s1bhat^2*(2*lam1 - 2*W2b^lam2*lam1))/(lam1*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h10_8=sum(- (covb*(s1bhat^4*(- 4*lam1^2*m2bhat^2*s2bhat + 2*lam1^2*s2bhat^3) - s1bhat^2*s2bhat^3*(4*lam1*m1bhat - 4*W1b^lam1 + 2*W1b^(2*lam1) + 2*lam1^2*m1bhat^2 - 4*W1b^lam1*lam1*m1bhat + 2)) - covb^3*(s2bhat*(4*lam1*m1bhat - 4*W1b^lam1 + 2*W1b^(2*lam1) + 2*lam1^2*m1bhat^2 - 4*W1b^lam1*lam1*m1bhat + 2) + 2*lam1^2*s1bhat^2*s2bhat) + s1bhat^4*s2bhat^3*(2*lam1*m2bhat + 2*lam1^2*m1bhat*m2bhat - 2*W1b^lam1*lam1*m2bhat) + covb^2*s1bhat^2*s2bhat*(6*lam1*m2bhat + 6*lam1^2*m1bhat*m2bhat - 6*W1b^lam1*lam1*m2bhat))/(lam1^2*(s1bhat^2*s2bhat^2 - covb^2)^3) - (lam2*((6*lam1 + 6*lam1^2*m1bhat - 6*W1b^lam1*lam1 - 6*W2b^lam2*lam1 + 6*W1b^lam1*W2b^lam2*lam1 - 6*W2b^lam2*lam1^2*m1bhat)*covb^2*s1bhat^2*s2bhat + (8*W2b^lam2*lam1^2*m2bhat - 8*lam1^2*m2bhat)*covb*s1bhat^4*s2bhat + (2*lam1 + 2*lam1^2*m1bhat - 2*W1b^lam1*lam1 - 2*W2b^lam2*lam1 + 2*W1b^lam1*W2b^lam2*lam1 - 2*W2b^lam2*lam1^2*m1bhat)*s1bhat^4*s2bhat^3) - covb*s1bhat^4*s2bhat*(4*W2b^(2*lam2)*lam1^2 - 8*W2b^lam2*lam1^2 + 4*lam1^2))/(lam1^2*lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3))
    h10_9=0
    h10_10=sum(- (s1bhat^2*(covb^2*(6*lam2*m2bhat - 6*W2b^lam2 + 3*W2b^(2*lam2) + 3*lam2^2*m2bhat^2 - 6*W2b^lam2*lam2*m2bhat + 3) - covb*(6*lam2*m1bhat*s2bhat^2 - 6*W2b^lam2*lam2*m1bhat*s2bhat^2 + 6*lam2^2*m1bhat*m2bhat*s2bhat^2) + lam2^2*m1bhat^2*s2bhat^4) - covb^3*(2*lam2*m1bhat + 2*lam2^2*m1bhat*m2bhat - 2*W2b^lam2*lam2*m1bhat) + s1bhat^4*(W2b^(2*lam2)*s2bhat^2 - 2*W2b^lam2*s2bhat^2 + s2bhat^2 - lam2^2*s2bhat^4 + 2*lam2*m2bhat*s2bhat^2 + lam2^2*m2bhat^2*s2bhat^2 - 2*W2b^lam2*lam2*m2bhat*s2bhat^2) + covb^4*lam2^2 + 3*covb^2*lam2^2*m1bhat^2*s2bhat^2)/(lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3) - (covb^2*(3*lam2^2*s2bhat^2 - 6*W1b^lam1*lam2^2*s2bhat^2 + 3*W1b^(2*lam1)*lam2^2*s2bhat^2) + s1bhat^2*(lam2^2*s2bhat^4 - 2*W1b^lam1*lam2^2*s2bhat^4 + W1b^(2*lam1)*lam2^2*s2bhat^4) - lam1*(covb^3*(2*lam2 + 2*lam2^2*m2bhat - 2*W1b^lam1*lam2 - 2*W2b^lam2*lam2 + 2*W1b^lam1*W2b^lam2*lam2 - 2*W1b^lam1*lam2^2*m2bhat) - covb^2*(6*lam2^2*m1bhat*s2bhat^2 - 6*W1b^lam1*lam2^2*m1bhat*s2bhat^2) + s1bhat^2*(covb*(6*lam2*s2bhat^2 + 6*lam2^2*m2bhat*s2bhat^2 - 6*W1b^lam1*lam2*s2bhat^2 - 6*W2b^lam2*lam2*s2bhat^2 - 6*W1b^lam1*lam2^2*m2bhat*s2bhat^2 + 6*W1b^lam1*W2b^lam2*lam2*s2bhat^2) - 2*lam2^2*m1bhat*s2bhat^4 + 2*W1b^lam1*lam2^2*m1bhat*s2bhat^4)))/(lam1^2*lam2^2*(s1bhat^2*s2bhat^2 - covb^2)^3))
    h10_11=sum(-((W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(covb^2*lam1 - 2*covb*lam2*s2bhat^2 + lam1*s1bhat^2*s2bhat^2 - W2b^lam2*covb^2*lam1 + 2*W1b^lam1*covb*lam2*s2bhat^2 + covb^2*lam1*lam2*m2bhat - W2b^lam2*lam1*s1bhat^2*s2bhat^2 + lam1*lam2*m2bhat*s1bhat^2*s2bhat^2 - 2*covb*lam1*lam2*m1bhat*s2bhat^2))/(lam1^3*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h10_12=sum(-((W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1)*(covb^2*lam2 - 2*covb*lam1*s1bhat^2 + lam2*s1bhat^2*s2bhat^2 - W1b^lam1*covb^2*lam2 + 2*W2b^lam2*covb*lam1*s1bhat^2 + covb^2*lam1*lam2*m1bhat - W1b^lam1*lam2*s1bhat^2*s2bhat^2 + lam1*lam2*m1bhat*s1bhat^2*s2bhat^2 - 2*covb*lam1*lam2*m2bhat*s1bhat^2))/(lam1*lam2^3*(s1bhat^2*s2bhat^2 - covb^2)^2))



    h11_1=sum((2*s2ahat^2*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1))/(lam1^2*(- 2*cova^2 + 2*s1ahat^2*s2ahat^2)))
    h11_2=sum((2*s1ahat*s2ahat^2*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(cova*lam1 - lam2*s2ahat^2 - W2a^lam2*cova*lam1 + W1a^lam1*lam2*s2ahat^2 + cova*lam1*lam2*m2ahat - lam1*lam2*m1ahat*s2ahat^2))/(lam1^3*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h11_3=sum((2*s2bhat^2*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1))/(lam1^2*(- 2*covb^2 + 2*s1bhat^2*s2bhat^2)))
    h11_4=sum((2*s1bhat*s2bhat^2*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(covb*lam1 - lam2*s2bhat^2 - W2b^lam2*covb*lam1 + W1b^lam1*lam2*s2bhat^2 + covb*lam1*lam2*m2bhat - lam1*lam2*m1bhat*s2bhat^2))/(lam1^3*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h11_5=sum(-(cova*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1))/(lam1^2*(- cova^2 + s1ahat^2*s2ahat^2)))
    h11_6=sum(-(2*cova*s2ahat*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(cova*lam2 - lam1*s1ahat^2 - W1a^lam1*cova*lam2 + W2a^lam2*lam1*s1ahat^2 + cova*lam1*lam2*m1ahat - lam1*lam2*m2ahat*s1ahat^2))/(lam1^3*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h11_7=sum(-(covb*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1))/(lam1^2*(- covb^2 + s1bhat^2*s2bhat^2)))
    h11_8=sum(-(2*covb*s2bhat*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(covb*lam2 - lam1*s1bhat^2 - W1b^lam1*covb*lam2 + W2b^lam2*lam1*s1bhat^2 + covb*lam1*lam2*m1bhat - lam1*lam2*m2bhat*s1bhat^2))/(lam1^3*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h11_9=sum(-((W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(cova^2*lam1 - 2*cova*lam2*s2ahat^2 + lam1*s1ahat^2*s2ahat^2 - W2a^lam2*cova^2*lam1 + 2*W1a^lam1*cova*lam2*s2ahat^2 + cova^2*lam1*lam2*m2ahat - W2a^lam2*lam1*s1ahat^2*s2ahat^2 + lam1*lam2*m2ahat*s1ahat^2*s2ahat^2 - 2*cova*lam1*lam2*m1ahat*s2ahat^2))/(lam1^3*lam2*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h11_10=sum(-((W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(covb^2*lam1 - 2*covb*lam2*s2bhat^2 + lam1*s1bhat^2*s2bhat^2 - W2b^lam2*covb^2*lam1 + 2*W1b^lam1*covb*lam2*s2bhat^2 + covb^2*lam1*lam2*m2bhat - W2b^lam2*lam1*s1bhat^2*s2bhat^2 + lam1*lam2*m2bhat*s1bhat^2*s2bhat^2 - 2*covb*lam1*lam2*m1bhat*s2bhat^2))/(lam1^3*lam2*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h11_11=sum(((2*((W1a^lam1-1)/lam1^2-(W1a^lam1*log(W1a))/lam1)^2)/s1ahat^2-(2*(m1ahat-(W1a^lam1-1)/lam1)*((2*W1a^lam1-2)/lam1^3-(2*W1a^lam1*log(W1a))/lam1^2+(W1a^lam1*log(W1a)^2)/lam1))/s1ahat^2+(2*cova*(m2ahat-(W2a^lam2-1)/lam2)*((2*W1a^lam1-2)/lam1^3-(2*W1a^lam1*log(W1a))/lam1^2+(W1a^lam1*log(W1a)^2)/lam1))/(s1ahat^2*s2ahat^2))/((2*cova^2)/(s1ahat^2*s2ahat^2)-2))+sum(((2*((W1b^lam1-1)/lam1^2-(W1b^lam1*log(W1b))/lam1)^2)/s1bhat^2-(2*(m1bhat-(W1b^lam1-1)/lam1)*((2*W1b^lam1-2)/lam1^3-(2*W1b^lam1*log(W1b))/lam1^2+(W1b^lam1*log(W1b)^2)/lam1))/s1bhat^2+(2*covb*(m2bhat-(W2b^lam2-1)/lam2)*((2*W1b^lam1-2)/lam1^3-(2*W1b^lam1*log(W1b))/lam1^2+(W1b^lam1*log(W1b)^2)/lam1))/(s1bhat^2*s2bhat^2))/((2*covb^2)/(s1bhat^2*s2bhat^2)-2))
    h11_12=sum((2*cova*(W1a^lam1*lam1*log(W1a)-W1a^lam1+1)*(W2a^lam2*lam2*log(W2a)-W2a^lam2+1))/(lam1^2*lam2^2*(-2*cova^2+2*s1ahat^2*s2ahat^2)))+sum((2*covb*(W1b^lam1*lam1*log(W1b)-W1b^lam1+1)*(W2b^lam2*lam2*log(W2b)-W2b^lam2+1))/(lam1^2*lam2^2*(-2*covb^2+2*s1bhat^2*s2bhat^2)))


    h12_1=sum(-(cova*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1))/(lam2^2*(- cova^2 + s1ahat^2*s2ahat^2)))
    h12_2=sum(-(2*cova*s1ahat*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1)*(cova*lam1 - lam2*s2ahat^2 - W2a^lam2*cova*lam1 + W1a^lam1*lam2*s2ahat^2 + cova*lam1*lam2*m2ahat - lam1*lam2*m1ahat*s2ahat^2))/(lam1*lam2^3*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h12_3=sum(-(covb*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1))/(lam2^2*(- covb^2 + s1bhat^2*s2bhat^2)))
    h12_4=sum(-(2*covb*s1bhat*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1)*(covb*lam1 - lam2*s2bhat^2 - W2b^lam2*covb*lam1 + W1b^lam1*lam2*s2bhat^2 + covb*lam1*lam2*m2bhat - lam1*lam2*m1bhat*s2bhat^2))/(lam1*lam2^3*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h12_5=sum((2*s1ahat^2*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1))/(lam2^2*(- 2*cova^2 + 2*s1ahat^2*s2ahat^2)))
    h12_6=sum((2*s1ahat^2*s2ahat*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1)*(cova*lam2 - lam1*s1ahat^2 - W1a^lam1*cova*lam2 + W2a^lam2*lam1*s1ahat^2 + cova*lam1*lam2*m1ahat - lam1*lam2*m2ahat*s1ahat^2))/(lam1*lam2^3*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h12_7=sum((2*s1bhat^2*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1))/(lam2^2*(- 2*covb^2 + 2*s1bhat^2*s2bhat^2)))
    h12_8=sum((2*s1bhat^2*s2bhat*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1)*(covb*lam2 - lam1*s1bhat^2 - W1b^lam1*covb*lam2 + W2b^lam2*lam1*s1bhat^2 + covb*lam1*lam2*m1bhat - lam1*lam2*m2bhat*s1bhat^2))/(lam1*lam2^3*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h12_9=sum(-((W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1)*(cova^2*lam2 - 2*cova*lam1*s1ahat^2 + lam2*s1ahat^2*s2ahat^2 - W1a^lam1*cova^2*lam2 + 2*W2a^lam2*cova*lam1*s1ahat^2 + cova^2*lam1*lam2*m1ahat - W1a^lam1*lam2*s1ahat^2*s2ahat^2 + lam1*lam2*m1ahat*s1ahat^2*s2ahat^2 - 2*cova*lam1*lam2*m2ahat*s1ahat^2))/(lam1*lam2^3*(s1ahat^2*s2ahat^2 - cova^2)^2))
    h12_10=sum(-((W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1)*(covb^2*lam2 - 2*covb*lam1*s1bhat^2 + lam2*s1bhat^2*s2bhat^2 - W1b^lam1*covb^2*lam2 + 2*W2b^lam2*covb*lam1*s1bhat^2 + covb^2*lam1*lam2*m1bhat - W1b^lam1*lam2*s1bhat^2*s2bhat^2 + lam1*lam2*m1bhat*s1bhat^2*s2bhat^2 - 2*covb*lam1*lam2*m2bhat*s1bhat^2))/(lam1*lam2^3*(s1bhat^2*s2bhat^2 - covb^2)^2))
    h12_11=sum((2*cova*(W1a^lam1*lam1*log(W1a)-W1a^lam1+1)*(W2a^lam2*lam2*log(W2a)-W2a^lam2+1))/(lam1^2*lam2^2*(-2*cova^2+2*s1ahat^2*s2ahat^2)))+sum((2*covb*(W1b^lam1*lam1*log(W1b)-W1b^lam1+1)*(W2b^lam2*lam2*log(W2b)-W2b^lam2+1))/(lam1^2*lam2^2*(-2*covb^2+2*s1bhat^2*s2bhat^2)))
    h12_12=sum(((2*((W2a^lam2-1)/lam2^2-(W2a^lam2*log(W2a))/lam2)^2)/s2ahat^2-(2*(m2ahat-(W2a^lam2-1)/lam2)*((2*W2a^lam2-2)/lam2^3-(2*W2a^lam2*log(W2a))/lam2^2+(W2a^lam2*log(W2a)^2)/lam2))/s2ahat^2+(2*cova*(m1ahat-(W1a^lam1-1)/lam1)*((2*W2a^lam2-2)/lam2^3-(2*W2a^lam2*log(W2a))/lam2^2+(W2a^lam2*log(W2a)^2)/lam2))/(s1ahat^2*s2ahat^2))/((2*cova^2)/(s1ahat^2*s2ahat^2)-2))+sum(((2*((W2b^lam2-1)/lam2^2-(W2b^lam2*log(W2b))/lam2)^2)/s2bhat^2-(2*(m2bhat-(W2b^lam2-1)/lam2)*((2*W2b^lam2-2)/lam2^3-(2*W2b^lam2*log(W2b))/lam2^2+(W2b^lam2*log(W2b)^2)/lam2))/s2bhat^2+(2*covb*(m1bhat-(W1b^lam1-1)/lam1)*((2*W2b^lam2-2)/lam2^3-(2*W2b^lam2*log(W2b))/lam2^2+(W2b^lam2*log(W2b)^2)/lam2))/(s1bhat^2*s2bhat^2))/((2*covb^2)/(s1bhat^2*s2bhat^2)-2))



    h1r=c(h1_1,h1_2,h1_3,h1_4,h1_5,h1_6,h1_7,h1_8,h1_9,h1_10,h1_11,h1_12)
    h2r=c(h2_1,h2_2,h2_3,h2_4,h2_5,h2_6,h2_7,h2_8,h2_9,h2_10,h2_11,h2_12)
    h3r=c(h3_1,h3_2,h3_3,h3_4,h3_5,h3_6,h3_7,h3_8,h3_9,h3_10,h3_11,h3_12)
    h4r=c(h4_1,h4_2,h4_3,h4_4,h4_5,h4_6,h4_7,h4_8,h4_9,h4_10,h4_11,h4_12)
    h5r=c(h5_1,h5_2,h5_3,h5_4,h5_5,h5_6,h5_7,h5_8,h5_9,h5_10,h5_11,h5_12)
    h6r=c(h6_1,h6_2,h6_3,h6_4,h6_5,h6_6,h6_7,h6_8,h6_9,h6_10,h6_11,h6_12)
    h7r=c(h7_1,h7_2,h7_3,h7_4,h7_5,h7_6,h7_7,h7_8,h7_9,h7_10,h7_11,h7_12)
    h8r=c(h8_1,h8_2,h8_3,h8_4,h8_5,h8_6,h8_7,h8_8,h8_9,h8_10,h8_11,h8_12)
    h9r=c(h9_1,h9_2,h9_3,h9_4,h9_5,h9_6,h9_7,h9_8,h9_9,h9_10,h9_11,h9_12)

    h10r=c(h10_1,h10_2,h10_3,h10_4,h10_5,h10_6,h10_7,h10_8,h10_9,h10_10,h10_11,h10_12)
    h11r=c(h11_1,h11_2,h11_3,h11_4,h11_5,h11_6,h11_7,h11_8,h11_9,h11_10,h11_11,h11_12)
    h12r=c(h12_1,h12_2,h12_3,h12_4,h12_5,h12_6,h12_7,h12_8,h12_9,h12_10,h12_11,h12_12)


    HCF=rbind(h1r,h2r,h3r,h4r,h5r,h6r,h7r,h8r,h9r,h10r,h11r,h12r)
    HCF[1,2]=0;HCF[1,6]=0;HCF[1,9]=0;
    HCF[2,1]=0;HCF[2,5]=0;
    HCF[3,4]=0;HCF[3,8]=0;HCF[3,10]=0;HCF[4,3]=0;
    HCF[4,7]=0;
    HCF[5,2]=0;HCF[5,6]=0;HCF[5,9]=0;
    HCF[6,1]=0;HCF[6,5]=0;
    HCF[7,4]=0;HCF[7,8]=0;HCF[7,10]=0;
    HCF[8,3]=0;HCF[8,7]=0;
    HCF[9,1]=0;HCF[9,5]=0;
    HCF[10,3]=0;
    HCF[10,7]=0;


    HCF[1,1]=na*HCF[1,1];
    HCF[1,5]=na*HCF[1,5];
    HCF[5,1]=na*HCF[5,1];

    HCF[3,3]=nb*HCF[3,3];
    HCF[3,7]=nb*HCF[3,7];
    HCF[7,3]=nb*HCF[7,3];

    HCF[7,7]=nb*HCF[7,7];
    HCF[5,5]=na*HCF[5,5];



    VV=inv(-HCF)

    #############################################################################
    #############################################################################


    #================ DIFFERENCES OF THE YOUDEN INDICES=======================
    m1a=m1ahat;
    s1a=s1ahat;
    m1b=m1bhat;
    s1b=s1bhat;

    m2a=m2ahat;
    s2a=s2ahat;
    m2b=m2bhat;
    s2b=s2bhat;


    a1=m1b-m1a;
    b1=s1b/s1a;
    c1=(m1a*(b1^2-1)-a1+b1*sqrt(a1^2+(b1^2-1)*s1a^2*log(b1^2)))/(b1^2-1);
    J1=pnorm((m1b-c1)/s1b)+pnorm((c1-m1a)/s1a)-1

    J1BCCC=J1

    a2=m2b-m2a;
    b2=s2b/s2a;
    c2=(m2a*(b2^2-1)-a2+b2*sqrt(a2^2+(b2^2-1)*s2a^2*log(b2^2)))/(b2^2-1);
    J2=pnorm((m2b-c2)/s2b)+pnorm((c2-m2a)/s2a)-1

    J2BCCC=J2



    dJ1_dm1a = (2^(1/2)*exp(-(m1a - (m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a)/(s1b^2/s1a^2 - 1))^2/(2*s1a^2))*((s1b^2/s1a^2 + (s1b*(2*m1a - 2*m1b))/(2*s1a*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2)))/(s1b^2/s1a^2 - 1) - 1))/(2*pi^(1/2)*s1a) - (2^(1/2)*exp(-(m1b - (m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a)/(s1b^2/s1a^2 - 1))^2/(2*s1b^2))*(s1b^2/s1a^2 + (s1b*(2*m1a - 2*m1b))/(2*s1a*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))))/(2*pi^(1/2)*s1b*(s1b^2/s1a^2 - 1));
    dJ1_ds1a = (2^(1/2)*exp(-(m1b - (m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a)/(s1b^2/s1a^2 - 1))^2/(2*s1b^2))*(((s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a^2 + (2*m1a*s1b^2)/s1a^3 + (s1b*(2*s1a*(s1b^2/s1a^2 - 1) + (2*s1b^2*log(s1b^2/s1a^2))/s1a - 2*s1a*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1)))/(2*s1a*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2)))/(s1b^2/s1a^2 - 1) - (2*s1b^2*(m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a))/(s1a^3*(s1b^2/s1a^2 - 1)^2)))/(2*pi^(1/2)*s1b) - (exp(-(m1a - (m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a)/(s1b^2/s1a^2 - 1))^2/(2*s1a^2))*((2^(1/2)*(((s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a^2 + (2*m1a*s1b^2)/s1a^3 + (s1b*(2*s1a*(s1b^2/s1a^2 - 1) + (2*s1b^2*log(s1b^2/s1a^2))/s1a - 2*s1a*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1)))/(2*s1a*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2)))/(s1b^2/s1a^2 - 1) - (2*s1b^2*(m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a))/(s1a^3*(s1b^2/s1a^2 - 1)^2)))/(2*s1a) - (2^(1/2)*(m1a - (m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a)/(s1b^2/s1a^2 - 1)))/(2*s1a^2)))/pi^(1/2);
    dJ1_dm1b = (2^(1/2)*exp(-(m1b - (m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a)/(s1b^2/s1a^2 - 1))^2/(2*s1b^2))*(((s1b*(2*m1a - 2*m1b))/(2*s1a*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2)) + 1)/(s1b^2/s1a^2 - 1) + 1))/(2*pi^(1/2)*s1b) - (2^(1/2)*exp(-(m1a - (m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a)/(s1b^2/s1a^2 - 1))^2/(2*s1a^2))*((s1b*(2*m1a - 2*m1b))/(2*s1a*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2)) + 1))/(2*pi^(1/2)*s1a*(s1b^2/s1a^2 - 1));
    dJ1_ds1b = (2^(1/2)*exp(-(m1a - (m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a)/(s1b^2/s1a^2 - 1))^2/(2*s1a^2))*((((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2)/s1a + (2*m1a*s1b)/s1a^2 + (s1b*(2*s1b*log(s1b^2/s1a^2) + (2*s1a^2*(s1b^2/s1a^2 - 1))/s1b))/(2*s1a*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2)))/(s1b^2/s1a^2 - 1) - (2*s1b*(m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a))/(s1a^2*(s1b^2/s1a^2 - 1)^2)))/(2*pi^(1/2)*s1a) - (exp(-(m1b - (m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a)/(s1b^2/s1a^2 - 1))^2/(2*s1b^2))*((2^(1/2)*((((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2)/s1a + (2*m1a*s1b)/s1a^2 + (s1b*(2*s1b*log(s1b^2/s1a^2) + (2*s1a^2*(s1b^2/s1a^2 - 1))/s1b))/(2*s1a*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2)))/(s1b^2/s1a^2 - 1) - (2*s1b*(m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a))/(s1a^2*(s1b^2/s1a^2 - 1)^2)))/(2*s1b) + (2^(1/2)*(m1b - (m1a - m1b + m1a*(s1b^2/s1a^2 - 1) + (s1b*((m1a - m1b)^2 + s1a^2*log(s1b^2/s1a^2)*(s1b^2/s1a^2 - 1))^(1/2))/s1a)/(s1b^2/s1a^2 - 1)))/(2*s1b^2)))/pi^(1/2);

    dJ1_dm2a =0;
    dJ1_ds2a =0;
    dJ1_dm2b =0;
    dJ1_ds2b =0;


    dJ2_dm1a =0;
    dJ2_ds1a =0;
    dJ2_dm1b =0;
    dJ2_ds1b =0;

    dJ2_dm2a = (2^(1/2)*exp(-(m2a - (m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a)/(s2b^2/s2a^2 - 1))^2/(2*s2a^2))*((s2b^2/s2a^2 + (s2b*(2*m2a - 2*m2b))/(2*s2a*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2)))/(s2b^2/s2a^2 - 1) - 1))/(2*pi^(1/2)*s2a) - (2^(1/2)*exp(-(m2b - (m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a)/(s2b^2/s2a^2 - 1))^2/(2*s2b^2))*(s2b^2/s2a^2 + (s2b*(2*m2a - 2*m2b))/(2*s2a*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))))/(2*pi^(1/2)*s2b*(s2b^2/s2a^2 - 1));
    dJ2_ds2a = (2^(1/2)*exp(-(m2b - (m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a)/(s2b^2/s2a^2 - 1))^2/(2*s2b^2))*(((s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a^2 + (2*m2a*s2b^2)/s2a^3 + (s2b*(2*s2a*(s2b^2/s2a^2 - 1) + (2*s2b^2*log(s2b^2/s2a^2))/s2a - 2*s2a*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1)))/(2*s2a*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2)))/(s2b^2/s2a^2 - 1) - (2*s2b^2*(m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a))/(s2a^3*(s2b^2/s2a^2 - 1)^2)))/(2*pi^(1/2)*s2b) - (exp(-(m2a - (m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a)/(s2b^2/s2a^2 - 1))^2/(2*s2a^2))*((2^(1/2)*(((s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a^2 + (2*m2a*s2b^2)/s2a^3 + (s2b*(2*s2a*(s2b^2/s2a^2 - 1) + (2*s2b^2*log(s2b^2/s2a^2))/s2a - 2*s2a*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1)))/(2*s2a*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2)))/(s2b^2/s2a^2 - 1) - (2*s2b^2*(m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a))/(s2a^3*(s2b^2/s2a^2 - 1)^2)))/(2*s2a) - (2^(1/2)*(m2a - (m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a)/(s2b^2/s2a^2 - 1)))/(2*s2a^2)))/pi^(1/2);
    dJ2_dm2b = (2^(1/2)*exp(-(m2b - (m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a)/(s2b^2/s2a^2 - 1))^2/(2*s2b^2))*(((s2b*(2*m2a - 2*m2b))/(2*s2a*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2)) + 1)/(s2b^2/s2a^2 - 1) + 1))/(2*pi^(1/2)*s2b) - (2^(1/2)*exp(-(m2a - (m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a)/(s2b^2/s2a^2 - 1))^2/(2*s2a^2))*((s2b*(2*m2a - 2*m2b))/(2*s2a*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2)) + 1))/(2*pi^(1/2)*s2a*(s2b^2/s2a^2 - 1));
    dJ2_ds2b = (2^(1/2)*exp(-(m2a - (m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a)/(s2b^2/s2a^2 - 1))^2/(2*s2a^2))*((((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2)/s2a + (2*m2a*s2b)/s2a^2 + (s2b*(2*s2b*log(s2b^2/s2a^2) + (2*s2a^2*(s2b^2/s2a^2 - 1))/s2b))/(2*s2a*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2)))/(s2b^2/s2a^2 - 1) - (2*s2b*(m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a))/(s2a^2*(s2b^2/s2a^2 - 1)^2)))/(2*pi^(1/2)*s2a) - (exp(-(m2b - (m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a)/(s2b^2/s2a^2 - 1))^2/(2*s2b^2))*((2^(1/2)*((((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2)/s2a + (2*m2a*s2b)/s2a^2 + (s2b*(2*s2b*log(s2b^2/s2a^2) + (2*s2a^2*(s2b^2/s2a^2 - 1))/s2b))/(2*s2a*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2)))/(s2b^2/s2a^2 - 1) - (2*s2b*(m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a))/(s2a^2*(s2b^2/s2a^2 - 1)^2)))/(2*s2b) + (2^(1/2)*(m2b - (m2a - m2b + m2a*(s2b^2/s2a^2 - 1) + (s2b*((m2a - m2b)^2 + s2a^2*log(s2b^2/s2a^2)*(s2b^2/s2a^2 - 1))^(1/2))/s2a)/(s2b^2/s2a^2 - 1)))/(2*s2b^2)))/pi^(1/2);



    VVpars=VV[1:8,1:8];
    VJ1=c(dJ1_dm1a, dJ1_ds1a, dJ1_dm1b, dJ1_ds1b, dJ1_dm2a, dJ1_ds2a, dJ1_dm2b, dJ1_ds2b)%*%VVpars%*%t(t(c(dJ1_dm1a, dJ1_ds1a, dJ1_dm1b, dJ1_ds1b, dJ1_dm2a, dJ1_ds2a, dJ1_dm2b, dJ1_ds2b)));
    VJ2=c(dJ2_dm1a, dJ2_ds1a, dJ2_dm1b, dJ2_ds1b, dJ2_dm2a, dJ2_ds2a, dJ2_dm2b, dJ2_ds2b)%*%VVpars%*%t(t(c(dJ2_dm1a, dJ2_ds1a, dJ2_dm1b, dJ2_ds1b, dJ2_dm2a, dJ2_ds2a, dJ2_dm2b, dJ2_ds2b)));
    COVJ12=c(dJ1_dm1a, dJ1_ds1a, dJ1_dm1b, dJ1_ds1b, dJ1_dm2a, dJ1_ds2a, dJ1_dm2b, dJ1_ds2b)%*%VVpars%*%t(t(c(dJ2_dm1a, dJ2_ds1a, dJ2_dm1b, dJ2_ds1b, dJ2_dm2a, dJ2_ds2a, dJ2_dm2b, dJ2_ds2b)));

    kk1=c(VJ1, COVJ12)
    kk2=c(COVJ12, VJ2)
    V12=rbind(kk1,kk2)

    Z=(J2-J1)/(sqrt(VJ1+VJ2-2*COVJ12));
    J1original=J1;
    J2original=J2;
    SE=sqrt(VJ1+VJ2-2*COVJ12);
    CIoriginal=c((J2-J1)-pmalpha*SE, (J2-J1)+pmalpha*SE);
    pval2tJ=2*pnorm(-abs(Z))




    #JT2=log(0.5*(J2+1)/(1-0.5*(J2+1)));
    #JT1=log(0.5*(J1+1)/(1-0.5*(J1+1)));
    #VJT1=4/(J1^2-1)^2*VJ1;
    #VJT2=4/(J2^2-1)^2*VJ2;
    #COVJT12=c(-2/(J1^2-1), 0)%*%V12%*%t(t(c(0, -2/(J2^2-1))))
    #Zstar=(JT2-JT1)/(sqrt(VJT1+VJT2-2*COVJT12))
    #SE=(sqrt(VJT1+VJT2-2*COVJT12))
    #CIZstar=c((JT2-JT1)-pmalpha*SE, (JT2-JT1)+pmalpha*SE)
    #CIZstar



    JT2=qnorm(J2);
    JT1=qnorm(J1);
    VJT1=(1/(dnorm(qnorm(J1))))^2*VJ1;
    VJT2=(1/(dnorm(qnorm(J2))))^2*VJ2;
    COVJT12=c((1/(dnorm(qnorm(J1))))*VJ1, 0)%*%V12%*%t(t(c(0, (1/(dnorm(qnorm(J2))))*VJ2)))
    Zstar=(JT2-JT1)/(sqrt(VJT1+VJT2-2*COVJT12))
    SE=(sqrt(VJT1+VJT2-2*COVJT12))
    CIZstar=c((JT2-JT1)-pmalpha*SE, (JT2-JT1)+pmalpha*SE)
    CIZstar


    pval2t=2*pnorm(-abs(Zstar))
    #CIoriginal=c(2*(exp(CIZstar[1])/(1+exp(CIZstar[1])))-1  , 2*(exp(CIZstar[2])/(1+exp(CIZstar[2])))-1)

    #====TWO ROC FUNCTIONS: ONE IS THE BOXCOX AND THE OTHER THE REFERENCE LINE================
    roc1<-function(t){
      na = length(W1alam)
      nb = length(W1blam)

      1-pnorm(qnorm(1-t,
                    mean=mean(W1alam),
                    sd=std(W1alam)*((na-1)/na)),
              mean=mean(W1blam),
              sd=std(W1blam)*((nb-1)/nb))
    }

    roc2<-function(t){
      na = length(W2alam)
      nb = length(W2blam)

      1-pnorm(qnorm(1-t,
                    mean=mean(W2alam),
                    sd=std(W2alam)*((na-1)/na)),
              mean=mean(W2blam),
              sd=std(W2blam)*((nb-1)/nb))
    }

    rocuseless<-function(t){
      1-pnorm(qnorm(1-t,mean=1,sd=1),mean=1,sd=1)
    }


    Sens1=1-pnorm(c1, mean=mean(W1blam), sd=sd(W1blam))
    Spec1=  pnorm(c1, mean=mean(W1alam), sd=sd(W1alam))

    Sens2=1-pnorm(c2, mean=mean(W2blam), sd=sd(W2blam))
    Spec2=  pnorm(c2, mean=mean(W2alam), sd=sd(W2alam))

    #================IF PLOTS ARE REQUESTED PLOT THE ROCS==
    if (plots=="on") {
      # x11()
      plot(linspace(0,1,1000),roc1(linspace(0,1,1000)),main="Box-Cox Based ROCs",xlab="FPR = 1 - Specificity",ylab="TPR = Sensitivity",type="l",col="red")
      lines(linspace(0,1,1000),roc2(linspace(0,1,1000)),main="Box-Cox Based ROCs",xlab="FPR = 1 - Specificity",ylab="TPR = Sensitivity",type="l",col="black")
      lines(linspace(0,1,1000),linspace(0,1,1000),type="l", lty=2)


      points(1-Spec1,Sens1, col = "red")
      points(1-Spec2,Sens2, col = "black")

      #line for the Youden:
      lines(c(1-Spec1,1-Spec1),c(rocuseless(1-Spec1),Sens1),col="blue")
      lines(c(1-Spec2,1-Spec2),c(rocuseless(1-Spec2),Sens2),col="green")

      legend("bottomright", legend=c(paste("ROC for Marker 1 with J =", round(J1original,4), " "),
                                     paste("ROC for Marker 2 with J =", round(J2original,4), " "),
                                     paste("P-value for the difference:",round(pval2t,4), " ")),

             col=c("red", "black", "white"), lty=c(1, 1, NA), pch = c(NA, NA, NA), cex=0.8)

    }






    # if (plots=="on"){
    #    points(1-Spec1,Sens1, col = "red")
    #    points(1-Spec2,Sens2, col = "black")
    #
    #    #line for the Youden:
    #    lines(c(1-Spec1,1-Spec1),c(rocuseless(1-Spec1),Sens1),col="blue")
    #    lines(c(1-Spec2,1-Spec2),c(rocuseless(1-Spec2),Sens2),col="green")
    #
    #  }

    #======================================================



    res <- matrix(c(J1original,J2original, pval2t,CIZstar[1],CIZstar[2]),ncol=5,byrow=TRUE)
    colnames(res) <- c("J1:","    J2:","    p-value (trans):","   CI trans (LL):","    CI trans (UL):")
    rownames(res) <- c("Estimates:")
    res <- as.table(res)
    res


    #return(list(AUCmarker1=AUC1original,AUCmarker2=AUC2original, pvalue_difference= pval2t, CI_difference= pnorm(CIdiff), rocbc1=roc1, rocbc2=roc2))
    #return(list(resultstable=res,J1=J1original,J2=J2original, pvalue_difference= pval2t, CI_difference= CIoriginal, rocbc1=roc1, rocbc2=roc2))
    return(list(resultstable=res,J1=J1original,J2=J2original, pvalue_probit_difference= pval2t, CI_probit_difference= CIZstar, pvalue_difference= pval2tJ, CI_difference= CIoriginal, roc1=roc1, roc2=roc2, transx1=W1alam, transy1=W1blam, transx2=W2alam, transy2=W2blam))

  }


}
