## ----eval=TRUE, error=FALSE, warning=FALSE, cache=FALSE, comment=FALSE, echo=FALSE, message=FALSE, results=FALSE----
# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

library(pracma)
library(clinfun)
library(splancs)
library(mvtnorm)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------

rocboxcox<-function(marker, D, alpha, plots, printProgress = FALSE){

  if (plots!="on"){plots="off"}
  if (length(marker) != length(D)) {
    stop("ERROR: The length of the 'marker' and 'D' inputs must be equal.")
  } else if (min(D) != 0 | max(D) != 1) {
    stop("ERROR: Controls must be assigned a value of 0; cases must be assigned a value of 1. Both controls and cases should be included in the dataset.")
  } else if (sum(is.na(marker)) > 0 | sum(is.na(D)) > 0) {
    stop("ERROR: Please remove all missing data before running this function.")
  } else if (alpha <= 0 | alpha >= 1) {
    stop("ERROR: The level of significance, alpha, should be set between 0 and 1. A common choice is 0.05.")
  } else {
  
  x=marker[D==0]
  y=marker[D==1]


  if (plots!="on"){plots="off"}
  n1=length(x)
  n2=length(y)
  xor=x;yor=y;
  Za=qnorm(1-alpha/2)

  likbox<-function(x,y,h){
    n=length(x);
    m=length(y);
    out=c();
    for (i in 1:length(h)){
      #print(i)
      if (h[i]==0){
        xh=log(x);
        yh=log(y);
      } else {
        xh=((x^h[i])-1)/h[i];
        yh=((y^h[i])-1)/h[i];
      }


      out[i]<-c( -n/2*log(sum((xh-sum(xh)/n)^2)/n)  -m/2*log(sum((yh-sum(yh)/m)^2)/m) +(h[i]-1)*(sum(log(x))+sum(log(y))))
    }
    return(out)
  }





  boxcoxleo<-function(x,y){
    init=1;
    logL<-function(h){
      -likbox(x,y,h)
    }
    #lam=fminsearch(logL,init)
    #lam=c(lam$optbase$xopt)
    lam<- optim(1,logL,gr=NULL,method="BFGS", control=list(maxit=10000))
    lam=c(lam$par)
    transx=((x^lam)-1)/lam
    transy=((y^lam)-1)/lam

    #test1=print(shapiro.test(transx))
    #test2=print(shapiro.test(transy))


    return(list(transformation.parameter=lam,transx=((x^lam)-1)/lam, transy=((y^lam)-1)/lam ))
  }


  boxcoxleo2<-function(x,y){
    init=1;
    logL<-function(h){
      -likbox(x,y,h)
    }
    lam<- optim(1,logL,gr=NULL,method="BFGS", control=list(maxit=10000))
    lam=c(lam$par)
    #lam=fminsearch(logL,init)
    #lam=c(lam$optbase$xopt)
    transx=((x^lam)-1)/lam
    transy=((y^lam)-1)/lam



    return(list(transformation.parameter=lam,transx=((x^lam)-1)/lam, transy=((y^lam)-1)/lam ))
  }


  cc=boxcoxleo(x,y)

  transx=cc$transx
  transy=cc$transy

  #x11()
  #qqnorm(x)
  #title(main="                                         for X")
  #x11()
  #qqnorm(y)
  #title(main="                                         for Y")

  roc<-function(t,x,y){
    1-pnorm(qnorm(1-t,mean=mean(transx),sd=std(transx)),mean=mean(transy),sd=std(transy))
  }

  rocuseless<-function(t){
    1-pnorm(qnorm(1-t,mean=1,sd=1),mean=1,sd=1)
  }

  if (plots=="on"){
    txt <- paste("Box-Cox Based ROC, alpha =", round(alpha,3) )  
    plot(linspace(0,1,1000),roc(linspace(0,1,1000),transx,transy),main=" ",xlab="FPR = 1 - Specificity",ylab="TPR = Sensitivity",type="l",col="red")
    lines(linspace(0,1,10),linspace(0,1,10),type="l", lty=2)
    title(main=txt)
    
  }
  rocfun=function(t){1-pnorm(qnorm(1-t,mean=mean(transx),sd=std(transx)),mean=mean(transy),sd=std(transy))}
  m1hat=mean(transx)
  m2hat=mean(transy)
  s1hat=std(transx)
  s2hat=std(transy)
  n1=length(x)
  n2=length(y)

  #NOTE THAT:  m1hat->m1
  #            m2hat->m2
  #            s1hat->s1hat
  #            s2hat->s2hat


  vardeltasp= (-((s2hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) + 1)/s1hat)^2*(s1hat^2/n1)+((m1hat + (m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2))/(s1hat^2 - s2hat^2))/s1hat^2 - ((s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - 2*m2hat*s1hat + (s1hat*s2hat*(2*s1hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s1hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) - (2*s1hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s1hat)^2*(s1hat^2/(2*n1-2))+((s1hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat*(s1hat^2 - s2hat^2)))^2*(s2hat^2/n2)+(-((2*m1hat*s2hat + s1hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - (s1hat*s2hat*(2*s2hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) + (2*s2hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s1hat)^2*(s2hat^2/(2*n2-2))

  #delta method for var of delta_se
  vardeltase= ((s2hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s2hat*(s1hat^2 - s2hat^2)))^2*(s1hat^2/n1)+(((s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - 2*m2hat*s1hat + (s1hat*s2hat*(2*s1hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s1hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) - (2*s1hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s2hat)^2*(s1hat^2/(2*n1-2))+(-((s1hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) - 1)/s2hat)^2*(s2hat^2/n2)+(((2*m1hat*s2hat + s1hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - (s1hat*s2hat*(2*s2hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) + (2*s2hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s2hat - (m2hat + (m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2))/(s1hat^2 - s2hat^2))/s2hat^2)^2*(s2hat^2/(2*n2-2))

  #delta method for delta_e, delta_p covariance
  covdedp= ((s2hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s2hat*(s1hat^2 - s2hat^2)))*(-((s2hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) + 1)/s1hat)*(s1hat^2/n1)+(((s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - 2*m2hat*s1hat + (s1hat*s2hat*(2*s1hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s1hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) - (2*s1hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s2hat)*((m1hat + (m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2))/(s1hat^2 - s2hat^2))/s1hat^2 - ((s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - 2*m2hat*s1hat + (s1hat*s2hat*(2*s1hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s1hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) - (2*s1hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s1hat)*(s1hat^2/(2*n1-2))+(-((s1hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) - 1)/s2hat)*((s1hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat*(s1hat^2 - s2hat^2)))*(s2hat^2/n2)+(((2*m1hat*s2hat + s1hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - (s1hat*s2hat*(2*s2hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) + (2*s2hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s2hat - (m2hat + (m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2))/(s1hat^2 - s2hat^2))/s2hat^2)*(-((2*m1hat*s2hat + s1hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - (s1hat*s2hat*(2*s2hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) + (2*s2hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s1hat)*(s2hat^2/(2*n2-2))

  cutoff= (s1hat^2*m2hat-s2hat^2*m1hat-s1hat*s2hat*sqrt((m1hat-m2hat)^2+(s1hat^2-s2hat^2)*log(s1hat^2/s2hat^2)))/(s1hat^2-s2hat^2)

  deltasp= (((s1hat^2*m2hat-s2hat^2*m1hat-s1hat*s2hat*sqrt((m1hat-m2hat)^2+(s1hat^2-s2hat^2)*log(s1hat^2./s2hat^2)))/(s1hat^2-s2hat^2))-m1hat)/(s1hat)
  deltase= (m2hat-((s1hat^2*m2hat-s2hat^2*m1hat-s1hat*s2hat*sqrt((m1hat-m2hat)^2+(s1hat^2-s2hat^2)*log(s1hat^2/s2hat^2)))/(s1hat^2-s2hat^2)))/s2hat

  Se=pnorm(deltase)
  Sp=pnorm(deltasp)

  Se
  Sp

  if (plots=="on"){
    points(1-Sp,Se, col = "red")
  }
  covdeltas = matrix( c(vardeltasp, covdedp, covdedp, vardeltase), nrow=2, ncol=2,  byrow = TRUE)        # fill m




  svdA=svd(inv(covdeltas))
  a=sqrt(qchisq(1-alpha,2))/sqrt(svdA$d[1])
  b=sqrt(qchisq(1-alpha,2))/sqrt(svdA$d[2])

  theta=seq(from=0,to=2*pi+1/40,by=1/40)
  state1=a*cos(theta)
  state2=b*sin(theta)
  states=rbind(state1,state2)

  X=svdA$v %*%states
  x1=X[1,]+deltasp
  x2=X[2,]+deltase
  px1=pnorm(x1);
  px2=pnorm(x2);
  yegg=px2
  xegg=1-px1
  if (plots=="on"){
    lines(xegg,yegg,col="green")
  }



  #----Area egg---------------
  bxegg=c(xegg,xegg[1])
  byegg=c(yegg,yegg[1])
  ell <- cbind(bxegg, byegg)
  areaegg=areapl(ell)
  #----End Area egg---------------



  margcisp=c(pnorm(deltasp-qnorm(1-alpha/2)*sqrt(vardeltasp)), pnorm(deltasp+qnorm(1-alpha/2)*sqrt(vardeltasp)))
  margcise=c(pnorm(deltase-qnorm(1-alpha/2)*sqrt(vardeltase)), pnorm(deltase+qnorm(1-alpha/2)*sqrt(vardeltase)))
  #-----Now deal with the rectangle-----
  cisp=c(deltasp-qnorm(1-alpha/4)*sqrt(vardeltasp), deltasp+qnorm(1-alpha/4)*sqrt(vardeltasp))
  cise=c(deltase-qnorm(1-alpha/4)*sqrt(vardeltase), deltase+qnorm(1-alpha/4)*sqrt(vardeltase))
  bsenslb=pnorm(cise[1])#exp(cibootsse[1])/(1+exp(cibootsse[1]))    #NOW trans back, going to construct a 95% rectangle for boots
  bsensub=pnorm(cise[2])#exp(cibootsse[2])/(1+exp(cibootsse[2]))    #NOW trans back, going to construct a 95% rectangle for boots
  bspeclb=pnorm(cisp[1])#exp(cibootssp[1])/(1+exp(cibootssp[1]))    #NOW trans back, going to construct a 95% rectangle for boots
  bspecub=pnorm(cisp[2])#exp(cibootssp[2])/(1+exp(cibootssp[2]))    #NOW trans back, going to construct a 95% rectangle for boots
  #-----End of the rectangle------------

  if (plots=="on"){
    #-----Start--plot---Rectangle--------------------
    lines(c(1-bspeclb,1-bspecub),c(bsenslb,bsenslb),col="black")
    lines(c(1-bspecub,1-bspecub),c(bsenslb,bsensub),col="black")
    lines(c(1-bspecub,1-bspeclb),c(bsensub,bsensub),col="black")
    lines(c(1-bspeclb,1-bspeclb),c(bsenslb,bsensub),col="black")
    #-----End--plot---Rectangle--------------------
  }

  brectx=c(1-bspeclb,1-bspecub, 1-bspecub,1-bspecub,1-bspecub,1-bspeclb,1-bspeclb,1-bspeclb)
  brecty=c(bsenslb,bsenslb ,bsenslb,bsensub, bsensub,bsensub ,bsenslb,bsensub)
  brectx=c(brectx, brectx[1]);
  brecty=c(brecty, brecty[1]);

  re <- cbind(brectx, brecty)
  arearect=areapl(re)

  #line for the Youden:
  if (plots=="on"){
    lines(c(1-Sp,1-Sp),c(rocuseless(1-Sp),Se),col="blue")
  }
  #================AUC================

  invAUC= ((m2hat-m1hat)/s2hat)/(sqrt(1+(s1hat/s2hat)^2));
  auc=pnorm(invAUC)


  #====================AUC1==============================
  dm1 =-1/(s2hat*sqrt(1+(s1hat/s2hat)^2));
  ds1 = (s1hat*(m1hat-m2hat))/(s2hat^3*(1+(s1hat/s2hat)^2)^(3/2))
  dm2 =1/(s2hat*sqrt(1+(s1hat/s2hat)^2));
  ds2 =((m1hat-m2hat))/(s2hat^2*(1+(s1hat/s2hat)^2)^(3/2))




  I=zeros(5,5);
  sh=sqrt(1/(length(transx))*sum((transx-mean(transx))^2))
  sd=sqrt(1/(length(transy))*sum((transy-mean(transy))^2))

  mh=m1hat;
  md=m2hat;
  n=n1;
  m=n2;
  xlam=transx;
  ylam=transy;
  lam=cc$transformation.parameter
  I[1,1]=n/sh^2;
  I[2,2]=-(n/sh^2-3/sh^4*sum((xlam-mh)^2));
  I[3,3]=m/sd^2;
  I[4,4]=-(m/sd^2-3/sd^4*sum((ylam-md)^2));

  kk=  sum(((mh - (x^lam - 1)/lam)*((2*(x^lam - 1))/lam^3 - (2*x^lam*log(x))/lam^2 + (x^lam*log(x)^2)/lam))/sh^2) + - sum(((y^lam - 1)/lam^2 - (y^lam*log(y))/lam)^2/sd^2) + - sum(((x^lam - 1)/lam^2 - (x^lam*log(x))/lam)^2/sh^2) + + sum(((md - (y^lam - 1)/lam)*((2*(y^lam - 1))/lam^3 - (2*y^lam*log(y))/lam^2 + (y^lam*log(y)^2)/lam))/sd^2);
  I[5,5]=-kk ;

  I[1,5]=sum(((2*(x^lam - 1))/lam^2 - (2*x^lam*log(x))/lam)/(2*sh^2));
  I[2,5]=-sum((2*((x^lam - 1)/lam^2 - (x^lam*log(x))/lam)*(mh - (x^lam - 1)/lam))/sh^3);
  I[3,5]=-sum(-((2*(y^lam - 1))/lam^2 - (2*y^lam*log(y))/lam)/(2*sd^2));


  I[4,5]=-sum((2*((y^lam - 1)/lam^2 - (y^lam*log(y))/lam)*(md - (y^lam - 1)/lam))/sd^3);



  I[5,1]=I[1,5]
  I[5,2]=I[2,5]
  I[5,3]=I[3,5]
  I[5,4]=I[4,5]

  S=inv(I)

  S=S[1:4,1:4]

  varauc=t(c(dm1,ds1,dm2,ds2))%*% S %*% t(t(c(dm1,ds1,dm2,ds2)))

  CIauc=c(invAUC-Za*sqrt(varauc), invAUC+Za*sqrt(varauc))
  CIauc=pnorm(CIauc)
  CIauc
  Zauc=invAUC/sqrt(varauc)
  pvalauc=2*pnorm(-abs(Zauc))





  #================= DeltaBClamJT================
  mx=m1hat;
  my=m2hat;
  sx=s1hat;
  sy=s2hat;
  c=cutoff;
  Jhat=pnorm((m2hat-cutoff)/s2hat)+pnorm((cutoff-m1hat)/s1hat)-1;
  radd=a^2+(b^2-1)*sx^2*log(b^2);
  zy=(my-c)/sy;
  zx=(c-mx)/sx;

  dcdmx= (b^2 + a*b*radd^(-1/2)*(-1))/(b^2-1);
  dcdsx= (-2*a*b^2)/((b^2-1)^2*sx) +  (((b*(b^2+1))*(radd)^(1/2))/((b^2-1)^2*sx)  -  (sx*b*radd^(-1/2))/(b^2-1)*(log(b^2)+b^(2)-1) );
  dcdmy= (-1 + a*b*radd^(-1/2))/(b^2-1);
  dcdsy= (2*a*b)/((b^2-1)^2*sx) +  (((-b^2-1)*(radd)^(1/2))/((b^2-1)^2*sx)  +  (sy*b*radd^(-1/2))/(b^2-1)*(log(b^2)+1-b^(-2)) );

  d1=-sx^(-1)*dnorm(zx)+dcdmx*(sx^(-1)*dnorm(zx)-(sy^(-1)*dnorm(zy) ));
  d2=-zx*sx^(-1)*dnorm(zx)+dcdsx*(sx^(-1)*dnorm(zx)-(sy^(-1)*dnorm(zy) )) ;
  d3=sy^(-1)*dnorm(zy)+dcdmy*(sx^(-1)*dnorm(zx)-(sy^(-1)*dnorm(zy) )) ;
  d4=-zy*sy^(-1)*dnorm(zy)+dcdsy*(sx^(-1)*dnorm(zx)-(sy^(-1)*dnorm(zy) ));


  #varJ=t(c(d1, d2, d3, d4))%*%S%*%[d1 d2 d3 d4]';
  varJ=t(c(d1,d2,d3,d4))%*% S %*% t(t(c(d1,d2,d3,d4)))

#  Jstar=(Jhat+1)/2;
#  Jstar=log(Jstar/(1-Jstar));
#  VarJstar=(4/(Jhat^2 - 1)^2)*varJ

    Jstar=qnorm(Jhat);
    varJstar=(1/(dnorm(qnorm(Jhat))))^2*varJ
  
  ZJstar=Jstar/sqrt(varJstar)
  pvalJ=2*pnorm(-abs(ZJstar))

  CIstar=c(Jstar-Za*sqrt(varJstar), Jstar+Za*sqrt(varJstar));
  CI=c(0,0)
  CI[1]=exp(CIstar[1])/(1+exp(CIstar[1]))*2-1
  CI[2]=exp(CIstar[2])/(1+exp(CIstar[2]))*2-1
  CIJ=CI;
  Jhat
  CIwidth=max(CI)-min(CI)




  boots=0;c1hat_boots=0;c1hat_boots_or=0;  c1hat=cutoff;
  #===========CIs for the cutoff==================
  for (boots in 1:1000){
    #if (boots>=2){print((m1hat_boots*(b^2-1)-a+b*sqrt(a^2+(b^2-1)*s1hat_boots^2*log(b^2)))/(b^2-1))}
    #tryCatch({

      atx = sample(1:n1,n1,replace=T)
      aty = sample(1:n2,n2,replace=T)
      
      if (printProgress) {
        if (boots==100){print(boots);print("bootstap samples out of 1000...")}
        if (boots==200){print(boots);print("bootstap samples out of 1000...")}
        if (boots==300){print(boots);print("bootstap samples out of 1000...")}
        if (boots==400){print(boots);print("bootstap samples out of 1000...")}
        if (boots==500){print(boots);print("bootstap samples out of 1000...")}
        if (boots==600){print(boots);print("bootstap samples out of 1000...")}
        if (boots==700){print(boots);print("bootstap samples out of 1000...")}
        if (boots==800){print(boots);print("bootstap samples out of 1000...")}
        if (boots==900){print(boots);print("bootstap samples out of 1000...")}
        if (boots==1000){print(boots);print("Bootstrapping complete. You can     request all results through the summary() function.")}
      }

      xboots=xor[atx]
      yboots=yor[aty]

      ccb=boxcoxleo2(xboots,yboots)

      lam_boots=ccb$transformation.parameter
      transx1_boots=ccb$transx;
      transx2_boots=ccb$transy;


      m1hat_boots=mean(transx1_boots)
      m2hat_boots=mean(transx2_boots)
      s1hat_boots=std(transx1_boots)
      s2hat_boots=std(transx2_boots)

      b=s2hat_boots/s1hat_boots
      a=m2hat_boots-m1hat_boots
      #print(b)
      #print(a)
      c1hat_boots[boots]=(m1hat_boots*(b^2-1)-a+b*sqrt(a^2+(b^2-1)*s1hat_boots^2*log(b^2)))/(b^2-1)
      kr=(m1hat_boots*(b^2-1)-a+b*sqrt(a^2+(b^2-1)*s1hat_boots^2*log(b^2)))/(b^2-1)
      c1hat_boots_or[boots]  = ((kr*lam_boots+1))^(1/lam_boots)

    #}, error=function(e){})
  }


  varboots=var(na.omit(c1hat_boots_or))
  ccc=((c1hat*lam+1))^(1/lam)
  CI_BCAN=c(ccc-Za*sqrt(var(na.omit(c1hat_boots_or))),  ccc+Za*sqrt(var(na.omit(c1hat_boots_or))))
  CIwidth=max(CI_BCAN)-min(CI_BCAN)
  CIcutoff=CI_BCAN;
  cutoff=ccc;

  allc1=c1hat_boots_or

  #CI_BCPB=c(quantile(na.omit(c1hat_boots_or),0.025,"type"=2), quantile(na.omit(c1hat_boots_or),0.975,"type"=2))
  #CIwidth=max(CI_BCPB)-min(CI_BCPB)
  #CI_BCPB

  #legend("bottomright", legend=c("ROC estimate", "max of Youden index (J)", "rectangular 95% confidence region of the optimal operating point", "rectangular 95% confidence region of the optimal operating point"),
  #       col=c("red", "blue", "black", "green"), lty=c(1,1,1,1), cex=0.8)

  if (plots=="on"){
    legend("bottomright", legend=c(paste("ROC estimate with AUC = ", 
                                         round(auc,4), 
                                         ", CI: (", 
                                         round(CIauc[1],4), 
                                         ", ", 
                                         round(CIauc[2],4) ,")", sep = ""),
                                   paste("Maximized Youden index = ",round(Jhat,4), 
                                         ", CI: (", 
                                         round((CIJ[1]),4), 
                                         ", ", 
                                         round((CIJ[2]),4), ")", sep = ""),
                                   paste("Optimal pair of (FPR,TPR): (",
                                         round(1-Sp,4),
                                         ", ",
                                         round(Se,4), ")", sep = ""),
                                   paste("Area of the rect. conf. region =",round(arearect,4)),
                                   paste("Area of the egg conf. region =",round(areaegg,4)),
                                   paste("Youden based cutoff: ",round(cutoff,4),", CI: (", 
                                         round(CIcutoff[1],4),
                                         ", ",
                                         round(CIcutoff[2],4), ")", sep = ""),
                                   paste("Sp Marginal CI: (", 
                                         c(round(margcisp[1],4)), ", ", round(margcisp[2],4),
                                         ")", sep = ""),
                                   paste("Se Marginal CI: (", 
                                         c(round(margcise[1],4)), ", ", round(margcise[2],4), ")", 
                                         sep = "")),
           col=c("red", "blue", "red", "black", "green", "white", "white", "white"), 
           lty=c(1,1,1,1,1,1,NA, NA), 
           pch = c(NA, NA, 19, NA, NA, NA, NA, NA), cex=0.6)

  }



  #return(list(transx=((x^lam)-1)/lam, transy=((y^lam)-1)/lam , transformation.parameter=lam, AUC=auc, AUCCI=CIauc, J=Jhat, JCI=CIJ, Sens=Se, CImarginalSens=margcise, Spec=Sp, CImarginalSpec=margcisp, cutoff=cutoff, CIcutoff=CIcutoff, areaegg=areaegg, arearect=arearect, mxlam=mean(transx), sxlam=std(transx), mylam=mean(transy), sylam=std(transy)))

  res <- matrix(c(auc,CIauc[1],CIauc[2],
                  Jhat,CIJ[1],CIJ[2],
                  cutoff,CIcutoff[1],CIcutoff[2],
                  Sp,margcisp[1],margcisp[2],
                  Se,margcise[1],margcise[2]),ncol=3,byrow=TRUE)
  colnames(res) <- c("Estimate","Confidence","Interval")
  rownames(res) <- c("AUC","Jhat","cutoff","Sp","Se")
  res <- as.table(res)
  res

  return(list(transx=((x^lam)-1)/lam, transy=((y^lam)-1)/lam , transformation.parameter=lam, AUC=auc, AUCCI=CIauc, pvalueAUC=pvalauc, J=Jhat, JCI=CIJ, pvalueJ=pvalJ, Sens=Se, CImarginalSens=margcise, Spec=Sp, CImarginalSpec=margcisp, cutoff=cutoff, CIcutoff=CIcutoff,  areaegg=areaegg, arearect=arearect, mxlam=mean(transx), sxlam=std(transx), mylam=mean(transy), sylam=std(transy), results=res , rocfun=rocfun))
  }
}
  

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
rocboxcoxCI<-function(marker, D, givenSP, givenSE, alpha, plots){

  if (plots!="on"){plots="off"}
  if (length(marker) != length(D)) {
    stop("ERROR: The length of the 'marker' and 'D' inputs must be equal.")
  } else if (min(D) != 0 | max(D) != 1) {
    stop("ERROR: Controls must be assigned a value of 0; cases must be assigned a value of 1. Both controls and cases should be included in the dataset.")
  } else if (sum(is.na(marker)) > 0 | sum(is.na(D)) > 0) {
    stop("ERROR: Please remove all missing data before running this function.")
  } else if (alpha <= 0 | alpha >= 1) {
    stop("ERROR: The level of significance, alpha, should be set between 0 and 1. A common choice is 0.05.")
  } else {
  
#graphics.off()
  if (plots!="on"){plots="off"}

  xor=marker[D==0]
  yor=marker[D==1]

  x = xor
  y = yor

likbox<-function(x,y,h){
  n=length(x);
  m=length(y);
  out=c();
  for (i in 1:length(h)){
    #print(i)
    if (h[i]==0){
      xh=log(x);
      yh=log(y);
    } else {
      xh=((x^h[i])-1)/h[i];
      yh=((y^h[i])-1)/h[i];
    }


    out[i]<-c( -n/2*log(sum((xh-sum(xh)/n)^2)/n)  -m/2*log(sum((yh-sum(yh)/m)^2)/m) +(h[i]-1)*(sum(log(x))+sum(log(y))))
  }
  return(out)
}





boxcoxleo<-function(x,y){
  init=1;
  logL<-function(h){
    -likbox(x,y,h)
  }
  #lam=fminsearch(logL,init)
  #lam=c(lam$optbase$xopt)
  lam<- optim(1,logL,gr=NULL,method="BFGS", control=list(maxit=10000))
  lam=c(lam$par)
  transx=((x^lam)-1)/lam
  transy=((y^lam)-1)/lam

  #test1=print(shapiro.test(transx))
  #test2=print(shapiro.test(transy))


  return(list(transformation.parameter=lam,transx=((x^lam)-1)/lam, transy=((y^lam)-1)/lam ))
}


boxcoxleo2<-function(x,y){
  init=1;
  logL<-function(h){
    -likbox(x,y,h)
  }
  lam<- optim(1,logL,gr=NULL,method="BFGS", control=list(maxit=10000))
  lam=c(lam$par)
  #lam=fminsearch(logL,init)
  #lam=c(lam$optbase$xopt)
  transx=((x^lam)-1)/lam
  transy=((y^lam)-1)/lam



  return(list(transformation.parameter=lam,transx=((x^lam)-1)/lam, transy=((y^lam)-1)/lam ))
}


cc=boxcoxleo2(x,y)

transx=cc$transx
transy=cc$transy

#x11()
#qqnorm(x)
#title(main="                                         for X")
#x11()
#qqnorm(y)
#title(main="                                         for Y")




roc<-function(t,x,y){
  1-pnorm(qnorm(1-t,mean=mean(transx),sd=std(transx)),mean=mean(transy),sd=std(transy))
}

rocuseless<-function(t){
  1-pnorm(qnorm(1-t,mean=1,sd=1),mean=1,sd=1)
}

# x11()
#plot.new()
txt <- paste("Box-Cox Based ROC, alpha =", round(alpha,3) ) 
plot(linspace(0,1,1000),roc(linspace(0,1,1000),transx,transy),main=" ",xlab="FPR",ylab="TPR",type="l",col="red")
lines(linspace(0,1,10),linspace(0,1,10),type="l", lty=2)
title(main=txt)
m1hat=mean(transx)
m2hat=mean(transy)
s1hat=std(transx)
s2hat=std(transy)
n1=length(x)
n2=length(y)



I=zeros(5,5);
sh=sqrt(1/(length(transx))*sum((transx-mean(transx))^2))
sd=sqrt(1/(length(transy))*sum((transy-mean(transy))^2))

mh=m1hat;
md=m2hat;
n=n1;
m=n2;
xlam=transx;
ylam=transy;
lam=cc$transformation.parameter
I[1,1]=n/sh^2;
I[2,2]=-(n/sh^2-3/sh^4*sum((xlam-mh)^2));
I[3,3]=m/sd^2;
I[4,4]=-(m/sd^2-3/sd^4*sum((ylam-md)^2));

kk=  sum(((mh - (x^lam - 1)/lam)*((2*(x^lam - 1))/lam^3 - (2*x^lam*log(x))/lam^2 + (x^lam*log(x)^2)/lam))/sh^2) + - sum(((y^lam - 1)/lam^2 - (y^lam*log(y))/lam)^2/sd^2) + - sum(((x^lam - 1)/lam^2 - (x^lam*log(x))/lam)^2/sh^2) + + sum(((md - (y^lam - 1)/lam)*((2*(y^lam - 1))/lam^3 - (2*y^lam*log(y))/lam^2 + (y^lam*log(y)^2)/lam))/sd^2);
I[5,5]=-kk ;

I[1,5]=sum(((2*(x^lam - 1))/lam^2 - (2*x^lam*log(x))/lam)/(2*sh^2));
I[2,5]=-sum((2*((x^lam - 1)/lam^2 - (x^lam*log(x))/lam)*(mh - (x^lam - 1)/lam))/sh^3);
I[3,5]=-sum(-((2*(y^lam - 1))/lam^2 - (2*y^lam*log(y))/lam)/(2*sd^2));


I[4,5]=-sum((2*((y^lam - 1)/lam^2 - (y^lam*log(y))/lam)*(md - (y^lam - 1)/lam))/sd^3);



I[5,1]=I[1,5]
I[5,2]=I[2,5]
I[5,3]=I[3,5]
I[5,4]=I[4,5]

S=inv(I)

S=S[1:4,1:4]



#=============================================
#GRID LINES:

kk=1;CIROC=c(1,1);CIROCvec1=0;CIROCvec2=0;Sehat=0;
tt=c(linspace(0,1,2000))
for (t in tt){
  norminvROC= ((m2hat-m1hat)/s2hat+s1hat/s2hat*qnorm(t))
  #====================ROC1==============================
  dm1 =-1/s2hat;
  ds1 =-(2^(1/2)*erfcinv(2*t))/s2hat;
  dm2 =1/s2hat;
  ds2 =(m1hat - m2hat)/s2hat^2 + (2^(1/2)*s1hat*erfcinv(2*t))/s2hat^2

  vROC=t(c(dm1,ds1,dm2,ds2))%*% S %*% t(t(c(dm1,ds1,dm2,ds2)))

  CIinvROC=c(norminvROC-qnorm(1-alpha/2)*sqrt(vROC), norminvROC+qnorm(1-alpha/2)*sqrt(vROC))
  CIROC=c( pnorm(CIinvROC[1]), pnorm(CIinvROC[2]))
  CIROCvec1[kk]=CIROC[1]
  CIROCvec2[kk]=CIROC[2]

  Sehat[kk]=pnorm(norminvROC)

  #points(t,Sehat[length(Sehat)], col = "red")
  #lines(c(t,t),c(CIROC[1],CIROC[1]),col="black",lty=2)
  #lines(c(t,t),c(CIROC[2],CIROC[2]),col="black")

  kk=kk+1
}

if (plots=="on"){
lines(c(tt),c(CIROCvec1),col="black",lty=2)
lines(c(tt),c(CIROCvec2),col="black",lty=2)


#====END GRID LINES=========================
#=============================================

# x11()
# plot.new()
txt <- paste("Box-Cox Based ROC, alpha =", round(alpha,3) ) 
plot(linspace(0,1,1000),roc(linspace(0,1,1000),transx,transy),main=" ",xlab="FPR",ylab="TPR",type="l",col="red")
lines(linspace(0,1,10),linspace(0,1,10),type="l", lty=2)
title(main=txt)
}
kk=1;CIROC=c(1,1);CIROCvec1=0;CIROCvec2=0;Sehat=0;
tt=c(givenSP);tt=1-givenSP;
#tt=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
for (t in tt){
norminvROC= ((m2hat-m1hat)/s2hat+s1hat/s2hat*qnorm(t))
#====================ROC1==============================
dm1 =-1/s2hat;
ds1 =-(2^(1/2)*erfcinv(2*t))/s2hat;
dm2 =1/s2hat;
ds2 =(m1hat - m2hat)/s2hat^2 + (2^(1/2)*s1hat*erfcinv(2*t))/s2hat^2

vROC=t(c(dm1,ds1,dm2,ds2))%*% S %*% t(t(c(dm1,ds1,dm2,ds2)))

CIinvROC=c(norminvROC-qnorm(1-alpha/2)*sqrt(vROC), norminvROC+qnorm(1-alpha/2)*sqrt(vROC))
CIROC=c( pnorm(CIinvROC[1]), pnorm(CIinvROC[2]))
CIROCvec1[kk]=CIROC[1]
CIROCvec2[kk]=CIROC[2]

Sehat[kk]=pnorm(norminvROC)
if (plots=="on"){
points(t,Sehat[length(Sehat)], col = "green")
lines(c(t,t),c(CIROC[1],CIROC[2]),col="green")
}
kk=kk+1
}


Spvalues=1-tt;FPRvalues=tt;
SEandCIs = matrix( c(Spvalues,FPRvalues, Sehat, CIROCvec1,CIROCvec2), nrow=length(CIROCvec1), ncol=5)
colnames(SEandCIs)  <- c("Given Sp","Given FPR", "Sehat","LL of 95%CI","UL of 95%CI")
SEandCIs

CIlowSe=CIROCvec1
CIuppSe=CIROCvec2
CIse=t(rbind(CIlowSe,CIuppSe))







#==================ROCinv1 (CIs for Sp)===================

mm1=m1hat;
mm2=m2hat;
ss1=s1hat;
ss2=s2hat;
m1hat=mm2;
m2hat=mm1;
s1hat=ss2;
s2hat=ss1;
kk=1;CIROC=c(1,1);CIROCvec1=0;CIROCvec2=0;Sphat=0;

mh=m1hat;
md=m2hat;
n=n1;
m=n2;
xlam=transx;
ylam=transy;
lam=cc$transformation.parameter
I[1,1]=n/sh^2;
I[2,2]=-(n/sh^2-3/sh^4*sum((xlam-mh)^2));
I[3,3]=m/sd^2;
I[4,4]=-(m/sd^2-3/sd^4*sum((ylam-md)^2));

kk=  sum(((mh - (x^lam - 1)/lam)*((2*(x^lam - 1))/lam^3 - (2*x^lam*log(x))/lam^2 + (x^lam*log(x)^2)/lam))/sh^2) + - sum(((y^lam - 1)/lam^2 - (y^lam*log(y))/lam)^2/sd^2) + - sum(((x^lam - 1)/lam^2 - (x^lam*log(x))/lam)^2/sh^2) + + sum(((md - (y^lam - 1)/lam)*((2*(y^lam - 1))/lam^3 - (2*y^lam*log(y))/lam^2 + (y^lam*log(y)^2)/lam))/sd^2);
I[5,5]=-kk ;

I[1,5]=sum(((2*(x^lam - 1))/lam^2 - (2*x^lam*log(x))/lam)/(2*sh^2));
I[2,5]=-sum((2*((x^lam - 1)/lam^2 - (x^lam*log(x))/lam)*(mh - (x^lam - 1)/lam))/sh^3);
I[3,5]=-sum(-((2*(y^lam - 1))/lam^2 - (2*y^lam*log(y))/lam)/(2*sd^2));


I[4,5]=-sum((2*((y^lam - 1)/lam^2 - (y^lam*log(y))/lam)*(md - (y^lam - 1)/lam))/sd^3);



I[5,1]=I[1,5]
I[5,2]=I[2,5]
I[5,3]=I[3,5]
I[5,4]=I[4,5]

S=inv(I)

S=S[1:4,1:4]


pp=c(givenSE);
#pp=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
kk=1;CIROC=c(1,1);CIROCvec1=0;CIROCvec2=0;Sphat=0;

for (t in pp){
  norminvROC= ((m2hat-m1hat)/s2hat+s1hat/s2hat*qnorm(t))
  #====================ROC1==============================
  dm1 =-1/s2hat;
  ds1 =-(2^(1/2)*erfcinv(2*t))/s2hat;
  dm2 =1/s2hat;
  ds2 =(m1hat - m2hat)/s2hat^2 + (2^(1/2)*s1hat*erfcinv(2*t))/s2hat^2

  vROC=t(c(dm1,ds1,dm2,ds2))%*% S %*% t(t(c(dm1,ds1,dm2,ds2)))

  CIinvROC=c(norminvROC-qnorm(1-alpha/2)*sqrt(vROC), norminvROC+qnorm(1-alpha/2)*sqrt(vROC))
  CIROC=c( pnorm(CIinvROC[1]), pnorm(CIinvROC[2]))
  CIROCvec1[kk]=CIROC[1]
  CIROCvec2[kk]=CIROC[2]

  Sphat[kk]=pnorm(norminvROC)
  if (plots=="on"){
  points(Sphat[length(Sphat)],t, col = "purple")
  lines(c(CIROC[1],CIROC[2]),c(t,t),col="purple")
  }
  kk=kk+1
}




#SPandCIs = matrix( c(Sphat, CIROCvec1,CIROCvec2), nrow=length(CIROCvec1), ncol=3)
#colnames(SPandCIs)  <- c("Sphat","LL of 95%CI","UL of 95%CI")
#SPandCIs


Sevalues=givenSE;
SPandCIs = matrix( c(Sevalues, Sphat, CIROCvec1,CIROCvec2), nrow=length(CIROCvec1), ncol=4)
colnames(SPandCIs)  <- c("Given Se", "Sphat","LL of 95%CI","UL of 95%CI")
SPandCIs

CIlowSp=CIROCvec1
CIuppSp=CIROCvec2
CIsp=t(rbind(CIlowSp,CIuppSp))
return(list(SPandCIs=SPandCIs, SEandCIs=SEandCIs, Sevalues=Sevalues, Sphat=Sphat, CIsp=CIsp, Spvalues=Spvalues, Sphat=Sphat, CIse=CIse ))


}

}


## ---- echo=FALSE, eval=TRUE---------------------------------------------------
#====THIS IS THE MAIN FUNCTION CHECKBOXCOX ====================================================
#INPUT ARGUMENTS: TWO GROUPS (x,y) and plots="on" or "off"====================================


checkboxcox<-function(marker, D, plots, printShapiro = FALSE){
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if (plots!="on"){plots="off"}
  if (length(marker) != length(D)) {
    stop("ERROR: The length of the 'marker' and 'D' inputs must be equal.")
  } else if (min(D) != 0 | max(D) != 1) {
    stop("ERROR: Controls must be assigned a value of 0; cases must be assigned a value of 1. Both controls and cases should be included in the dataset.")
  } else if (sum(is.na(marker)) > 0 | sum(is.na(D)) > 0) {
    stop("ERROR: Please remove all missing data before running this function.")
  } else {
    
    statusD=D;
    xor=marker[D==0]
    yor=marker[D==1]
    
    x = xor
    y = yor
    scores = c(x, y)
  
  #====THIS IS THE LIKELIHOOD OF THE BOXCOX TRANSFORMATION=======================================
  #INPUT ARGUMENTS: TWO GROUPS (x,y) and plots="on" or "off"====================================
  
  likbox<-function(x,y,h){
    n=length(x);
    m=length(y);
    out=c();
    for (i in 1:length(h)){
      #print(i)
      if (h[i]==0){
        xh=log(x);
        yh=log(y);
      } else {
        xh=((x^h[i])-1)/h[i];
        yh=((y^h[i])-1)/h[i];
      }
  
  
      out[i]<-c( -n/2*log(sum((xh-sum(xh)/n)^2)/n)  -m/2*log(sum((yh-sum(yh)/m)^2)/m) +(h[i]-1)*(sum(log(x))+sum(log(y))))
    }
    return(out)
  }
  
  
  
  
  #====THIS THE FUNCTION NEEDED TO APPLY THE BOXCOX TRANSFORMATION INPUT ARGUMENTS:
  #====TWO GROUPS (x,y) and plots="on" or "off"====================================
  boxcoxleo<-function(x,y,plots){
    if (plots!="on"){plots="off"}
    init=1;
    logL<-function(h){
      -likbox(x,y,h)
    }
    #lam=fminsearch(logL,init)
    lam <- optimize(logL, c(-100, 100), tol = 0.0001)
    lam=c(lam$minimum) #c(lam$optbase$xopt)
    transx=((x^lam)-1)/lam
    transy=((y^lam)-1)/lam
  
  
  
    if (!printShapiro) {
      test1=shapiro.test(x)
      test2=shapiro.test(y)
      test3=shapiro.test(transx)
      test4=shapiro.test(transy)
    }
  
  
    if (printShapiro) {
      test1=print(shapiro.test(x))
      test2=print(shapiro.test(y))
  
      test3=print(shapiro.test(transx))
      test4=print(shapiro.test(transy))
    }
  
    pval1=test1$p.value
    pval2=test2$p.value
    pval3=test3$p.value
    pval4=test4$p.value
    pvalues=c(pval1,pval2,pval3,pval4)
  
    return(list(transformation.parameter=lam,transx=((x^lam)-1)/lam, transy=((y^lam)-1)/lam,pvalues=pvalues))
  }
  
  
  #====APPLY THE BOXCOX TRANSFORMATION WITH OR WITHOUT THE PLOTS================
  cc=boxcoxleo(x,y,plots)
  
  lam=cc$transformation.parameter
  pvalues=cc$pvalues
  
  transx=cc$transx #transformed scores for healthy
  transy=cc$transy #transformed scores for diseased
  
  if (plots!="on"){plots="off"}
  
  #====IF PLOTS ARE REQUESTED PROVIDE THE HISTOGRAMS AND QQ-PLOTS FOR EACH GROUP
  #====BEFORE AND AFTER THE TRANSFORMATION=======================================
  if (plots=="on") {
    # x11()
  par(mfrow=c(2,2))
    
  par(mar = c(4.1, 3.1, 3.1, 1.1))
    
  hist(x, main = "Histogram of X", xlab = "X")
  qqnorm(x, main = "Normal Q-Q Plot of X")
  qqline(x, col = 2,lwd=2,lty=2)
  
  hist(transx, main = "Histogram of Transformed X", xlab = "Transformed X")
  qqnorm(transx, main = "Normal Q-Q Plot of Transformed X")
  qqline(transx, col = 2,lwd=2,lty=2)
  
  par(mfrow=c(2,2))
  
  hist(y, main = "Histogram of Y", xlab = "Y")
  qqnorm(y, main = "Normal Q-Q Plot of Y")
  qqline(y, col = 2,lwd=2,lty=2)
  
  hist(transy, main = "Histogram of Transformed Y", xlab = "Transformed Y")
  qqnorm(transy, main = "Normal Q-Q Plot of Transformed Y")
  qqline(transy, col = 2,lwd=2,lty=2)

  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  }
  
  
  #====TWO ROC FUNCTIONS: ONE IS THE BOXCOX AND THE OTHER THE REFERENCE LINE================
  roc<-function(t){
    1-pnorm(qnorm(1-t,mean=mean(transx),sd=std(transx)),mean=mean(transy),sd=std(transy))
  }
  
  rocuseless<-function(t){
    1-pnorm(qnorm(1-t,mean=1,sd=1),mean=1,sd=1)
  }
  
  
  
  #================IF PLOTS ARE REQUESTED PLOT THE ROCS==
  if (plots=="on") {
  # x11()
    par(mfrow = c(1, 1))
  plot(linspace(0,1,1000),roc(linspace(0,1,1000)),main="Empirical and Box-Cox Based ROC",xlab="FPR = 1 - Specificity",ylab="TPR = Sensitivity",type="l",col="red")
  lines(linspace(0,1,1000),linspace(0,1,1000),type="l", lty=2)
  lines(roc.curve(scores, statusD), col=1)
  
  #================LEGEND OF THE ROC PLOT================
  legend("bottomright", legend=c(paste("Box-Cox ROC estimate "),
                                 paste("Empirical ROC estimate ")),
                                 col=c("red", "black"), lty=c(1, 1), pch = c(NA, NA), cex=0.8)
  
  }
  
  #================OUTPUT ARGUMENTS======================
  return(list(transformation.parameter=lam,transx=((x^lam)-1)/lam, transy=((y^lam)-1)/lam, pval_x1=pvalues[1], pval_x2=pvalues[2],pval_x1trans=pvalues[3], pval_x2trans=pvalues[4], rocbc=roc ))
  
  }
}


## ---- echo=FALSE, eval=TRUE---------------------------------------------------

comparebcAUC <-function (marker1, marker2, D, alpha, plots){
  
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
  
  init=c(0.8,0.8)
  lam=fminsearch(logL,init)
  lam=c(lam$optbase$xopt)
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
  h11_11=sum(((2*((W1a^lam1 - 1)/lam1^2 - (W1a^lam1*log(W1a))/lam1)^2)/s1ahat^2 - ((2*m1ahat - (2*(W1a^lam1 - 1))/lam1)*((2*W1a^lam1 - 2)/lam1^3 - (2*W1a^lam1*log(W1a))/lam1^2 + (W1a^lam1*log(W1a)^2)/lam1))/s1ahat^2 + (2*cova*(m2ahat - (W2a^lam2 - 1)/lam2)*((2*W1a^lam1 - 2)/lam1^3 - (2*W1a^lam1*log(W1a))/lam1^2 + (W1a^lam1*log(W1a)^2)/lam1))/(s1ahat^2*s2ahat^2))/((2*cova^2)/(s1ahat^2*s2ahat^2) - 2) + ((2*((W1b^lam1 - 1)/lam1^2 - (W1b^lam1*log(W1b))/lam1)^2)/s1bhat^2 - ((2*m1bhat - (2*(W1b^lam1 - 1))/lam1)*((2*W1b^lam1 - 2)/lam1^3 - (2*W1b^lam1*log(W1b))/lam1^2 + (W1b^lam1*log(W1b)^2)/lam1))/s1bhat^2 + (2*covb*(m2bhat - (W2b^lam2 - 1)/lam2)*((2*W1b^lam1 - 2)/lam1^3 - (2*W1b^lam1*log(W1b))/lam1^2 + (W1b^lam1*log(W1b)^2)/lam1))/(s1bhat^2*s2bhat^2))/((2*covb^2)/(s1bhat^2*s2bhat^2) - 2))
  h11_12=sum((2*cova*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1))/(lam1^2*lam2^2*(- 2*cova^2 + 2*s1ahat^2*s2ahat^2)) + (2*covb*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1))/(lam1^2*lam2^2*(- 2*covb^2 + 2*s1bhat^2*s2bhat^2)))
  
  
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
  h12_11=sum((2*cova*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1))/(lam1^2*lam2^2*(- 2*cova^2 + 2*s1ahat^2*s2ahat^2)) + (2*covb*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1))/(lam1^2*lam2^2*(- 2*covb^2 + 2*s1bhat^2*s2bhat^2)))
  h12_12=sum(((2*((W2a^lam2 - 1)/lam2^2 - (W2a^lam2*log(W2a))/lam2)^2)/s2ahat^2 - ((2*m2ahat - (2*(W2a^lam2 - 1))/lam2)*((2*W2a^lam2 - 2)/lam2^3 - (2*W2a^lam2*log(W2a))/lam2^2 + (W2a^lam2*log(W2a)^2)/lam2))/s2ahat^2 + (2*cova*(m1ahat - (W1a^lam1 - 1)/lam1)*((2*W2a^lam2 - 2)/lam2^3 - (2*W2a^lam2*log(W2a))/lam2^2 + (W2a^lam2*log(W2a)^2)/lam2))/(s1ahat^2*s2ahat^2))/((2*cova^2)/(s1ahat^2*s2ahat^2) - 2) + ((2*((W2b^lam2 - 1)/lam2^2 - (W2b^lam2*log(W2b))/lam2)^2)/s2bhat^2 - ((2*m2bhat - (2*(W2b^lam2 - 1))/lam2)*((2*W2b^lam2 - 2)/lam2^3 - (2*W2b^lam2*log(W2b))/lam2^2 + (W2b^lam2*log(W2b)^2)/lam2))/s2bhat^2 + (2*covb*(m1bhat - (W1b^lam1 - 1)/lam1)*((2*W2b^lam2 - 2)/lam2^3 - (2*W2b^lam2*log(W2b))/lam2^2 + (W2b^lam2*log(W2b)^2)/lam2))/(s1bhat^2*s2bhat^2))/((2*covb^2)/(s1bhat^2*s2bhat^2) - 2))
  
  
  
  
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
  #DEALING WITH ORIGINAL AUCS
  #############################################################################
  AUC1 =  erfc((2^(1/2)*(m1ahat - m1bhat))/(2*(s1ahat^2 + s1bhat^2)^(1/2)))/2
  dm1a =  -(2^(1/2)*exp(-(m1ahat - m1bhat)^2/(2*(s1ahat^2 + s1bhat^2))))/(2*pi^(1/2)*(s1ahat^2 + s1bhat^2)^(1/2))
  dm1b =  (2^(1/2)*exp(-(m1ahat - m1bhat)^2/(2*(s1ahat^2 + s1bhat^2))))/(2*pi^(1/2)*(s1ahat^2 + s1bhat^2)^(1/2))
  ds1a =  (2^(1/2)*s1ahat*exp(-(m1ahat - m1bhat)^2/(2*(s1ahat^2 + s1bhat^2)))*(m1ahat - m1bhat))/(2*pi^(1/2)*(s1ahat^2 + s1bhat^2)^(3/2))
  ds1b =  (2^(1/2)*s1bhat*exp(-(m1ahat - m1bhat)^2/(2*(s1ahat^2 + s1bhat^2)))*(m1ahat - m1bhat))/(2*pi^(1/2)*(s1ahat^2 + s1bhat^2)^(3/2))
  
  AUC2 =  erfc((2^(1/2)*(m2ahat - m2bhat))/(2*(s2ahat^2 + s2bhat^2)^(1/2)))/2
  dm2a =  -(2^(1/2)*exp(-(m2ahat - m2bhat)^2/(2*(s2ahat^2 + s2bhat^2))))/(2*pi^(1/2)*(s2ahat^2 + s2bhat^2)^(1/2))
  dm2b =  (2^(1/2)*exp(-(m2ahat - m2bhat)^2/(2*(s2ahat^2 + s2bhat^2))))/(2*pi^(1/2)*(s2ahat^2 + s2bhat^2)^(1/2))
  ds2a =  (2^(1/2)*s2ahat*exp(-(m2ahat - m2bhat)^2/(2*(s2ahat^2 + s2bhat^2)))*(m2ahat - m2bhat))/(2*pi^(1/2)*(s2ahat^2 + s2bhat^2)^(3/2))
  ds2b =  (2^(1/2)*s2bhat*exp(-(m2ahat - m2bhat)^2/(2*(s2ahat^2 + s2bhat^2)))*(m2ahat - m2bhat))/(2*pi^(1/2)*(s2ahat^2 + s2bhat^2)^(3/2))
  
  varAUC1=c(dm1a, ds1a, dm1b, ds1b)%*%VV[1:4,1:4]%*%t(t(c(dm1a, ds1a, dm1b, ds1b)));
  varAUC2=c(dm2a, ds2a, dm2b, ds2b)%*%VV[5:8,5:8]%*%t(t(c(dm2a, ds2a, dm2b, ds2b)));
  
  varAUC1=c(varAUC1)
  varAUC2=c(varAUC2)
  
  V=zeros(8,8);
  covAUC=c(dm1a, ds1a, dm1b, ds1b, 0, 0, 0, 0)%*%VV[1:8,1:8]%*%t(t(c(0, 0, 0, 0, dm2a, ds2a, dm2b, ds2b)))
  covAUC=c(covAUC)
  
  
  ZdAUC=(AUC2-AUC1)/(sqrt(varAUC1+varAUC2-2*covAUC))
  SEdAUC=sqrt(varAUC1+varAUC2-2*covAUC)
  CIdiffdAUC=c(AUC2-AUC1-pmalpha*SEdAUC, AUC2-AUC1+pmalpha*SEdAUC)
  pval2tdAUC=2*pnorm(-abs(ZdAUC))
  
  
  #############################################################################
  #PROBIT VERSIONS
  #############################################################################
  
  
  #norminv versions
  AUC1= ((m1bhat-m1ahat)/s1bhat)/(sqrt(1+(s1ahat/s1bhat)^2));
  AUC2= ((m2bhat-m2ahat)/s2bhat)/(sqrt(1+(s2ahat/s2bhat)^2));
  AUC1MC=AUC1
  AUC2MC=AUC2
  AUC1original=pnorm(AUC1)
  AUC2original=pnorm(AUC2)
  
  
  #====================norminv AUC1==============================
  dm1a =-1/(s1bhat*sqrt(1+(s1ahat/s1bhat)^2));
  ds1a = (s1ahat*(m1ahat-m1bhat))/(s1bhat^3*(1+(s1ahat/s1bhat)^2)^(3/2))
  dm1b =1/(s1bhat*sqrt(1+(s1ahat/s1bhat)^2));
  ds1b =((m1ahat-m1bhat))/(s1bhat^2*(1+(s1ahat/s1bhat)^2)^(3/2))
  
  varAUC1=c(dm1a, ds1a, dm1b, ds1b)%*%VV[1:4,1:4]%*%t(t(c(dm1a, ds1a, dm1b, ds1b)))
  
  
  #====================norminv AUC2==============================
  dm2a =-1/(s2bhat*sqrt(1+(s2ahat/s2bhat)^2));
  ds2a = (s2ahat*(m2ahat-m2bhat))/(s2bhat^3*(1+(s2ahat/s2bhat)^2)^(3/2))
  dm2b =1/(s2bhat*sqrt(1+(s2ahat/s2bhat)^2));
  ds2b =((m2ahat-m2bhat))/(s2bhat^2*(1+(s2ahat/s2bhat)^2)^(3/2))
  
  varAUC2=c(dm2a, ds2a, dm2b, ds2b)%*%VV[5:8,5:8]%*%t(t(c(dm2a, ds2a, dm2b, ds2b)));
  
  
  #====================norminv cov(AUC1,AUC2)==============================
  
  V=zeros(8,8);
  covAUC=c(dm1a, ds1a, dm1b, ds1b, 0, 0, 0, 0)%*%VV[1:8,1:8]%*%t(t(c(0, 0, 0, 0, dm2a, ds2a, dm2b, ds2b)))
  covAUC=c(covAUC)
  
  
  Z=(AUC2MC-AUC1MC)/(sqrt(varAUC1+varAUC2-2*covAUC))
  SE=sqrt(varAUC1+varAUC2-2*covAUC)
  CIdiff=c(AUC2MC-AUC1MC-pmalpha*SE, AUC2MC-AUC1MC+pmalpha*SE)
  
  pval2t=2*pnorm(-abs(Z))
  
  
  #################### IF PLOTS ARE REQUESTED #################################
  
  #====TWO ROC FUNCTIONS: ONE IS THE BOXCOX AND THE OTHER THE REFERENCE LINE================
  roc1<-function(t){
    1-pnorm(qnorm(1-t,mean=mean(W1alam),sd=std(W1alam)),mean=mean(W1blam),sd=std(W1blam))
  }
  
  roc2<-function(t){
    1-pnorm(qnorm(1-t,mean=mean(W2alam),sd=std(W2alam)),mean=mean(W2blam),sd=std(W2blam))
  }
  
  rocuseless<-function(t){
    1-pnorm(qnorm(1-t,mean=1,sd=1),mean=1,sd=1)
  }
  
  
  
  #================IF PLOTS ARE REQUESTED PLOT THE ROCS==
  if (plots=="on") {
    # x11()
    plot(linspace(0,1,1000),roc1(linspace(0,1,1000)),main="Box-Cox Based ROCs",xlab="FPR = 1 - Specificity",ylab="TPR = Sensitivity",type="l",col="red")
    lines(linspace(0,1,1000),roc2(linspace(0,1,1000)),main="Box-Cox Based ROCs",xlab="FPR = 1 - Specificity",ylab="TPR = Sensitivity",type="l",col="black")
    lines(linspace(0,1,1000),linspace(0,1,1000),type="l", lty=2)
  
  
  
    legend("bottomright", legend=c(paste("ROC for Marker 1 with AUC =", round(AUC1original,4)),
                                   paste("ROC for Marker 2 with AUC =", round(AUC2original,4)),
                                   paste("P-value for the probit AUC difference:",round(pval2t,4), " ")),
  
           col=c("red", "black", "white"), lty=c(1, 1, NA), pch = c(NA, NA, NA), cex=0.8)
  
  }
  
  
  ###############################################################################
  ###############################################################################
  #================OUTPUT ARGUMENTS======================
  
  res <- matrix(c(round(AUC1original,4),round(AUC2original,7), round(pval2t,7) ,CIdiff[1],CIdiff[2]),ncol=5,byrow=TRUE)
  colnames(res) <- c("AUC 1:","    AUC 2:","    p-value (probit):","   CI probit (LL):","    CI probit (UL):")
  rownames(res) <- c("Estimates:")
  res <- as.table(res)
  res
  
  
  #return(list(AUCmarker1=AUC1original,AUCmarker2=AUC2original, pvalue_difference= pval2t, CI_difference= pnorm(CIdiff), rocbc1=roc1, rocbc2=roc2))
  return(list(resultstable=res,AUCmarker1=AUC1original,AUCmarker2=AUC2original, pvalue_probit_difference= pval2t, CI_probit_difference= CIdiff, pvalue_difference= pval2tdAUC, CI_difference= CIdiffdAUC, roc1=roc1, roc2=roc2, transx1=W1alam, transy1=W1blam, transx2=W2alam, transy2=W2blam))
  
  }
}





## ---- echo=FALSE, eval=TRUE---------------------------------------------------

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
      
        init=c(0.8,0.8)
        lam=fminsearch(logL,init)
        lam=c(lam$optbase$xopt)
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
        h11_11=sum(((2*((W1a^lam1 - 1)/lam1^2 - (W1a^lam1*log(W1a))/lam1)^2)/s1ahat^2 - ((2*m1ahat - (2*(W1a^lam1 - 1))/lam1)*((2*W1a^lam1 - 2)/lam1^3 - (2*W1a^lam1*log(W1a))/lam1^2 + (W1a^lam1*log(W1a)^2)/lam1))/s1ahat^2 + (2*cova*(m2ahat - (W2a^lam2 - 1)/lam2)*((2*W1a^lam1 - 2)/lam1^3 - (2*W1a^lam1*log(W1a))/lam1^2 + (W1a^lam1*log(W1a)^2)/lam1))/(s1ahat^2*s2ahat^2))/((2*cova^2)/(s1ahat^2*s2ahat^2) - 2) + ((2*((W1b^lam1 - 1)/lam1^2 - (W1b^lam1*log(W1b))/lam1)^2)/s1bhat^2 - ((2*m1bhat - (2*(W1b^lam1 - 1))/lam1)*((2*W1b^lam1 - 2)/lam1^3 - (2*W1b^lam1*log(W1b))/lam1^2 + (W1b^lam1*log(W1b)^2)/lam1))/s1bhat^2 + (2*covb*(m2bhat - (W2b^lam2 - 1)/lam2)*((2*W1b^lam1 - 2)/lam1^3 - (2*W1b^lam1*log(W1b))/lam1^2 + (W1b^lam1*log(W1b)^2)/lam1))/(s1bhat^2*s2bhat^2))/((2*covb^2)/(s1bhat^2*s2bhat^2) - 2))
        h11_12=sum((2*cova*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1))/(lam1^2*lam2^2*(- 2*cova^2 + 2*s1ahat^2*s2ahat^2)) + (2*covb*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1))/(lam1^2*lam2^2*(- 2*covb^2 + 2*s1bhat^2*s2bhat^2)))
      
      
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
        h12_11=sum((2*cova*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1))/(lam1^2*lam2^2*(- 2*cova^2 + 2*s1ahat^2*s2ahat^2)) + (2*covb*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1))/(lam1^2*lam2^2*(- 2*covb^2 + 2*s1bhat^2*s2bhat^2)))
        h12_12=sum(((2*((W2a^lam2 - 1)/lam2^2 - (W2a^lam2*log(W2a))/lam2)^2)/s2ahat^2 - ((2*m2ahat - (2*(W2a^lam2 - 1))/lam2)*((2*W2a^lam2 - 2)/lam2^3 - (2*W2a^lam2*log(W2a))/lam2^2 + (W2a^lam2*log(W2a)^2)/lam2))/s2ahat^2 + (2*cova*(m1ahat - (W1a^lam1 - 1)/lam1)*((2*W2a^lam2 - 2)/lam2^3 - (2*W2a^lam2*log(W2a))/lam2^2 + (W2a^lam2*log(W2a)^2)/lam2))/(s1ahat^2*s2ahat^2))/((2*cova^2)/(s1ahat^2*s2ahat^2) - 2) + ((2*((W2b^lam2 - 1)/lam2^2 - (W2b^lam2*log(W2b))/lam2)^2)/s2bhat^2 - ((2*m2bhat - (2*(W2b^lam2 - 1))/lam2)*((2*W2b^lam2 - 2)/lam2^3 - (2*W2b^lam2*log(W2b))/lam2^2 + (W2b^lam2*log(W2b)^2)/lam2))/s2bhat^2 + (2*covb*(m1bhat - (W1b^lam1 - 1)/lam1)*((2*W2b^lam2 - 2)/lam2^3 - (2*W2b^lam2*log(W2b))/lam2^2 + (W2b^lam2*log(W2b)^2)/lam2))/(s1bhat^2*s2bhat^2))/((2*covb^2)/(s1bhat^2*s2bhat^2) - 2))
      
      
      
      
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
          1-pnorm(qnorm(1-t,mean=mean(W1alam),sd=std(W1alam)),mean=mean(W1blam),sd=std(W1blam))
        }
      
        roc2<-function(t){
          1-pnorm(qnorm(1-t,mean=mean(W2alam),sd=std(W2alam)),mean=mean(W2blam),sd=std(W2blam))
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
        return(list(resultstable=res,J1=J1original,J2=J2original, pvalue_transJ_difference= pval2t, CI_transJ_difference= CIZstar, pvalue_difference= pval2tJ, CI_difference= CIoriginal, roc1=roc1, roc2=roc2, transx1=W1alam, transy1=W1blam, transx2=W2alam, transy2=W2blam))

}


}



## ---- echo=FALSE, eval=TRUE---------------------------------------------------

comparebcSens <-function (marker1, marker2, D, atSpec, alpha, plots){

  if ((length(marker1) != length(D)) | (length(marker2) != length(D))) {
    stop("ERROR: The length of the 'marker' and 'D' inputs must be equal.")
  } else if (min(D) != 0 | max(D) != 1) {
    stop("ERROR: Controls must be assigned a value of 0; cases must be assigned a value of 1. Both controls and cases should be included in the dataset.")
  } else if (alpha <= 0 | alpha >= 1) {
    stop("ERROR: The level of significance, alpha, should be set between 0 and 1. A common choice is 0.05.")
  } else if (sum(is.na(marker1)) > 0 | sum(is.na(marker2)) > 0 | sum(is.na(D)) > 0 | sum(is.na(atSpec)) > 0) {
    stop("ERROR: Please remove all missing data before running this function.")
  } else if (!is.numeric(atSpec)) {
    stop("ERROR: 'atSpec' must be a numeric vector.")
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

    cova= cov(Walam)*(na-1)/na;
    covb= cov(Wblam)*(nb-1)/nb;
    mualam=c(mean(W1alam), mean(W2alam))
    mublam=c(mean(W1blam), mean(W2blam))

    out=  -sum(log(dmvnorm(Walam,mualam,cova)))-sum(log(dmvnorm(Wblam,mublam,covb))) -(lam[1]-1)*sum(log(W1a))-(lam[2]-1)*sum(log(W2a))-(lam[1]-1)*sum(log(W1b))-(lam[2]-1)*sum(log(W2b));

    return(out)
  }

  lam=c(2,2)
  likbox2D(W1a,W1b,W2a,W2b,lam)

  logL<-function(h){
    likbox2D(W1a,W1b,W2a,W2b,h)
  }

  init=c(0.8,0.8)
  lam=fminsearch(logL,init)
  lam=c(lam$optbase$xopt)
  lam

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
  h11_11=sum(((2*((W1a^lam1 - 1)/lam1^2 - (W1a^lam1*log(W1a))/lam1)^2)/s1ahat^2 - ((2*m1ahat - (2*(W1a^lam1 - 1))/lam1)*((2*W1a^lam1 - 2)/lam1^3 - (2*W1a^lam1*log(W1a))/lam1^2 + (W1a^lam1*log(W1a)^2)/lam1))/s1ahat^2 + (2*cova*(m2ahat - (W2a^lam2 - 1)/lam2)*((2*W1a^lam1 - 2)/lam1^3 - (2*W1a^lam1*log(W1a))/lam1^2 + (W1a^lam1*log(W1a)^2)/lam1))/(s1ahat^2*s2ahat^2))/((2*cova^2)/(s1ahat^2*s2ahat^2) - 2) + ((2*((W1b^lam1 - 1)/lam1^2 - (W1b^lam1*log(W1b))/lam1)^2)/s1bhat^2 - ((2*m1bhat - (2*(W1b^lam1 - 1))/lam1)*((2*W1b^lam1 - 2)/lam1^3 - (2*W1b^lam1*log(W1b))/lam1^2 + (W1b^lam1*log(W1b)^2)/lam1))/s1bhat^2 + (2*covb*(m2bhat - (W2b^lam2 - 1)/lam2)*((2*W1b^lam1 - 2)/lam1^3 - (2*W1b^lam1*log(W1b))/lam1^2 + (W1b^lam1*log(W1b)^2)/lam1))/(s1bhat^2*s2bhat^2))/((2*covb^2)/(s1bhat^2*s2bhat^2) - 2))
  h11_12=sum((2*cova*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1))/(lam1^2*lam2^2*(- 2*cova^2 + 2*s1ahat^2*s2ahat^2)) + (2*covb*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1))/(lam1^2*lam2^2*(- 2*covb^2 + 2*s1bhat^2*s2bhat^2)))


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
  h12_11=sum((2*cova*(W1a^lam1*lam1*log(W1a) - W1a^lam1 + 1)*(W2a^lam2*lam2*log(W2a) - W2a^lam2 + 1))/(lam1^2*lam2^2*(- 2*cova^2 + 2*s1ahat^2*s2ahat^2)) + (2*covb*(W1b^lam1*lam1*log(W1b) - W1b^lam1 + 1)*(W2b^lam2*lam2*log(W2b) - W2b^lam2 + 1))/(lam1^2*lam2^2*(- 2*covb^2 + 2*s1bhat^2*s2bhat^2)))
  h12_12=sum(((2*((W2a^lam2 - 1)/lam2^2 - (W2a^lam2*log(W2a))/lam2)^2)/s2ahat^2 - ((2*m2ahat - (2*(W2a^lam2 - 1))/lam2)*((2*W2a^lam2 - 2)/lam2^3 - (2*W2a^lam2*log(W2a))/lam2^2 + (W2a^lam2*log(W2a)^2)/lam2))/s2ahat^2 + (2*cova*(m1ahat - (W1a^lam1 - 1)/lam1)*((2*W2a^lam2 - 2)/lam2^3 - (2*W2a^lam2*log(W2a))/lam2^2 + (W2a^lam2*log(W2a)^2)/lam2))/(s1ahat^2*s2ahat^2))/((2*cova^2)/(s1ahat^2*s2ahat^2) - 2) + ((2*((W2b^lam2 - 1)/lam2^2 - (W2b^lam2*log(W2b))/lam2)^2)/s2bhat^2 - ((2*m2bhat - (2*(W2b^lam2 - 1))/lam2)*((2*W2b^lam2 - 2)/lam2^3 - (2*W2b^lam2*log(W2b))/lam2^2 + (W2b^lam2*log(W2b)^2)/lam2))/s2bhat^2 + (2*covb*(m1bhat - (W1b^lam1 - 1)/lam1)*((2*W2b^lam2 - 2)/lam2^3 - (2*W2b^lam2*log(W2b))/lam2^2 + (W2b^lam2*log(W2b)^2)/lam2))/(s1bhat^2*s2bhat^2))/((2*covb^2)/(s1bhat^2*s2bhat^2) - 2))




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


  att=1-atSpec;
  t=att;
  s1alam=sqrt( var(W1alam)*(na-1)/na )
  s1blam=sqrt( var(W1blam)*(nb-1)/nb )
  s2alam=sqrt( var(W2alam)*(na-1)/na )
  s2blam=sqrt( var(W2blam)*(nb-1)/nb )

  m1alam=mean(W1alam)
  m1blam=mean(W1blam)
  m2alam=mean(W2alam)
  m2blam=mean(W2blam)

  #############################################################################
  #############################################################################
  #==================AT A GIVEN T FOR THE ORIGINAL===================

  ROC1or= qnorm(1-pnorm(qnorm(1-t,m1alam,s1alam),m1blam,s1blam));
  ROC2or= qnorm(1-pnorm(qnorm(1-t,m2alam,s2alam),m2blam,s2blam));


  dm1a =    -(2^(1/2)*exp(-((m1ahat - m1bhat)/s1bhat + (2^(1/2)*s1ahat*erfcinv(2*t))/s1bhat)^2/2))/(2*s1bhat*pi^(1/2))
  dm1b =    (2^(1/2)*exp(-((m1ahat - m1bhat)/s1bhat + (2^(1/2)*s1ahat*erfcinv(2*t))/s1bhat)^2/2))/(2*s1bhat*pi^(1/2))
  ds1a =    -(erfcinv(2*t)*exp(-((m1ahat - m1bhat)/s1bhat + (2^(1/2)*s1ahat*erfcinv(2*t))/s1bhat)^2/2))/(s1bhat*pi^(1/2))
  ds1b =    (2^(1/2)*exp(-((m1ahat - m1bhat)/s1bhat + (2^(1/2)*s1ahat*erfcinv(2*t))/s1bhat)^2/2)*((m1ahat - m1bhat)/s1bhat^2 + (2^(1/2)*s1ahat*erfcinv(2*t))/s1bhat^2))/(2*pi^(1/2))
  der1=c(dm1a, ds1a, dm1b, ds1b);
  vROC1=c(der1)%*%VV[1:4,1:4]%*%t(t(c(der1)))
  vROC1=c(vROC1)

  dm2a =    -(2^(1/2)*exp(-((m2ahat - m2bhat)/s2bhat + (2^(1/2)*s2ahat*erfcinv(2*t))/s2bhat)^2/2))/(2*s2bhat*pi^(1/2))
  dm2b =    (2^(1/2)*exp(-((m2ahat - m2bhat)/s2bhat + (2^(1/2)*s2ahat*erfcinv(2*t))/s2bhat)^2/2))/(2*s2bhat*pi^(1/2))
  ds2a =    -(erfcinv(2*t)*exp(-((m2ahat - m2bhat)/s2bhat + (2^(1/2)*s2ahat*erfcinv(2*t))/s2bhat)^2/2))/(s2bhat*pi^(1/2))
  ds2b =    (2^(1/2)*exp(-((m2ahat - m2bhat)/s2bhat + (2^(1/2)*s2ahat*erfcinv(2*t))/s2bhat)^2/2)*((m2ahat - m2bhat)/s2bhat^2 + (2^(1/2)*s2ahat*erfcinv(2*t))/s2bhat^2))/(2*pi^(1/2))
  der2=c(dm2a, ds2a, dm2b, ds2b);

  vROC2=c(der2)%*%VV[5:8,5:8]%*%t(t(c(der2)))
  covROC=c(der1, 0, 0, 0, 0)%*%VV[1:8,1:8]%*%t(t(c(0, 0, 0, 0, der2)))

  Z=(ROC2or-ROC1or)/(sqrt(vROC2+vROC1-2*covROC))
  SEZ=sqrt(vROC1+vROC2-2*covROC)

  pval2tZ=2*pnorm(-abs(Z));

  CIoriginal=c((ROC2or-ROC1or)-pmalpha*SEZ, (ROC2or-ROC1or)+pmalpha*SEZ)
  CIoriginal

  #############################################################################
  #############################################################################
  #==================AT A GIVEN T FOR THE TRANSORMED===================


  ROC1MC= qnorm(1-pnorm(qnorm(1-t,m1alam,s1alam),m1blam,s1blam));
  ROC2MC= qnorm(1-pnorm(qnorm(1-t,m2alam,s2alam),m2blam,s2blam));


  dm1a =-1/s1blam;
  ds1a =-(2^(1/2)*erfcinv(2*t))/s1blam;
  dm1b =1/s1blam;
  ds1b =(m1alam - m1blam)/s1blam^2 + (2^(1/2)*s1alam*erfcinv(2*t))/s1blam^2;
  der1=c(dm1a, ds1a, dm1b, ds1b);
  vROC1=c(der1)%*%VV[1:4,1:4]%*%t(t(c(der1)))
  vROC1=c(vROC1)


  dm2a =-1/s2blam;
  ds2a =-(2^(1/2)*erfcinv(2*t))/s2blam;
  dm2b =1/s2blam;
  ds2b =(m2alam - m2blam)/s2blam^2 + (2^(1/2)*s2alam*erfcinv(2*t))/s2blam^2;
  der2=c(dm2a, ds2a, dm2b, ds2b);
  vROC2=c(der2)%*%VV[5:8,5:8]%*%t(t(c(der2)))
  vROC2=c(vROC2)


  covROC=c(der1, 0, 0, 0, 0)%*%VV[1:8,1:8]%*%t(t(c(0, 0, 0, 0, der2)))
  covROC=c(covROC)


  Zroct=(ROC2MC-ROC1MC)/(sqrt(vROC2+vROC1-2*covROC))
  SEroct=sqrt(vROC1+vROC2-2*covROC)

  pval2t_roct=2*pnorm(-abs(Zroct));pval2t=pval2t_roct;
  CIZstar=c((ROC2MC-ROC1MC)-pmalpha*SEroct, (ROC2MC-ROC1MC)+pmalpha*SEroct)
  CIZstar

  SE1=pnorm(ROC1MC)
  SE2=pnorm(ROC2MC)

  res <- matrix(c(SE1,SE2, pval2t, CIZstar[1],CIZstar[2]),ncol=5,byrow=TRUE)
  colnames(res) <- c("Sens 1:","    Sens 2:","    p-value (probit):","   CI probit (LL):","    CI probit (UL):")
  rownames(res) <- c("Estimates:")
  res <- as.table(res)
  res


  #====TWO ROC FUNCTIONS: ONE IS THE BOXCOX AND THE OTHER THE REFERENCE LINE================
  roc1<-function(t){
    1-pnorm(qnorm(1-t,mean=mean(W1alam),sd=std(W1alam)),mean=mean(W1blam),sd=std(W1blam))
  }

  roc2<-function(t){
    1-pnorm(qnorm(1-t,mean=mean(W2alam),sd=std(W2alam)),mean=mean(W2blam),sd=std(W2blam))
  }

  rocuseless<-function(t){
    1-pnorm(qnorm(1-t,mean=1,sd=1),mean=1,sd=1)
  }


  #================IF PLOTS ARE REQUESTED PLOT THE ROCS==
  if (plots=="on") {
    # x11()
    plot(linspace(0,1,1000),roc1(linspace(0,1,1000)),main="Box-Cox Based ROCs",xlab="FPR = 1 - Specificity",ylab="TPR = Sensitivity",type="l",col="red")
    lines(linspace(0,1,1000),roc2(linspace(0,1,1000)),main="Box-Cox Based ROCs",xlab="FPR = 1 - Specificity",ylab="TPR = Sensitivity",type="l",col="black")
    lines(linspace(0,1,1000),linspace(0,1,1000),type="l", lty=2)


    points(1-atSpec,SE1, col = "red")
    points(1-atSpec,SE2, col = "black")


    legend("bottomright", legend=c(paste("ROC for Marker 1 with Sens =", round(SE1,4), "at Spec =", atSpec),
                                   paste("ROC for Marker 2 with Sens =", round(SE2,4), "at Spec =", atSpec),
                                   paste("P-value for the Sens difference:",round(pval2t,4))),

           col=c("red", "black", "white"), lty=c(1, 1, NA), pch = c(NA, NA, NA), cex=0.8)

  }




  return(list(resultstable=res,Sens1=SE1,Sens2=SE2, pvalue_probit_difference= pval2t, CI_probit_difference= CIZstar, pvalue_difference= pval2tZ, CI_difference= CIoriginal, roc1=roc1, roc2=roc2, transx1=W1alam, transy1=W1blam, transx2=W2alam, transy2=W2blam))




  }
}


## ---- echo=FALSE, eval=TRUE---------------------------------------------------
comparebcSpec <-function (marker1, marker2, D, atSens, alpha, plots){
  
    if ((length(marker1) != length(D)) | (length(marker2) != length(D))) {
    stop("ERROR: The length of the 'marker' and 'D' inputs must be equal.")
  } else if (min(D) != 0 | max(D) != 1) {
    stop("ERROR: Controls must be assigned a value of 0; cases must be assigned a value of 1. Both controls and cases should be included in the dataset.")
  } else if (alpha <= 0 | alpha >= 1) {
    stop("ERROR: The level of significance, alpha, should be set between 0 and 1. A common choice is 0.05.")
  } else if (sum(is.na(marker1)) > 0 | sum(is.na(marker2)) > 0 | sum(is.na(D)) > 0 | sum(is.na(atSens)) > 0) {
    stop("ERROR: Please remove all missing data before running this function.")
  } else if (!is.numeric(atSens)) {
    stop("ERROR: 'atSpec' must be a numeric vector.")
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
  
  Dflip=abs(D-1);
  atSpec=1-atSens;
  #out=comparebcSens(marker1group1=W1b, marker1group2=W1a ,marker2group1=W2b, marker2group2=W2a, atSpec=atSens, plots="off")
  out=comparebcSens(marker1, marker2, Dflip, atSpec=1-atSens, alpha, plots="off")
  #return(list(resultstable=res,Sens1=SE1,Sens2=SE2, pvalue_probit_difference= pval2t, CI_probit_difference= CIZstar, pvalue_difference= pval2tZ, CI_difference= CIoriginal, roc1=roc1, roc2=roc2, transx1=W1alam, transy1=W1blam, transx2=W2alam, transy2=W2blam))
  
FPR1=out$Sens1
FPR2=out$Sens2

CIZstar=out$CI_probit_difference
CIoriginal=out$CI_difference

pval2t=out$pvalue_probit_difference
pval2tZ=out$pvalue_difference



W1alam=out$transx1
W1blam=out$transy1
W2alam=out$transx2
W2blam=out$transy2
#====TWO ROC FUNCTIONS: ONE IS THE BOXCOX AND THE OTHER THE REFERENCE LINE================
roct1<-function(t){
  1-pnorm(qnorm(1-t,mean=mean(W1blam),sd=std(W1blam)),mean=mean(W1alam),sd=std(W1alam))
}

roct2<-function(t){
  1-pnorm(qnorm(1-t,mean=mean(W2blam),sd=std(W2blam)),mean=mean(W2alam),sd=std(W2alam))
}

rocuseless<-function(t){
  1-pnorm(qnorm(1-t,mean=1,sd=1),mean=1,sd=1)
}


#================IF PLOTS ARE REQUESTED PLOT THE ROCS==
if (plots=="on") {
  # x11()
  plot(linspace(0,1,1000),roct1(linspace(0,1,1000)),main="Box-Cox Based ROCs",xlab="FPR = 1 - Specificity",ylab="TPR = Sensitivity",type="l",col="red")
  lines(linspace(0,1,1000),roct2(linspace(0,1,1000)),main="Box-Cox Based ROCs",xlab="FPR = 1 - Specificity",ylab="TPR = Sensitivity",type="l",col="black")
  lines(linspace(0,1,1000),linspace(0,1,1000),type="l", lty=2)
  
  
  points(FPR1, 1-atSpec, col = "red")
  points(FPR2, 1-atSpec, col = "black")
  
  #line for the Youden:
  #lines(c(1-atSpec,1-SE1),c(1-atSpec,1-SE2),col="green")
  
  legend("bottomright", legend=c(paste("ROC for Marker 1 with FPR =", round(FPR1,4), "at Sens =", atSens), 
                                 paste("ROC for Marker 2 with FPR =", round(FPR2,4), "at Sens =", atSens),
                                 paste("P-value for the Spec difference:",round(pval2t,4))),
         
         col=c("red", "black", "white"), lty=c(1, 1, NA), pch = c(NA, NA, NA), cex=0.8)
  
}


res <- matrix(c(FPR1,FPR2, pval2t, CIZstar[1],CIZstar[2]),ncol=5,byrow=TRUE)
colnames(res) <- c("FPR 1:","    FPR 2:","    p-value (probit):","   CI probit (LL):","    CI probit (UL):")
rownames(res) <- c("Estimates:")
res <- as.table(res)
res


#return(list(resultstable=res,Sens1=SE1,Sens2=SE2, pvalue_probit_difference= pval2t, CI_probit_difference= CIZstar, pvalue_difference= pval2tZ, CI_difference= CIoriginal, roc1=roc1, roc2=roc2, transx1=W1alam, transy1=W1blam, transx2=W2alam, transy2=W2blam))
 return(list(resultstable=res,FPR1=FPR1,FPR2=FPR2, pvalue_probit_difference= pval2t, CI_probit_difference= CIZstar, pvalue_difference= pval2tZ, CI_difference= CIoriginal, roc1=roct1, roc2=roct2, transx1=W1alam, transy1=W1blam, transx2=W2alam, transy2=W2blam))
#return(list(FPR1=FPR1))


  }
}



## ---- echo=TRUE, eval=TRUE----------------------------------------------------
#DATA GENERATION
set.seed(123)
x=rgamma(100, shape=2, rate = 8) # Generates biomarker data from a gamma distribution 
                                 # for the healthy group.
y=rgamma(100, shape=2, rate = 4) # Generates biomarker data from a gamma distribution 
                                 # for the diseased group.
scores=c(x,y)
D=c(zeros(1,100), ones(1,100))

out=checkboxcox(marker=scores, D, plots="on", printShapiro = TRUE)

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
summary(out)

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
set.seed(123)
x=rgamma(100, shape=2, rate = 8)
y=rgamma(100, shape=2, rate = 4)
scores=c(x,y)
D=c(zeros(1,100), ones(1,100))
out=rocboxcox(marker=scores,D, 0.05, plots="on", printProgress=FALSE)

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
givenSP=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
givenSE=NA
out=rocboxcoxCI(marker=scores ,D, givenSP=givenSP, givenSE=NA, alpha=0.05, plots="on")
#

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
#GENERATE SOME BIVARIATE DATA===
set.seed(123)

nx <- 100
 Sx <- matrix(c(1,   0.5,
                0.5,  1), 
            nrow = 2, ncol = 2)

mux <- c(X = 10, Y = 12)
X=rmvnorm(nx, mean = mux, sigma = Sx)

ny <- 100
Sy <- matrix(c(1.1,   0.6,
               0.6,  1.1), 
             nrow = 2, ncol = 2)

muy <- c(X = 11, Y = 13.7)
Y=rmvnorm(ny, mean = muy, sigma = Sy)

dx=zeros(nx,1)
dy=ones(ny,1)

markers=rbind(X,Y);
marker1=markers[,1]
marker2=markers[,2]
D=c(rbind(dx,dy))

#==DATA HAVE BEEN GENERATED====


#===COMPARE THE AUCs of Marker 1 vs Maker 2

out=comparebcAUC(marker1, marker2, D, alpha=0.05,  plots="on")



## ---- echo=TRUE, eval=TRUE----------------------------------------------------


#==DATA HAVE BEEN GENERATED====
#===COMPARE THE Js of Marker 1 vs Maker 2

out=comparebcJ(marker1, marker2, D, alpha=0.05,  plots="on")

#


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

out=comparebcSens(marker1=marker1, marker2=marker2, D=D, alpha =0.05, atSpec=0.8, plots="on")
summary(out)
#  


## ---- echo=TRUE, eval=TRUE----------------------------------------------------

out=comparebcSpec(marker1=marker1, marker2=marker2, D=D, alpha =0.05, atSens=0.8, plots="on")
summary(out)
#  


