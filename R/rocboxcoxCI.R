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
  } else if (!((sum(is.na(givenSP)) == 0 & sum(is.na(givenSE) > 0)) |
               (sum(is.na(givenSE)) == 0 & sum(is.na(givenSP) > 0)))) {
    stop("ERROR: Exactly one of 'givenSP' and 'givenSE' must be set to NA.")
  } else if (sum(marker < 0) > 0) {
    stop("ERROR: To use the Box-Cox transformation, all marker values must be positive.")
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




    roc<-function(t,trans_x,trans_y){

      na = length(trans_x)
      nb = length(trans_y)

      stransx=sqrt(var(trans_x)*(na-1)/na)
      stransy=sqrt(var(trans_y)*(nb-1)/nb)

      1-pnorm(qnorm(1-t,
                    mean=mean(trans_x),
                    sd=stransx),
              mean=mean(trans_y),
              sd=stransy)
    }

    rocuseless<-function(t){
      1-pnorm(qnorm(1-t,mean=1,sd=1),mean=1,sd=1)
    }

    # x11()
    #plot.new()

    if (plots == "on") {
      txt <- paste("Box-Cox Based ROC, alpha =", round(alpha,3) )
      plot(linspace(0,1,1000),roc(linspace(0,1,1000),transx,transy),main=" ",xlab="FPR",ylab="TPR",type="l",col="red")
      lines(linspace(0,1,10),linspace(0,1,10),type="l", lty=2)
      title(main=txt)
    }

    m1hat=mean(transx)
    m2hat=mean(transy)
    # s1hat=std(transx)
    # s2hat=std(transy)
    n1=length(x)
    n2=length(y)



    I=zeros(5,5);
    sh=sqrt(1/(length(transx))*sum((transx-mean(transx))^2))
    sd=sqrt(1/(length(transy))*sum((transy-mean(transy))^2))

    s1hat = sh
    s2hat = sd

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
    SEandCIs = matrix( c(Spvalues, Sehat, CIROCvec1,CIROCvec2), nrow=length(CIROCvec1), ncol=4)
    colnames(SEandCIs)  <- c("Given Sp", "Sehat","LL of 95%CI","UL of 95%CI")
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

    return(list(SPandCIs=formattable(as.matrix(SPandCIs), digits = 4, format = "f"),
                SEandCIs=formattable(as.matrix(SEandCIs), digits = 4, format = "f"),
                Sevalues=Sevalues,
                Sphat=Sphat,
                CIsp=formattable(as.matrix(CIsp), digits = 4, format = "f"),
                Spvalues=Spvalues,
                Sehat=Sehat,
                CIse=formattable(as.matrix(CIse), digits = 4, format = "f")))


  }

}
