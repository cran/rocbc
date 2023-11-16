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
  } else if (min(marker<= 0)) {
    stop("ERROR: All marker scores need to be positive")
  } else if (sum(marker < 0) > 0) {
    stop("ERROR: To use the Box-Cox transformation, all marker values must be positive.")
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

    COVlam <-function (x,y){

      #%Accepts X<Y Yields as an output the Variance Covariance matrix
      #%of the mh, sh, md, sd, and lam, (with this order) where lam is the estimated boxcox
      #%parameter Thus, the output will be a 5x5 matrix

      #%---FOR TWO POPULATIONS:

      n=length(x); #original scores
      m=length(y);

      #  %mh-->g(1)
      #  %sh-->g(2)
      #  %md-->g(3)
      #  %sd-->g(4)
      #  %lam-->g(5)
      #  %  logL=@(g) -n*log(sh)-1/(2*sh^2)*sum((xlam-mh)^2)+(lam-1)*sum(log(x)) -
      #  %             m*log(sd)-1/(2*sd^2)*sum((ylam-md)^2)+(lam-1)*sum(log(y))

      outt=boxcoxleo(x,y); #the original scores
      transx=outt$transx
      transy=outt$transy
      lam=outt$transformation.parameter


      mh=mean(transx); xlam=transx;
      md=mean(transy); ylam=transy;
      sh=sqrt(1/(length(transx))*sum((transx-mean(transx))^2));
      sd=sqrt(1/(length(transy))*sum((transy-mean(transy))^2));


      I=zeros(5,5);
      I[1,1]=n/sh^2;
      I[2,2]=-(n/sh^2-3/sh^4*sum((xlam-mh)^2));
      I[3,3]=m/sd^2;
      I[4,4]=-(m/sd^2-3/sd^4*sum((ylam-md)^2));

      kk=  sum(((mh - (x^lam - 1)/lam)*((2*(x^lam - 1))/lam^3 - (2*x^lam*log(x))/lam^2 + (x^lam*log(x)^2)/lam))/sh^2) +  - sum(((y^lam - 1)/lam^2 - (y^lam*log(y))/lam)^2/sd^2) +   - sum(((x^lam - 1)/lam^2 - (x^lam*log(x))/lam)^2/sh^2) +   + sum(((md - (y^lam - 1)/lam)*((2*(y^lam - 1))/lam^3 - (2*y^lam*log(y))/lam^2 + (y^lam*log(y)^2)/lam))/sd^2);
      I[5,5]=-kk ;

      I[1,5]=sum(((2*(x^lam - 1))/lam^2 - (2*x^lam*log(x))/lam)/(2*sh^2));
      I[2,5]=-sum((2*((x^lam - 1)/lam^2 - (x^lam*log(x))/lam)*(mh - (x^lam - 1)/lam))/sh^3);
      I[3,5]=-sum(-((2*(y^lam - 1))/lam^2 - (2*y^lam*log(y))/lam)/(2*sd^2));


      I[4,5]=-sum((2*((y^lam - 1)/lam^2 - (y^lam*log(y))/lam)*(md - (y^lam - 1)/lam))/sd^3);

      I[5,1]=I[1,5];
      I[5,2]=I[2,5];
      I[5,3]=I[3,5];
      I[5,4]=I[4,5];


      S= inv(I); #%that's better

      return(list(out=S))

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

    if (plots=="on"){
      txt <- paste("Box-Cox Based ROC, alpha =", formattable(alpha, digits = 3))
      plot(linspace(0,1,1000),roc(linspace(0,1,1000),transx,transy),main=" ",xlab="FPR = 1 - Specificity",ylab="TPR = Sensitivity",type="l",col="red")
      lines(linspace(0,1,10),linspace(0,1,10),type="l", lty=2)
      title(main=txt)

    }
    rocfun=function(t){1-pnorm(qnorm(1-t,mean=mean(transx),sd=sqrt(var(transx)*(length(transx)-1)/length(transx))),mean=mean(transy),sd=sqrt(var(transy)*(length(transy)-1)/length(transy)))}
    m1hat=mean(transx)
    m2hat=mean(transy)
    s1hat=sqrt(var(transx)*(length(transx)-1)/length(transx))
    s2hat=sqrt(var(transy)*(length(transy)-1)/length(transy))
    n1=length(x)
    n2=length(y)

    #NOTE THAT:  m1hat->m1
    #            m2hat->m2
    #            s1hat->s1hat
    #            s2hat->s2hat

    #########################################

    derde =c( (s2hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s2hat*(s1hat^2 - s2hat^2)), ((s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - 2*m2hat*s1hat + (s1hat*s2hat*(2*s1hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s1hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) - (2*s1hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s2hat, -((s1hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) - 1)/s2hat, ((2*m1hat*s2hat + s1hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - (s1hat*s2hat*(2*s2hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) + (2*s2hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s2hat - (m2hat + (m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2))/(s1hat^2 - s2hat^2))/s2hat^2);
    derdp =c( -((s2hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) + 1)/s1hat, (m1hat + (m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2))/(s1hat^2 - s2hat^2))/s1hat^2 - ((s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - 2*m2hat*s1hat + (s1hat*s2hat*(2*s1hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s1hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) - (2*s1hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s1hat, (s1hat^2 + (s1hat*s2hat*(2*m1hat - 2*m2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat*(s1hat^2 - s2hat^2)), -((2*m1hat*s2hat + s1hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2) - (s1hat*s2hat*(2*s2hat*log(s1hat^2/s2hat^2) + (2*(s1hat^2 - s2hat^2))/s2hat))/(2*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2) + (2*s2hat*(m1hat*s2hat^2 - m2hat*s1hat^2 + s1hat*s2hat*(log(s1hat^2/s2hat^2)*(s1hat^2 - s2hat^2) + (m1hat - m2hat)^2)^(1/2)))/(s1hat^2 - s2hat^2)^2)/s1hat);

    outt=COVlam(xor,yor);

    S=outt$out;
    S=S[1:4,1:4];

    if (sum(is.infinite(S)) > 0 |
        sum(is.nan(S)) > 0 |
        sum(is.na(S)) > 0 |
        sum(diag(S) < 0) > 0) {
      stop("ERROR: The information matrix cannot be inverted - this might be due to the scale of the marker. Try rescaling the marker measurements by subtracting a constant or by dividing with a constant all marker measurements (of both groups) before you try again.")
    }

    varde=t(derde)%*%S%*%t(t(derde));
    vardp=t(derdp)%*%S%*%t(t(derdp));
    covdedp=t(derde)%*%S%*%t(t(derdp));
    ##########################################

    vardeltase=varde;
    vardeltasp=vardp;

    #########################

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

    if (sum(is.infinite(covdeltas)) > 0 |
        sum(is.nan(covdeltas)) > 0) {
      stop("ERROR: The information matrix cannot be inverted - this might be due to the scale of the marker. Try rescaling the marker measurements by subtracting a constant or by dividing with a constant all marker measurements (of both groups) before you try again.")
    }

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

    ## EVERYTHING ELSE ================================================================================

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

    b=sy/sx;
    a=my-mx;

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
    #CI[1]=exp(CIstar[1])/(1+exp(CIstar[1]))*2-1
    #CI[2]=exp(CIstar[2])/(1+exp(CIstar[2]))*2-1

    CI[1]=pnorm(CIstar[1])
    CI[2]=pnorm(CIstar[2])

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
      s1hat_boots=sqrt(var(transx1_boots)*(length(transx1_boots)-1)/length(transx1_boots))
      s2hat_boots=sqrt(var(transx2_boots)*(length(transx2_boots)-1)/length(transx2_boots))

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
                                           formattable(auc, digits = 4, format = "f"),
                                           ", CI: (",
                                           formattable(CIauc[1], digits = 4, format = "f"),
                                           ", ",
                                           formattable(CIauc[2], digits = 4, format = "f") ,")", sep = ""),
                                     paste("Maximized Youden index = ", formattable(Jhat, digits = 4, format = "f"),
                                           ", CI: (",
                                           formattable((CIJ[1]), digits = 4, format = "f"),
                                           ", ",
                                           formattable((CIJ[2]), digits = 4, format = "f"), ")", sep = ""),
                                     paste("Optimal pair of (FPR,TPR): (",
                                           formattable(1-Sp, digits = 4, format = "f"),
                                           ", ",
                                           formattable(Se, digits = 4, format = "f"), ")", sep = ""),
                                     paste("Area of the rect. conf. region =", formattable(arearect, digits = 4, format = "f")),
                                     paste("Area of the egg conf. region =", formattable(areaegg, digits = 4, format = "f")),
                                     paste("Youden based cutoff: ", formattable(cutoff, digits = 4, format = "f"),", CI: (",
                                           formattable(CIcutoff[1], digits = 4, format = "f"),
                                           ", ",
                                           formattable(CIcutoff[2], digits = 4, format = "f"), ")", sep = ""),
                                     paste("Sp Marginal CI: (",
                                           c(formattable(margcisp[1], digits = 4, format = "f")), ", ", formattable(margcisp[2], digits = 4, format = "f"),
                                           ")", sep = ""),
                                     paste("Se Marginal CI: (",
                                           c(formattable(margcise[1], digits = 4, format = "f")), ", ", formattable(margcise[2], digits = 4, format = "f"), ")",
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
                    Se,margcise[1],margcise[2],
                    lam,NA,NA),ncol=3,byrow=TRUE)
    colnames(res) <- c("Estimate","Confidence","Interval")
    rownames(res) <- c("AUC","Jhat","cutoff","Sp","Se", "Lambda")
    res <- data.frame(res)
    res <- formattable(as.matrix(res), digits = 4, format = "f")
    res

    if (auc < 0.5) {
      print("NOTE: AUC < 0.5; the ordering of the two groups may need to be reversed.")
    }

    return(list(transx=((x^lam)-1)/lam, transy=((y^lam)-1)/lam , transformation.parameter=lam, AUC=auc, AUCCI=CIauc, pvalueAUC=pvalauc, J=Jhat, JCI=CIJ, pvalueJ=pvalJ, Sens=Se, CImarginalSens=margcise, Spec=Sp, CImarginalSpec=margcisp, cutoff=cutoff, CIcutoff=CIcutoff,  areaegg=areaegg, arearect=arearect, mxlam=mean(transx), sxlam=sqrt(var(transx)*(length(transx)-1)/length(transx)), mylam=mean(transy), sylam=sqrt(var(transy)*(length(transy)-1)/length(transy)), results=res , roc=rocfun))
  }
}
