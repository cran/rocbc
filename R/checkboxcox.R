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
