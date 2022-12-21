comparebcSpec <-function (marker1, marker2, D, atSens, alpha, plots){

  if ((length(marker1) != length(D)) | (length(marker2) != length(D))) {
    stop("ERROR: The length of the 'marker' and 'D' inputs must be equal.")
  } else if (min(D) != 0 | max(D) != 1) {
    stop("ERROR: Controls must be assigned a value of 0; cases must be assigned a value of 1. Both controls and cases should be included in the dataset.")
  } else if (alpha <= 0 | alpha >= 1) {
    stop("ERROR: The level of significance, alpha, should be set between 0 and 1. A common choice is 0.05.")
  } else if (sum(is.na(marker1)) > 0 | sum(is.na(marker2)) > 0 | sum(is.na(D)) > 0 | sum(is.na(atSens)) > 0) {
    stop("ERROR: Please remove all missing data before running this function.")
  } else if ((!is.numeric(atSens)) | (length(atSens) != 1)) {
    stop("ERROR: 'atSens' must be a single numeric value.")
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
      na = length(W1alam)
      nb = length(W1blam)

      s1alam=sqrt(var(W1alam)*(na-1)/na)
      s1blam=sqrt(var(W1blam)*(nb-1)/nb)

      1-pnorm(qnorm(1-t,
                    mean=mean(W1blam),
                    sd=s1blam),
              mean=mean(W1alam),
              sd=s1alam)
    }

    roct2<-function(t){
      na = length(W2alam)
      nb = length(W2blam)

      s2alam=sqrt(var(W2alam)*(na-1)/na)
      s2blam=sqrt(var(W2blam)*(nb-1)/nb)

      1-pnorm(qnorm(1-t,
                    mean=mean(W2blam),
                    sd=s2blam),
              mean=mean(W2alam),
              sd=s2alam)
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

      legend("bottomright", legend=c(paste("ROC for Marker 1 with FPR =", formattable(FPR1, digits = 4, format = "f"), "at Sens =", formattable(atSens, digits = 4, format = "f")),
                                     paste("ROC for Marker 2 with FPR =", formattable(FPR2, digits = 4, format = "f"), "at Sens =", formattable(atSens, digits = 4, format = "f")),
                                     paste("P-value for the Spec difference:", formattable(pval2t, digits = 4, format = "f"))),

             col=c("red", "black", "white"), lty=c(1, 1, NA), pch = c(NA, NA, NA), cex=0.8)

    }

    res <- data.frame(FPR1 = formattable(FPR1, digits = 4, format = "f"),
                      FPR2 = formattable(FPR2, digits = 4, format = "f"),
                      p_value_probit = formattable(pval2t, digits = 4, format = "f"),
                      p_value = formattable(pval2tZ, digits = 4, format = "f"),
                      ci_ll = formattable(CIoriginal[1], digits = 4, format = "f"),
                      ci_ul = formattable(CIoriginal[2], digits = 4, format = "f"))
    rownames(res) <- c("Estimates:")
    colnames(res) <- c("AUC 1", "AUC 2", "P-Val (Probit)", "P-Val", "CI (LL)", "CI (UL)")
    res <- formattable(as.matrix(res), digits = 4, format = "f")
    res


    #return(list(resultstable=res,Sens1=SE1,Sens2=SE2, pvalue_probit_difference= pval2t, CI_probit_difference= CIZstar, pvalue_difference= pval2tZ, CI_difference= CIoriginal, roc1=roc1, roc2=roc2, transx1=W1alam, transy1=W1blam, transx2=W2alam, transy2=W2blam))
    return(list(resultstable=res,FPR1=FPR1,FPR2=FPR2, pvalue_probit_difference= pval2t, pvalue_difference= pval2tZ, CI_difference= CIoriginal, roc1=roct1, roc2=roct2, transx1=W1alam, transy1=W1blam, transx2=W2alam, transy2=W2blam))
    #return(list(FPR1=FPR1))


  }
}
