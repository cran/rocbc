checkboxcox2 <-function (marker1, marker2, D, plots, printShapiro = FALSE){

  if ((length(marker1) != length(D)) | (length(marker2) != length(D))) {
    stop("ERROR: The length of the 'marker' and 'D' inputs must be equal.")
  } else if (min(D) != 0 | max(D) != 1) {
    stop("ERROR: Controls must be assigned a value of 0; cases must be assigned a value of 1. Both controls and cases should be included in the dataset.")
  } else if ((sum(marker1 < 0) > 0) | (sum(marker2 < 2)) > 0) {
    stop("ERROR: To use the Box-Cox transformation, all marker values must be positive.")
  } else {

    erf <- function (x) 2 * pnorm(x * sqrt(2)) - 1
    erfc <- function (x) 2 * pnorm(x * sqrt(2), lower.tail = FALSE)
    erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
    erfcinv <- function (x) qnorm(x/2, lower.tail = FALSE)/sqrt(2)

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

    # init=c(0.8,0.8)
    # lam=fminsearch(logL,init)
    # lam=c(lam$optbase$xopt)
    lam<- optim(c(1,1),logL,gr=NULL,method="BFGS", control=list(maxit=10000))
    lam=c(lam$par)
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

    # ==== PERFORM SHAPIRO-WILK TESTS =============================================================

    x1 = W1a; y1 = W1b; transx1 = W1alam; transy1 = W1blam
    x2 = W2a; y2 = W2b; transx2 = W2alam; transy2 = W2blam

    if (!printShapiro) {
      test_x1=shapiro.test(x1)
      test_y1=shapiro.test(y1)
      test_t_x1=shapiro.test(transx1)
      test_t_y1=shapiro.test(transy1)

      test_x2=shapiro.test(x2)
      test_y2=shapiro.test(y2)
      test_t_x2=shapiro.test(transx2)
      test_t_y2=shapiro.test(transy2)
    }

    if (printShapiro) {
      test_x1=print(shapiro.test(x1))
      test_y1=print(shapiro.test(y1))
      test_t_x1=print(shapiro.test(transx1))
      test_t_y1=print(shapiro.test(transy1))

      test_x2=print(shapiro.test(x2))
      test_y2=print(shapiro.test(y2))
      test_t_x2=print(shapiro.test(transx2))
      test_t_y2=print(shapiro.test(transy2))
    }

    pval_x1=test_x1$p.value
    pval_y1=test_y1$p.value
    pval_t_x1=test_t_x1$p.value
    pval_t_y1=test_t_y1$p.value

    pval_x2=test_x2$p.value
    pval_y2=test_y2$p.value
    pval_t_x2=test_t_x2$p.value
    pval_t_y2=test_t_y2$p.value

    pvalues=c(pval_x1, pval_y1, pval_t_x1, pval_t_y1,
              pval_x2, pval_y2, pval_t_x2, pval_t_y2)

    res_shapiro = list(test_x1 = test_x1,
                       test_y1 = test_y1,
                       test_t_x1 = test_t_x1,
                       test_t_y1 = test_t_y1,
                       test_x2 = test_x2,
                       test_y2 = test_y2,
                       test_t_x2 = test_t_x2,
                       test_t_y2 = test_t_y2)

    #====TWO ROC FUNCTIONS: ONE IS THE BOXCOX AND THE OTHER THE REFERENCE LINE================

    roc1<-function(t){
      na = length(W1alam)
      nb = length(W1blam)

      s1alam=sqrt( var(W1alam)*(na-1)/na )
      s1blam=sqrt( var(W1blam)*(nb-1)/nb )

      1-pnorm(qnorm(1-t,
                    mean=mean(W1alam),
                    sd=s1alam),
              mean=mean(W1blam),
              sd=s1blam)
    }

    roc2<-function(t){
      na = length(W2alam)
      nb = length(W2blam)

      s2alam=sqrt( var(W2alam)*(na-1)/na )
      s2blam=sqrt( var(W2blam)*(nb-1)/nb )

      1-pnorm(qnorm(1-t,
                    mean=mean(W2alam),
                    sd=s2alam),
              mean=mean(W2blam),
              sd=s2blam)
    }

    rocuseless<-function(t){
      1-pnorm(qnorm(1-t,mean=1,sd=1),mean=1,sd=1)
    }

    #================IF PLOTS ARE REQUESTED PLOT ALL OF THE THINGS====

    if (plots=="on") {

      par(mfrow=c(2,2))
      par(mar = c(4.1, 3.1, 3.1, 1.1))

      # HISTOGRAMS & QQ, X1

      hist(W1a, main = "Histogram of X1", xlab = "X1")
      qqnorm(W1a, main = "Normal Q-Q Plot of X1")
      qqline(W1a, col = 2,lwd=2,lty=2)

      hist(W1alam, main = "Histogram of Transformed X1", xlab = "Transformed X1")
      qqnorm(W1alam, main = "Normal Q-Q Plot of Transformed X1")
      qqline(W1alam, col = 2,lwd=2,lty=2)

      # HISTOGRAMS & QQ, Y1

      hist(W1b, main = "Histogram of Y1", xlab = "Y1")
      qqnorm(W1b, main = "Normal Q-Q Plot of Y1")
      qqline(W1b, col = 2,lwd=2,lty=2)

      hist(W1blam, main = "Histogram of Transformed Y1", xlab = "Transformed Y1")
      qqnorm(W1blam, main = "Normal Q-Q Plot of Transformed Y1")
      qqline(W1blam, col = 2,lwd=2,lty=2)

      # HISTOGRAMS & QQ, X2

      hist(W2a, main = "Histogram of X2", xlab = "X2")
      qqnorm(W2a, main = "Normal Q-Q Plot of X2")
      qqline(W2a, col = 2,lwd=2,lty=2)

      hist(W2alam, main = "Histogram of Transformed X2", xlab = "Transformed X2")
      qqnorm(W2alam, main = "Normal Q-Q Plot of Transformed X2")
      qqline(W2alam, col = 2,lwd=2,lty=2)

      # HISTOGRAMS & QQ, Y2

      hist(W2b, main = "Histogram of Y2", xlab = "Y2")
      qqnorm(W2b, main = "Normal Q-Q Plot of Y2")
      qqline(W2b, col = 2,lwd=2,lty=2)

      hist(W2blam, main = "Histogram of Transformed Y2", xlab = "Transformed Y2")
      qqnorm(W2blam, main = "Normal Q-Q Plot of Transformed Y2")
      qqline(W2blam, col = 2,lwd=2,lty=2)

      par(mar = c(5.1, 4.1, 4.1, 2.1))

      # ROC PLOT 1

      par(mfrow = c(1, 2))

      plot(linspace(0,1,1000),
           roc1(linspace(0,1,1000)),
           main="Marker1",
           xlab="FPR = 1 - Specificity",
           ylab="TPR = Sensitivity",
           type="l",
           col="red")

      lines(linspace(0,1,1000),
            linspace(0,1,1000),
            type="l",
            lty=2)

      lines(roc.curve(marker1, D), col="red", lty = 3, lwd = 1.5)

      # LEGEND

      legend("bottomright", legend=c(paste("Box-Cox ROC Estimate"),
                                     paste("Empirical ROC Estimate")),
             col=c("red", "red"), lty=c(1, 3), lwd = c(1, 1.5), pch = c(NA, NA), cex=0.8)

      # ROC PLOT 2

      plot(linspace(0,1,1000),
           roc2(linspace(0,1,1000)),
           main="Marker2",
           xlab="FPR = 1 - Specificity",
           ylab="TPR = Sensitivity",
           type="l",
           col="black")

      lines(linspace(0,1,1000),
            linspace(0,1,1000),
            type="l",
            lty=2)

      lines(roc.curve(marker2, D), col="black", lty = 3, lwd = 1.5)

      # LEGEND

      legend("bottomright", legend=c(paste("Box-Cox ROC estimate"),
                                     paste("Empirical ROC estimate")),
             col=c("black", "black"), lty=c(1, 3), lwd = c(1, 1.5), pch = c(NA, NA), cex=0.8)

      par(mfrow = c(1, 1))

    }

    return(list(res_shapiro = res_shapiro,
                transformation.parameter.1 = lam[1],
                transx1 = W1alam,
                transy1 = W1blam,
                transformation.parameter.2 = lam[2],
                transx2 = W2alam,
                transy2 = W2blam,
                pval_x1 = pvalues[1],
                pval_y1 = pvalues[2],
                pval_transx1 = pvalues[3],
                pval_transy1 = pvalues[4],
                pval_x2 = pvalues[5],
                pval_y2 = pvalues[6],
                pval_transx2 = pvalues[7],
                pval_transy2 = pvalues[8],
                roc1 = roc1,
                roc2 = roc2))

  }
}
