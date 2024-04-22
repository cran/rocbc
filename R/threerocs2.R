threerocs2 <- function (marker1, marker2, D, plots) {


  erf <- function (x) 2 * pnorm(x * sqrt(2)) - 1
  erfc <- function (x) 2 * pnorm(x * sqrt(2), lower.tail = FALSE)
  erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
  erfcinv <- function (x) qnorm(x/2, lower.tail = FALSE)/sqrt(2)




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

  lam1=lam[1]
  lam2=lam[2]

  cova= cov(Walam)*(na-1)/na;cova=cova[1,2];
  covb= cov(Wblam)*(nb-1)/nb;covb=covb[1,2];


  s1alam=sqrt( var(W1alam)*(na-1)/na )
  s1blam=sqrt( var(W1blam)*(nb-1)/nb )
  s2alam=sqrt( var(W2alam)*(na-1)/na )
  s2blam=sqrt( var(W2blam)*(nb-1)/nb )

  m1alam=mean(W1alam)
  m1blam=mean(W1blam)
  m2alam=mean(W2alam)
  m2blam=mean(W2blam)

  roc1bc<-function(t){
    1-pnorm(qnorm(1-t,mean=mean(W1alam),sd=s1alam),mean=mean(W1blam),sd=s1blam)
  }

  roc2bc<-function(t){
    1-pnorm(qnorm(1-t,mean=mean(W2alam),sd=s2alam),mean=mean(W2blam),sd=s2blam)
  }


  #######################  MRM ###################


  rating=c(marker1,marker2)
  truth=c(D,D)
  reader=c(linspace(1,1,length(marker1)),linspace(2,2,length(marker2)))
  curves_binorm <- roc_curves(truth, rating,
                              groups = list(Reader = reader),
                              method = "binormal")

  params_binorm <- MRMCaov::parameters(curves_binorm)
  t=linspace(0,1,1000);

  roct_mrm_1=pnorm(c(params_binorm$a[1])+c(params_binorm$b[1])*qnorm(t))
  roct_mrm_2=pnorm(c(params_binorm$a[2])+c(params_binorm$b[2])*qnorm(t))


  # Plot the empirical ROC curve
  roc_obj1 <- roc(D, marker1, quiet = TRUE)
  auc_emp1 <- auc(roc_obj1)
  tpr1=roc_obj1$sensitivities
  fpr1=1-roc_obj1$specificities

  roc_obj2 <- roc(D, marker2, quiet = TRUE)
  auc_emp2 <- auc(roc_obj2)
  tpr2=roc_obj2$sensitivities
  fpr2=1-roc_obj2$specificities

  AUC_empirical_marker1=auc_emp1
  AUC_empirical_marker2=auc_emp2

  if (plots == "on") {

    plot(t, roc1bc(t), type="l", lty = "dotted", col="indianred2", lwd=1.5, ylab="TPR = Sensitivity", xlab="FPR = 1 - Sensitivity")
    lines(t, roc2bc(t), type="l", lty = "dotted", col="red4", lwd=1.5)
    lines(t, roct_mrm_1, type="l", lty = "dashed", lwd=1.5, col="green2")
    lines(t, roct_mrm_2, type="l", lty = "dashed", lwd=1.5, col="darkgreen")
    lines(fpr1, tpr1, type="l", lty = "solid", lwd=1, col="deepskyblue")
    lines(fpr2, tpr2, type="l",lty = "solid", lwd=1, col="blue4")

  }


  #plot(fpr1, tpr1, type = "l", col = "magenta", lwd = 2,
  #     xlim = c(0, 1), ylim = c(0, 1), xlab = "False Positive Rate (FPR)",
  #     ylab = "True Positive Rate (TPR)", main = "Three Estimates of the ROC Curve")

  AUC1bc=trapz(t,roc1bc(t))
  AUC2bc=trapz(t,roc2bc(t))
  AUCmetz1=trapz(t,roct_mrm_1)
  AUCmetz2=trapz(t,roct_mrm_2)
  AUCemp1=abs(trapz(fpr1,tpr1))
  AUCemp2=abs(trapz(fpr2,tpr2))




  #legend("bottomright", legend = c("Box-Cox marker 1", "Box-Cox marker 2", "Metz marker 1",, "Metz marker 2","Empirical marker 1","Empirical marker 2"),
  #         col = c("black","red","blue","green", "magenta", "cyan"), lty = c("solid", "solid", "dashed","dashed","dashed","dashed"),
  #         lwd = 2, bty = "n")


  #CORRECT ONE:

  if (plots == "on") {

    legend("bottomright",
           legend = c(paste("Box-Cox Marker 1 (AUC = ", formattable(AUC1bc, digits = 4, format = "f"), ")", sep = ""),
                      paste("Box-Cox Marker 2 (AUC = ", formattable(AUC2bc, digits = 4, format = "f"), ")", sep = ""),
                      paste("Metz Marker 1 (AUC = ", formattable(AUCmetz1, digits = 4, format = "f"), ")", sep = ""),
                      paste("Metz Marker 2 (AUC = ", formattable(AUCmetz2, digits = 4, format = "f"), ")", sep = ""),
                      paste("Empirical Marker 1 (AUC = ", formattable(AUCemp1, digits = 4, format = "f"), ")", sep = ""),
                      paste("Empirical Marker 2 (AUC = ", formattable(AUCemp2, digits = 4, format = "f"), ")", sep = "")),
           col = c("indianred2", "red4", "green2", "darkgreen", "deepskyblue", "blue4"),
           lty = c("dotted", "dotted", "dashed", "dashed", "solid", "solid"),
           cex = 0.8, lwd = 2)


    title(main = "Three Estimates of the ROC Curve (Two Markers)")

  }

  return(list(AUC_BoxCox1 = AUC1bc,
              AUC_BoxCox2 = AUC2bc,
              AUC_Metz1 = AUCmetz1,
              AUC_Metz2 = AUCmetz2,
              AUC_Empirical1 = AUCemp1,
              AUC_Empirical2 = AUCemp2))


}
