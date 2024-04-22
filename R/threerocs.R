threerocs=function(marker, D, plots)
{

  x=marker[D==0];
  y=marker[D==1];


  # Plot the empirical ROC curve
  roc_obj <- roc(D, marker, quiet = TRUE)
  auc_value <- auc(roc_obj)
  tpr=roc_obj$sensitivities
  fpr=1-roc_obj$specificities

  if (plots == "on") {
    plot(fpr, tpr, type = "l", col = "blue",
         xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR = 1 - Specificity",
         ylab = "TPR = Sensitivity", main = "Three Estimates of the ROC Curve")
  }

  t=linspace(0,1,10000) #set a grid for the FPR

  ############## BOX-COX #############################
  roctbc=checkboxcox(marker,D, plots="off", printShapiro = FALSE);
  myroc=rocboxcox(marker,D,0.05,plots="off",FALSE)
  rocboxc=myroc$rocbc
  AUCbc=myroc$AUC
  if (plots == "on") {lines(t,roctbc$roc(t),col="red")}

  ############## Metz    #############################
  curve2= roc_curves(D,marker,method="binormal")
  rocmrm=points(curve2, metric="specificity", values=1-c(t))
  rocmrm$TPR
  AUCmrm=trapz(t,rocmrm$TPR)
  if (plots == "on") {lines(t,rocmrm$TPR,col="forestgreen", lty = "dashed")}

  #legend(x = "right", legend=c("one","two","three"))

  if (plots == "on") {
    legend("bottomright", legend = c(paste("Empirical (AUC = ", formattable(auc_value, digits = 4, format = "f"), ")", sep = ""),
                                     paste("Box-Cox (AUC = ", formattable(AUCbc, digits = 4, format = "f"), ")", sep = ""),
                                     paste("Metz (AUC = ", formattable(AUCmrm, digits = 4, format = "f"), ")", sep = "")),
           col = c("blue", "red", "forestgreen"),
           lty = c("solid", "dotted", "dashed"),
           cex = 0.8)
  }

  return(list(AUC_Empirical = as.numeric(auc_value),
              AUC_Metz = AUCmrm,
              AUC_BoxCox = AUCbc))

}
