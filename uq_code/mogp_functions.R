# Saving a list of MOGP emulators
SaveMulti <- function(ems, filename){
  ell <- length(ems)
  dir.create(filename)
  for (i in 1:ell){
    tmp_em <- ems[[i]]
    save_ExUQmogp(tmp_em, filename = paste0(filename, '/', filename, '_', i))
  }
}

# Loading a list of MOGP emulators
# directory is location, filename is common tag
LoadMulti <- function(directory, filename){
  ell <- length(list.files(paste0(directory, '/'))) / 2
  tmp_list <- NULL
  for (i in 1:ell){
    tmp_list[[i]] <- load_ExUQmogp(filename = paste0(directory, '/', filename, '_', i))
  }
  return(tmp_list)
}

# Create basis (DataBasis object) and other preliminary calculations (inverting W)
CreateBasis <- function(data, type = 'SVD', W = NULL, RemoveMean = TRUE){
  # Invert W, creating Winv and tagging with attributes to enable fast calculations later on
  if (is.null(W)){
    W <- diag(nrow(data))
  }
  
  Winv <- GetInverse(W)
  
  if (type %in% c('SVD', 'svd', 'L2')){
    DataBasis <- MakeDataBasis(data, weightinv = NULL, RemoveMean = RemoveMean, StoreEigen = TRUE)
  }
  
  if (type %in% c('WSVD', 'wsvd', 'wSVD')){
    DataBasis <- MakeDataBasis(data, weightinv = Winv, RemoveMean = RemoveMean, StoreEigen = TRUE)
  }
  
  DataBasis$Type <- type
  
  if (DataBasis$Type %in% c('svd', 'L2')){
    DataBasis$Type <- 'SVD'
  }
  
  if (DataBasis$Type %in% c('WSVD', 'wsvd')){
    DataBasis$Type <- 'wSVD'
  }
  
  return(DataBasis)
}





AssessBasis <- function(DataBasis, Obs){
  
  max_q <- dim(DataBasis$tBasis)[2]
  PlotData <- data.frame(Size = 1:max_q, Error = numeric(max_q), Explained = numeric(max_q))
  
  if (!is.null(DataBasis$scaling)){
    PlotData$Error <- errors(DataBasis$tBasis, Obs - DataBasis$EnsembleMean, DataBasis$Winv)*DataBasis$scaling^2
  }
  else {
    PlotData$Error <- errors(DataBasis$tBasis, Obs - DataBasis$EnsembleMean, DataBasis$Winv)
  }
  
  if (is.null(DataBasis$Winv)){
    var_sum <- crossprod(c(DataBasis$CentredField))
  }
  else {
    var_sum <- sum(diag(t(DataBasis$CentredField) %*% DataBasis$Winv %*% DataBasis$CentredField))
  }
  
  for (i in 1:max_q){
    PlotData$Explained[i] <- VarExplained(DataBasis$tBasis[,1:i], DataBasis$CentredField, DataBasis$Winv, total_sum = var_sum)
  }
  
  PlotData$Explained <- round(PlotData$Explained, 10) # to make sure plots when = 1
  
  chi_bound <- qchisq(0.995, nrow(DataBasis$tBasis)) / nrow(DataBasis$tBasis)
  max_y <- max(c(PlotData$Error, chi_bound + 0.05))
  
  var_plot <- ggplot(data = PlotData, aes(x = Size)) +
    geom_line(aes(y = Error), col = 'red') +
    geom_line(aes(y = Explained * max_y), col = 'blue') +
    xlab('Basis size') +
    scale_y_continuous(
      #name = 'Error',
      limits = c(0,max_y),
      #breaks = seq(from = 0, to = 1*max_y, by = 0.25*max_y),
      #labels = NULL,
      sec.axis = sec_axis(trans=~./max_y, name="Explained")) +
    #geom_vline(xintercept = q) +
    geom_hline(yintercept = chi_bound, linetype = 'dashed', col = 'red') +
    geom_hline(yintercept = 0.9*max_y, linetype = 'dashed', col = 'blue')
  #theme(panel.grid.major = element_blank(), 
  #      panel.grid.minor = element_blank())
  
  # Then don't need to do ExplainT, as already have the error/explained combinations
  # Give list of (threshold, q)
  thresholds <- c(0.8, 0.85, 0.9, 0.95, 0.99, 0.999)
  q <- numeric(length(thresholds))
  for (j in 1:length(q)){
    q[j] <- which(PlotData$Explained >= thresholds[j])[1]
  }
  
  return(list(plot = var_plot,
              Errors = PlotData,
              Truncations = data.frame(Threshold = thresholds,
                                       q = q)))
}


# Use this to create tData, projected obs, etc.
GetEmulatableData <- function(Design, DataBasis, BasisSize = NULL, Obs, Noise = TRUE){
  if(Noise){
    Noise <- runif(length(Design[,1]),-1,1)
    Design <- cbind(Design, Noise)
  }
  tData <- Project(DataBasis$CentredField, DataBasis$tBasis[,1:BasisSize], weightinv = DataBasis$Winv)
  colnames(tData) <- paste0('C', 1:ncol(tData)) 
  tData <- cbind(Design, tData)
  return(tData)
}


MVImplausibilityMOGP <- function(NewData, Emulator, DataBasis, Discrepancy, Obs, ObsErr){
  tEmulator <- Emulator$mogp$predict(as.matrix(NewData), deriv=FALSE)
  
  HistoryMatch(DataBasis, Obs, t(tEmulator$mean), t(tEmulator$unc), ObsErr, Discrepancy, weightinv = DataBasis$Winv)
  #n <- ncol(tEmulator$mean)
  #tpreds <- t(tEmulator$mean)
  #tunc <- t(tEmulator$unc)
  #as.numeric(mclapply(1:n, function(k) ImplCoeff(tpreds[k,], tunc[k,], Obs, ObsErr, Discrepancy)))
}


BuildNewEmulatorsFourier <- function(tData, HowManyEmulators, 
                                     additionalVariables=NULL, 
                                     Choices = lapply(1:HowManyEmulators,
                                                      function(k) choices.default), meanFun = "linear", 
                                     kernel = c("Gaussian"),
                                     SpecFormula = NULL,
                                     Active = NULL,
                                     ...){
  #'@description Builds MO_GP emulators for a data frame 
  #'@param tData This is the data frame you wish to use. The format should be D+1+Q where Q is the number of targets each occupying one of the last Q columns, D is the number of inputs occupying the first D columns. The D+1th column should be a vector of random normal numbers (between -1 and 1) called "Noise". We use it in our LM code to stop our algorithms overfitting by selecting signals by chance. All variables should have names.
  #'@param HowManyEmulators How many emulators are required. The code will fit the first HowManyEmulators GPs up to Q.
  #'@param additionalVariables Are there variables that must be "active" (e.g. you really want to use them in a decision problem or similar) or should be included? Often all variables are specified here, but it defaults to NULL.
  #'@param Choices A list of choices with each of the HowManyEmulators elements being a choice list compatible with GetPriors() (see the documentation for GetPriors)
  #'@param meanFun Currently a single string either "fitted", or "linear" ("constant" to come). If "fitted", our custom global mean functions are fitted and then used to fit the GP. Recommended for higher dimensions and for history matching. If "linear", a linear mean function is fitted to all emulators and using additionalVariables. A list implementation will be considered in future versions. Could also make a list where it could be a formula (is.formula)
  #'@param kernel A vector of strings that corresponds to the type of kernel either "Gaussian", or "Matern52"
  #'Default is to use Gaussian kernel for all emulators.
  #'@details If mean functions are not given (an option that will be added soon) our automatic LM code fits a global mean function for each metric. Get Priors is then used to extract the priors before we establish an MOGP and fit the parameters by MAP estimation. MAP improves on MLE here as we avoid the ridge on the likelihood surface for GPs.
  #'@return A list with 2 elements. 1 the mogp, 2 a list containing the elements used for fitting: the mean functions (containing element linModel as the lm object), the design, a list of active input indices for each emulator (to be used for subsetting the design to produce diagnostics), and the prior choices.
  ###Mean function selection###
  lastCand <- which(names(tData)=="Noise")
  if(length(lastCand)<1)
    stop("tData should have a column called 'Noise' separating the inputs and outputs.")
  if(is.null(HowManyEmulators))
    HowManyEmulators <- length(names(tData)) - lastCand
  if(!(HowManyEmulators == length(names(tData)) - lastCand)){
    tData <- tData[,c(1:lastCand,(lastCand+1):(lastCand+HowManyEmulators))]
  }
  if(meanFun =="fitted"){
    tdfs <- DegreesFreedomDefault(Choices, N=length(tData[,1]))
    lm.list = lapply(1:HowManyEmulators, function(k) 
      try(EMULATE.lm(Response=names(tData)[lastCand+k],
                     tData=tData, tcands=names(tData)[1:lastCand],
                     tcanfacs=NULL, 
                     TryFouriers=Choices[[k]]$lm.tryFouriers, 
                     maxOrder=Choices[[k]]$lm.maxOrder,
                     maxdf = tdfs[k])))
  }
  else if(meanFun == "linear"){
    if(is.null(additionalVariables))
      stop("When specifying linear meanFun, please pass the active inputs into additionalVariables")
    linPredictor <- paste(additionalVariables,collapse="+")
    lm.list = lapply(1:HowManyEmulators, function(k) list(linModel=eval(parse(text=paste("lm(", paste(names(tData[lastCand+k]), linPredictor, sep="~"), ", data=tData)", sep="")))))
  }
  else if(meanFun == 'spec'){
    lm.list = lapply(1:HowManyEmulators, function(k) list(linModel = eval(parse(text=paste("lm(", paste(names(tData[lastCand+k]), SpecFormula[k], sep="~"), ", data=tData)", sep="")))))
  }
  
  
  else{
    stop("meanFun must either be 'fitted' or 'linear' in this version")
  }
  ###Prepare the data for MOGP### 
  tfirst <- lastCand + 1
  target_names <- names(tData)[tfirst:length(names(tData))]
  target_list <- extract_targets(tData[,-which(names(tData)=="Noise")], target_names)
  inputs <- target_list[[1]]
  targets <- target_list[[2]]
  inputdict <- target_list[[3]]
  d <- dim(inputs)[2]
  
  if(meanFun=="fitted"){
    ActiveVariableIndices <- lapply(lm.list, function(tlm) which((names(tData)%in%additionalVariables) | (names(tData)%in%tlm$Names) | (names(tData) %in% names(tlm$Fouriers))))
  }
  else if(meanFun == "linear"){
    ActiveVariableIndices <- lapply(lm.list, function(tlm) which(names(tData)%in%additionalVariables))
  }
  else if(meanFun == "spec"){
    ActiveVariableIndices <- lapply(Active, function(tlm) which(names(tData)%in%tlm))
  }
  ###Prepare the mean functions for MOGP### 
  #mean_func.list.MGP <- lapply(lm.list, function(e) FormulaToString(formula(e$linModel)))
  
  mean_func.list.MGP <- lapply(lm.list, function(e) FormulaToStringConstant(formula(e$linModel)))
  
  ###Establish the priors for the emulators###
  Priors <- lapply(1:HowManyEmulators, function(l) GetPriorsConstant(lm.list[[l]], d=d, Choices[[l]], ActiveVariableIndices[[l]]))
  
  ###Establish the kernel types for MOGP###
  if(length(kernel) == 1) {
    Kernels <- lapply(1:HowManyEmulators, function(l) GetKernel(kernel))
  } else {
    Kernels <- lapply(1:HowManyEmulators, function(l) GetKernel(kernel[l])) 
  }
  
  if (HowManyEmulators == 1){
    targets <- targets - predict(lm.list[[1]]$linModel, tData)
  }
  else {
    for (cc in 1:HowManyEmulators){
      targets[cc,] <- targets[cc,] - predict(lm.list[[cc]]$linModel, tData)
    }
  }
  
  ###Establish and fit the MOGP###
  Emulators <- mogp_emulator$MultiOutputGP(inputs, targets, mean = mean_func.list.MGP,
                                           priors = Priors, inputdict = inputdict,
                                           nugget = lapply(Choices,function(e) e$Nugget), 
                                           kernel = Kernels)
  Emulators <- mogp_emulator$fit_GP_MAP(Emulators)
  
  ###Prepare return objects###
  Design <- tData[,1:(lastCand-1), drop = FALSE]
  fitting <- list(lm.object = lm.list,
                  Design = Design, 
                  ActiveIndices = ActiveVariableIndices,
                  PriorChoices = Choices)
  
  return(list(mogp = Emulators, # call mogp
              fitting.elements= fitting))
}


FormulaToStringConstant <- function(tformula){
  #'@description Parse and R formula into a string
  #'@param tformula an R formula of the form y~I(x1)+I(x_2)+ etc
  #'@return A string version to be passed into MOGP
  f = as.character(tformula)
  paste(f[2], f[1], 1, sep="")
}

GetPriorsConstant <- function(lm.emulator, d, Choices, ActiveVariables){
  #'@description This function constructs a subjective prior for the parameters of a GP emulator to be fit by MO_GP. 
  #'@param lm.emulator A lm emulator list (see AutoLMCode). The only required element of this list is the linModel component. lm.emulator$linModel is an lm object fitted to the target data. Custom lm objects can be specified using lm.emulator = list(linModel=lm(...), entering your own formulae). It is suggested that this is done elsewhere and passed here.
  #'@param d the number of input parameters
  #'@param Choices A list containing hyperprior choices that control the subjective prior.
  #'@param NonInformativeRegression If TRUE, a uniform prior is used for the regression parameters.
  #'@param NonInformativeCorrelationLengths If TRUE, a uniform prior is used for the correlation parameters. 
  #'@param NonInformativeSigma If TRUE, a uniform prior is used for the sigma parameter.
  #'@param NonInformativeNugget If TRUE, a uniform prior is used for the regression parameters.
  #'@param BetaRegressMean Prior mean for the regression coefficients. The intercept is given a uniform prior, all other regression terms get a Normal prior with mean BetaRegressMean and variance "BetaRegressSigma"
  #'@param BetaRegressSigma Prior variance for the regression coefficients.
  #'@param DeltaActiveMean The mean of a lognormal prior for the active inputs (see details). 
  #'@param DeltaActiveSigma The variance of a lognormal prior for the active inputs (see details).
  #'@param DeltaInactiveMean The mean of a lognormal prior for the inactive inputs (see details).
  #'@param DeltaInactiveSigma The variance of a lognormal prior for the inactive inputs (see details).
  #'@param NuggetProportion What proportion of the data variability is suspected to be nugget. Only used if Nugget="fit"
  #'@param Nugget either a number or "fit", "fixed" or "adaptive". This is seen by MOGP which will only fit the nugget (and hence require a prior) if "fit" is chosen. If the others, the nugget is fixed to the value given or fitted to avoid numerical instabilities adaptively. (See MOGP documentation.)
  #'@param ActiveVariables Indices indicating which parameters in the data are active
  #'@return A list of priors in the following order. A prior for the intercept (defaults to NULL), p priors for the regression coefficients, d priors for the correlation lengths, a prior for sigma squared and a prior for the nugget (NULL if Choices$Nugget != "fit")
  #'@details The linear mean has p coefficients not associated with the intercept. Each coefficient is Beta ~ Normal(BetaRegressMean, BetaRegressSigma). Our default is based on mean 0 variance 10 and is weakly informative.
  #'@details Input parameters are classified either as active or inactive. These definitions depend on whether the parameters were highlighted by preliminary fitting algorithms in our linear modelling code (see AutoLMcode) and whether the user asks for them specifically to be included (so note they may do nothing, but still get an "active" prior. There are 2 classes of prior we use for the correlation lengths. The first is for the "active" ones: log(delta) ~ Normal(DeltaActiveMean, DeltaActiveSigma) (so that delta is lognormal). default this to N(0,0.5) on the basis that we expect there to be correlation and so that they are not too long (we penalise the rigde on the likelihood surface). Inactive parameters get the same type of prior with very strong default values N(5,0.005), giving the whole distribtion at contributions to the correlation near 1. As our covariance functions are separable, the inactive variables therefore don't alter the fit. All values for these parameters are controllable by the user.
  #'@details. Sigma^2 ~ InvGamma(M,V) where we use a parameterisation of the inverse gamma based on the mode and variance. The mode is chosen so reflect the variance we can explain with simple linear fits and the variance is set using invgamMode with a bound based on the total variability in the original data. The idea behind the prior is that we allow Sigma^2 to be as large as the variance of the data so that our emulator explains none of the variability, but we expect to explain as much as could be explained by a preliminary fit of a basic model.
  #'@details The nugget distribution, if required, is found as with sigma but using the choice NuggetProportion to reflect what percentage of variability we might expect to be nugget. We may add the ability to add a user prior here. Only really important for internal variability models like climate models. Most deterministic models should use "adaptive" or "fixed" nugget. 
  p <- 0
  Priors <- lapply(1:(1+p+d+1+1), function(e) NULL)
  if(!is.null(Choices$intercept))
    print("NULL intercept fitted by default: change code")
  
  Priors[[1]] <- mogp_priors$NormalPrior(0., 0.0001)
  #Regression priors
  #if(!(Choices$NonInformativeRegression)){
  #  Betas <- lapply(1:p, function(e) mogp_priors$NormalPrior(Choices$BetaRegressMean, Choices$BetaRegressSigma))
  #  for(i in 2:(p+1)){#the first element is NULL for the intercept term
  #    Priors[[i]] <- Betas[[i-1]]
  #  }
  #}
  #Correlation length priors
  if(!(Choices$NonInformativeCorrelationLengths)){
    Deltas <- lapply(1:d, function(k) {if(k %in% ActiveVariables) 
    {mogp_priors$NormalPrior(Choices$DeltaActiveMean,Choices$DeltaActiveSigma)} 
      else {mogp_priors$NormalPrior(Choices$DeltaInactiveMean,Choices$DeltaInactiveSigma)}})
    for(j in (p+2):(p+1+d)){
      Priors[[j]] <- Deltas[[j-p-1]]
    }
  }
  #Sigma and nugget priors
  ModeSig <- var(lm.emulator$linModel$residuals)
  boundSig <- ModeSig/(1-summary(lm.emulator$linModel)$r.squared)
  if(!(Choices$NonInformativeSigma)){
    SigmaParams <- invgamMode((1-Choices$NuggetProportion)*boundSig,ModeSig)
    Sigma <- mogp_priors$InvGammaPrior(SigmaParams$alpha,SigmaParams$beta)
    Priors[[d+p+2]] <- Sigma
  }
  if(Choices$Nugget=="fit"){
    if(!Choices$NonInformativeNugget){
      NuggetParams <- invgamMode(Choices$NuggetProportion*boundSig ,Choices$NuggetProportion*ModeSig)
      Nugget <- mogp_priors$InvGammaPrior(NuggetParams$alpha, NuggetParams$beta)
      Priors[[d+p+3]] <- Nugget
    }
  }
  return(Priors)
}


PredictMOGP <- function(emulator, design){
  n_em <- emulator$mogp$n_emulators
  gp_pred <- emulator$mogp$predict(as.matrix(design), deriv=FALSE)
  for (i in 1:n_em){
    gp_pred$mean[i,] <- gp_pred$mean[i,] + predict(emulator$fitting.elements$lm.object[[i]]$linModel, design)
  }
  return(gp_pred)
}

ValidateMOGP <- function(emulator, ValidationData, which.emulator=1, IndivPars = FALSE, SDs = 2){
  require(sfsmisc)
  
  preds <- PredictMOGP(emulator, ValidationData[,1:(which(colnames(ValidationData) == 'Noise') - 1)])
  response <- ValidationData[,which(colnames(ValidationData) == 'Noise') + which.emulator]
  
  preds$upper95 <- preds$mean + SDs*sqrt(preds$unc)
  preds$lower95 <- preds$mean - SDs*sqrt(preds$unc)
  upp <- max(c(preds$upper95[which.emulator,], response))
  low <- min(c(preds$lower95[which.emulator,], response))
  
  if (IndivPars == TRUE){
    p <- length(emulator$fitting.elements$ActiveIndices[[which.emulator]]) + 1
    if(p<2){
      par(mfrow = c(1, 1), mar=c(4, 4, 1, 1))
    }
    else if(p<3){
      par(mfrow = c(1, 2), mar=c(4, 4, 1, 1))
    }
    else if(p<4){
      par(mfrow = c(1, 3), mar=c(4, 4, 1, 1))
    }
    else if(p <5){
      par(mfrow = c(2, 2), mar=c(4, 4, 1, 1))
    }
    else if(p<7){
      par(mfrow = c(2, 3), mar=c(4, 4, 1, 1))
    }
  }
  errbar(preds$mean[which.emulator,], preds$mean[which.emulator,], preds$upper95[which.emulator,], preds$lower95[which.emulator,], cap = 0.015, pch=20, 
         ylim=c(low,upp), main="",xlab = "Prediction",ylab="Data")
  points(preds$mean[which.emulator,], response, pch=19,
         col = ifelse(response > preds$upper95[which.emulator,] | response < preds$lower95[which.emulator,], "red", "green"))
  if (IndivPars == TRUE){
    active <- emulator$fitting.elements$ActiveIndices[[which.emulator]]
    
    for (i in 1:length(active)){
      errbar(ValidationData[,active[i]], preds$mean[which.emulator,], preds$upper95[which.emulator,], preds$lower95[which.emulator,], cap = 0.015, pch=20, 
             ylim=c(low,upp), xlab = "Input",ylab="Prediction", main = paste(colnames(ValidationData)[active[i]]))
      points(ValidationData[,active[i]], response, pch=19,
             col = ifelse(response > preds$upper95[which.emulator,] | response < preds$lower95[which.emulator,], "red", "green"))
    }
  }
}

#' Matrix inversion via cholesky decomposition
#'
#' Inverts matrix W, assigning attributes for whether W is diagonal, to speed up other calculations.
#'
#' @param W square positive definite variance matrix
#'
#' @return Inverse of W, with attributes 'identity' and 'diagonal', used by other functions in the package to make calculations more efficient.
#'
#' @examples Winv <- GetInverse(diag(100))
#' attributes(Winv) # diagonal = TRUE, identity = TRUE
#'
#' Winv2 <- GetInverse(runif(100,0.1,1)*diag(100))
#' attributes(Winv2) # diagonal = TRUE, identity = FALSE
#'
#' Winv3 <- GetInverse(seq(0.1,1,length=100) %*% t(seq(0.1,1,length=100)) + 0.1*diag(100))
#' attributes(Winv3) # diagonal = FALSE, identity = FALSE
#'
#' @export
GetInverse <- function(W){
  diagmat <- all(W[lower.tri(W)] == 0, W[upper.tri(W)] == 0)
  if (diagmat == TRUE){
    InvW <- diag(1 / diag(W))
  }
  else {
    Q <- chol(W)
    y <- backsolve(Q, diag(dim(W)[1]), transpose = TRUE)
    InvW <- crossprod(y, y)
  }
  attr(InvW, 'diagonal') <- diagmat
  if (all(diag(InvW) == 1) & diagmat == TRUE){
    attr(InvW, 'identity') <- TRUE
  }
  else {
    attr(InvW, 'identity') <- FALSE
  }
  return(InvW)
}

#' Formatting data
#'
#' Formats data so that it is in the correct form for use in other functions, and calculates the (weighted) SVD basis of the ensemble
#'
#' @param data a matrix containing individual fields in the columns (i.e. the matrix has dimension lxn)
#' @param weightinv the inverse of lxl positive definite weight matrix W. If NULL, the identity matrix is used
#' @param RemoveMean if TRUE, centres the data prior to calculating the basis
#' @param StoreEigen if TRUE, stores Q, lambda from eigendecomposition of W (in order to make later calculations more efficient)
#'
#' @return \item{tBasis}{The (weighted) SVD basis of the centred ensemble if RemoveMean = TRUE, of the original data otherwise}
#' \item{CentredField}{The centred data if RemoveMean = TRUE, the original data otherwise.}
#' \item{EnsembleMean}{The mean across the columns of the data. A zero vector if RemoveMean = FALSE}
#' \item{}
#'
#' @export
MakeDataBasis <- function(data, weightinv = NULL, W = NULL, RemoveMean = TRUE, StoreEigen = TRUE){
  if (RemoveMean == TRUE){
    EnsembleMean <- apply(data, 1, mean)
    CentredField <- 0*data
    for (i in 1:dim(data)[2]){
      CentredField[,i] <- data[,i] - EnsembleMean
    }
  }
  else {
    EnsembleMean <- c(rep(0, dim(data)[1]))
    CentredField <- data
  }
  #if (is.null(weightinv)){
  #  weightinv <- diag(dim(data)[1])
  #}
  if (is.null(W)){
    tSVD <- wsvd(t(CentredField), weightinv = weightinv)
    tBasis <- tSVD$v
    if (StoreEigen == TRUE){
      Q <- tSVD$Q
      Lambda <- tSVD$Lambda
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Q = Q, Lambda = Lambda))
    }
    else {
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean))
    }
  }
  else if (!is.null(W) & is.null(weightinv)){
    eig <- eigen(W)
    Q <- eig$vectors
    Lambda <- 1 / eig$values
    Winv <- Q %*% diag(Lambda) %*% t(Q)
    attr(Winv, 'diagonal') <- FALSE
    attr(Winv, 'identity') <- FALSE
    tSVD <- wsvd(t(CentredField), weightinv = Winv, Q = Q, Lambda = Lambda)
    tBasis <- tSVD$v
    if (StoreEigen == TRUE){
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Q = Q, Lambda = Lambda, Winv = Winv))
    }
    else {
      return(list(tBasis = tBasis, CentredField = CentredField, EnsembleMean = EnsembleMean, Winv = Winv))
    }
  }
}


#' Weighted singular value decomposition
#'
#' Calculates the SVD basis across the output, given the inverse of W.
#'
#' @param data n x l matrix to calculate basis from (i.e. rows are output fields).
#' @param weightinv l x l inverse of W. If NULL, calculates standard SVD.
#' @param Q l x l matrix from eigen decomposition of W^{-1}, if already have this then speeds up calculation of basis
#' @param Lambda vector from eigen decomposition of W^{-1}, if already have this then speeds up calculation of basis
#'
#' @return The weighted SVD of the data.
#'
wsvd <- function(data, weightinv = NULL, Q = NULL, Lambda = NULL){
  if (is.null(weightinv)){
    svd_output <- svd(data)
  }
  else {
    stopifnot(dim(data)[2] == dim(weightinv)[1])
    if (is.null(Q) & attributes(weightinv)$diagonal == FALSE){
      eig <- eigen(weightinv)
      Q <- eig$vectors
      Lambda <- eig$values
      data_w <- data %*% Q %*% diag(sqrt(Lambda)) %*% t(Q)
      svd_output <- svd(data_w)
      svd_output$v <- t(t(svd_output$v) %*% Q %*% diag(1 / sqrt(Lambda)) %*% t(Q))
      svd_output$Q <- Q
      svd_output$Lambda <- Lambda
    }
    else if (is.null(Q) & attributes(weightinv)$diagonal == TRUE){
      diag_values <- diag(weightinv)
      data_w <- data %*% diag(sqrt(diag_values))
      svd_output <- svd(data_w)
      svd_output$v <- t(t(svd_output$v) %*% diag(1 / sqrt(diag_values)))
    }
    else if (!is.null(Q)){
      data_w <- data %*% Q %*% diag(sqrt(Lambda)) %*% t(Q)
      svd_output <- svd(data_w)
      svd_output$v <- t(t(svd_output$v) %*% Q %*% diag(1 / sqrt(Lambda)) %*% t(Q))
      svd_output$Q <- Q
      svd_output$Lambda <- Lambda
    }
  }
  return(svd_output)
}

errors <- function(basis, obs, weightinv=NULL){
  p <- dim(basis)[2]
  err <- numeric(p)
  if (is.null(weightinv)){
    weightinv <- diag(dim(basis)[1])
    attr(weightinv, 'diagonal') <- attr(weightinv, 'identity') <- TRUE
  }
  for (i in 1:p){
    err[i] <- ReconError(obs, basis[,1:i], weightinv)
  }
  return(err)
}

#' Calculating the proportion of data explained by a basis
#'
#' Calculates the proportion of the data that is explained by projection onto a basis.
#'
#' @param basis The basis
#' @param data The data to be explained
#' @param weightinv Inverse of W (identity if NULL)
#' @param total_sum The total sum of squares of the data with respect to W
#' @param psi t(original_basis) %*% weightinv %*% original_basis, where the new basis is a linear combination of some original basis
#' @param basis_lincom Vector of linear combinations (if new basis is a linear combination of some original basis)
#'
#' @return The proportion of variability in the data that is explained by the basis
#'
#' @export
VarExplained <- function(basis, data, weightinv = NULL, total_sum = NULL, psi = NULL, basis_lincom = NULL){
  coeffs <- t(CalcScores(data, basis, weightinv))
  recon <- basis %*% coeffs
  if (is.null(weightinv)){
    explained <- crossprod(c(recon))/crossprod(c(data))
  }
  else {
    if (is.null(psi)){
      if (attributes(weightinv)$diagonal == TRUE){
        explained_num <- sum(t(recon)^2 %*% diag(weightinv))
      }
      else {
        explained_num <- sum(diag(t(recon) %*% weightinv %*% recon))
      }
    }
    else {
      stopifnot(!is.null(basis_lincom))
      explained_num <- t(coeffs) %*% t(basis_lincom) %*% psi %*%
        basis_lincom %*% coeffs
      explained_num <- sum(diag(explained_num))
    }
    #explained_num <- 0
    #for (i in 1:dim(data)[2]){
    #  explained_num <- explained_num + t(recon[,i]) %*% weightinv %*% recon[,i]
    #}
    #explained_den <- 0
    #for (i in 1:dim(data)[2]){
    #  explained_den <- explained_den + t(data[,i]) %*% weightinv %*% data[,i]
    #}
    if (is.null(total_sum)){
      if (attributes(weightinv)$diagonal == TRUE){
        explained_den <- sum(t(data)^2 %*% diag(weightinv))
      }
      else {
        explained_den <- sum(diag(t(data) %*% weightinv %*% data))
      }
    }
    else {
      explained_den <- total_sum
    }
    explained <- explained_num / explained_den
  }
  return(explained)
}

#' Reconstruction error
#'
#' Calculates the reconstruction error, R_W(basis, obs), of the observations given a basis and W.
#'
#' @param obs The observations
#' @param basis Basis to project and reconstruct the observations with
#' @param weightinv Inverse of weight matrix W. If NULL (default), calculates the mean squared error
#' @param scale If TRUE, scales by the dimension (so analogous to mean squared error)
#'
#' @return The reconstruction error
#'
#' @export
ReconError <- function(obs, basis, weightinv = NULL, scale = TRUE){
  if (is.null(weightinv)){
    weightinv <- 0
    attr(weightinv, 'diagonal') <- attr(weightinv, 'identity') <- TRUE
  }
  field <- ReconObs(obs, basis, weightinv)
  A <- c(obs) - field
  mask <- which(is.na(A))
  if(length(mask)>0){
    A <- A[-mask]
  }
  if (scale == TRUE){
    s <- length(c(obs))-length(mask)
  }
  else {
    s <- 1
  }
  if (attributes(weightinv)$diagonal == FALSE){
    if(length(mask)>0){
      warning("Implicit assumption that weight specified on the full field even though applying a mask to missing obs/ensemble grid boxes")
      weightinv <- weightinv[-mask,-mask]
    }
    wmse <- (t(A) %*% weightinv %*% A)/ s
  }
  else {
    if (attributes(weightinv)$identity == TRUE){
      wmse <- crossprod(A)/ s
    }
    else {
      wmse <- crossprod(A/(1/diag(weightinv)), A)/ s
    }
  }
  return(as.numeric(wmse))
}

#' Project and reconstruct a given field
#'
#' Gives the reconstruction of a field using a basis, by projecting and back-projecting on this basis.
#'
#' @param field Vector over original field
#' @param basis Basis matrix
#'
#' @return Reconstruction of the original field.
#'
#' @examples
#'
#' @export
ReconField <- function(field, basis, ...){
  nb <- is.null(dim(basis))
  if(!nb)
    basis1 <- basis[,1]
  else
    basis1 <- basis
  field <- c(field)
  mask <- which(is.na(field-basis1))
  if(length(mask)>0){
    recons <- rep(NA, length(field))
    field <- field[-mask]
    if(nb)
      basis <- basis[-mask]
    else
      basis <- basis[-mask,]
    proj <- CalcScores(field, basis, ...)
    recons.partial <- Recon(proj, basis)
    recons[-mask] <- recons.partial
  }
  else{
    proj <- CalcScores(field, basis, ...)
    recons <- Recon(proj, basis)
  }
  return(recons)
}

ReconObs <- ReconField

#' Projection onto a basis
#'
#' Calculates the coefficients given by projecting data onto a basis
#'
#' @param data Data matrix to be projected, where each column is a representation on the original field
#' @param basis Basis matrix
#' @param weightinv If NULL, uses standard SVD projection. Otherwise, uses weighted projection.
#'
#' @return Matrix of basis coefficients
#'
#' @examples # First generate some data
#'
#' l <- 100 # dimension of output
#' n <- 10 # number of runs
#' DataBasis <- MakeDataBasis(data = matrix(runif(l*n), nrow=l, ncol=n), RemoveMean = TRUE) # data is 100x10
#'
#' # Project the (centred) ensemble onto the first 3 vectors of the SVD basis
#'
#' Coefficients <- CalcScores(data = DataBasis$CentredField, basis = DataBasis$tBasis[,1:3])
#'
#' # Instead of projecting using W = I, define a W with varying diagonal
#'
#' W <- runif(l, 1, 5) * diag(l) # 100x100 diagonal matrix
#' W_inv <- GetInverse(W) # inverse needed for projection
#' Coefficients_weighted <- CalcScores(data = DataBasis$CentredField, basis = DataBasis$tBasis[,1:3], weightinv = W_inv)
#'
#' @export
Project <- function(data, basis, weightinv = NULL){
  d <- dim(data)[2]
  if (is.null(d)){
    d <- 1
  }
  p <- dim(basis)[2]
  l <- dim(basis)[1]
  if (is.null(p)){
    p <- 1
  }
  if (d == 1){
    data <- as.vector(data)
  }
  if (is.null(weightinv)){
    weightinv <- 0 # just need to set as something that isn't NULL so can give attribute
    attr(weightinv, 'diagonal') <- attr(weightinv, 'identity') <- TRUE
  }
  if (attributes(weightinv)$identity == TRUE){
    V <- t(basis) %*% basis
    Q <- chol(V)
    y <- backsolve(Q, diag(p), transpose = TRUE)
    x <- backsolve(Q, t(basis) %*% data, transpose = TRUE)
    scores <- crossprod(y, x)
  }
  else if (attributes(weightinv)$diagonal == TRUE) {
    V <- t(basis) %*% (diag(weightinv) * basis)
    Q <- chol(V)
    y <- backsolve(Q, diag(p), transpose = TRUE)
    tmp <- t(basis) %*% (diag(weightinv) * data)
    x <- backsolve(Q, tmp, transpose = TRUE)
    scores <- crossprod(y, x)
  }
  else {
    V <- t(basis) %*% weightinv %*% basis
    Q <- chol(V)
    y <- backsolve(Q, diag(p), transpose = TRUE)
    x <- backsolve(Q, t(basis) %*% weightinv %*% data, transpose = TRUE)
    scores <- crossprod(y, x)
  }
  return(t(scores))
}

CalcScores <- Project

#' Field reconstructions from coefficients
#'
#' Given a vector of coefficients for a basis, calculates the field
#'
#' @param coeffs Coefficient vector
#' @param basis Basis matrix
#' @return Reconstructed field.
#'
#' @export
Recon <- function(coeffs, basis){
  if (is.null(dim(basis)[2])){
    q <- 1
  }
  else {
    q <- dim(basis)[2]
  }
  stopifnot(length(coeffs) == q)
  if (is.null(dim(basis)[2])){
    reconstruction <- basis*as.numeric(coeffs)
  }
  else {
    reconstruction <- basis%*%as.numeric(coeffs)
  }
  return(reconstruction)
}

#' Matrix projection
#'
#' Projects a variance matrix onto a given basis
#'
#' @param mat A square matrix to be projected onto the basis
#' @param basis The basis to project with
#' @param weightinv The inverse of positive definite matrix W. If NULL, uses the standard projection, otherwise projects in the norm given by W.
#'
#' @return The projection of the original matrix on the basis.
#'
#' @export
VarProj <- function(mat, basis, weightinv = NULL){
  if (is.null(weightinv)){
    proj <- t(basis) %*% mat %*% basis
  }
  else {
    V <- t(basis) %*% weightinv %*% basis
    Q <- chol(V)
    y <- backsolve(Q, diag(dim(basis)[2]), transpose = TRUE)
    x <- backsolve(Q, t(basis) %*% weightinv, transpose = TRUE)
    comp <- crossprod(y, x)
    proj <- comp %*% mat %*% t(comp)
  }
  return(proj)
}

#' History matching
#'
#' Calculates the implausibility (for the \ell-dimensional field) for a sample of expectations and variances from basis emulators
#'
#' @param DataBasis object containing the basis used in emulation ($tBasis)
#' @param Obs observation vector (length \ell), must be centred
#' @param Expectation a matrix containing emulator expectations, where a given row contains the expectations for the q emulated basis vectors, for some x
#' @param Variance a matrix containing emulatorvariances, where a given row contains the variances for the q emulated basis vectors, for some x
#' @param Error observation error variance matrix
#' @param Disc discrepancy variance matrix
#' @param weightinv if not NULL, the inverse of W = var_err + var_disc, used for projection
#' #' 
#' @return \item{impl}{Vector of implausibilities corresponding to the rows of Expectation and Variance}
#' \item{bound}{The chi-squared bound for an \ell-dimensional field}
#' \item{nroy}{Proportion of parameter settings that are not ruled out, using bound}
#' \item{inNROY}{Vector indicating whether a parameter setting is ruled out}
#'
#' @export
HistoryMatch <- function(DataBasis, Obs, Expectation, Variance, Error, Disc, weightinv = NULL, BasisUncertainty = FALSE){
  q <- dim(Expectation)[2]
  Basis <- DataBasis$tBasis[,1:q]
  l <- dim(Basis)[1]
  stopifnot(q == dim(Variance)[2])
  W <- Error + Disc
  if (is.null(weightinv)){
    weightinv <- GetInverse(W)
  }
  R_W <- ReconError(Obs, Basis, weightinv = weightinv, scale = FALSE)
  # Project observations onto basis if required
  if (length(Obs) == l){
    ObsProj <- CalcScores(Obs, Basis, weightinv = weightinv)
  }
  # Add uncertainty from discarded basis vectors?
  if (BasisUncertainty == TRUE){
    BasisVar <- DiscardedBasisVariance(DataBasis, q, weightinv)
    W <- W + BasisVar
  }
  # Project variance matrices onto basis if required
  if (dim(Disc)[1] == l){
    WProj <- VarProj(W, Basis, weightinv = weightinv)
  }
  nn <- dim(Expectation)[1]
  impl <- as.numeric(mclapply(1:nn, function(i) ImplCoeff(Expectation[i,], Variance[i,], ObsProj, WProj, 0*WProj)))
  impl <- impl + rep(R_W, nn) 
  bound <- qchisq(0.995, l)
  nroy <- sum(impl < bound)/nn
  inNROY <- impl < bound
  return(list(impl = impl, bound = bound, nroy = nroy, inNROY = inNROY))
}

#' Coefficient implausibility
#'
#' Calculates the coefficient implausibility for a single x, given projected quantities
#'
#' @param Expectation length q vector with emulator expectations
#' @param Variance length q vector with emulator variances
#' @param Obs projected observations
#' @param Error projected observation error variance matrix
#' @param Disc projected discrepancy variance matrix
#'  
#' @return The coefficient implausibility (given the matrix used in projection)
#'
#' @export
ImplCoeff <- function(Expectation, Variance, Obs, Error, Disc){
  V <- Error + Disc + diag(Variance)
  Q <- chol(V)
  proj.output <- Expectation
  y <- backsolve(Q, as.vector(Obs - proj.output), transpose = TRUE)
  impl <- crossprod(y,y)
  return(impl)
}





