getLL <- function(object, ...){
  models <- rlang::list2(...)
  
  if(length(models)>0) { # if several models are provided
    modNames <- all.vars(match.call()) # store the names of the models given as arguments
    
    # include "object" in "models"
    modcopy <- list()
    modcopy[[1]] <- object
    for(i in 1:length(models))
      modcopy[[i+1]] <- models[[i]]
    models <- modcopy
    
    for(i in 1:length(models)){
      if(!is.null(models[[i]]$modelName)) modNames[i] <- models[[i]]$modelName
    }
    
    if(any(!unlist(lapply(models,function(x) "momentuHMM" %in% class(x))))) stop("all models must be momentuHMM objects")
    
    LLs <- data.frame(model = modNames, log_lik = NA)
    for (i in c(1:length(models))){
      if (is.null(models[[i]])){LLs[i, "log_lik"] <- NA
      }else{
      LLs[i, "log_lik"] <- -models[[i]]$mod$minimum
    }
    arrange(LLs, log_lik)
    }
  }else{
    if (is.null(object)){
      LLs <- NA
      }else{
      LLs <- -object$mod$minimum
      }
  }
  return(LLs)
}