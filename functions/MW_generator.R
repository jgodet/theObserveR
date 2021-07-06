# MW_generator 
#' @author  TaiohY
# 06.07.2021

# Generate matching weights from database with categorical variable
#' @param formula regression formula
#' @param data dataframe imput
#' @param id cluster identification
#' @return dataframe
#' @examples MW_generator(formula = A ~ X1 + X2 + X3, data = DF)
#' @examples MW_generator(formula = A ~ X1 + X2 + X3, data = DF, id = "ID") 


MW_generator <- function(formula, data, id=NULL){
  
  # Requirements
  lapply(c("lme4","varhandle", "parallel"),
         require, 
         character.only = T)
  
  # Settings
  data <- as.data.frame(data)
  VarTrt <- paste(formula[[2]])
  VarAdj <- paste(formula[[3]][2])
  
  # Dummy variables
  df.dummy <- as.data.frame(to.dummy(data[,VarTrt], "tt"))
  df <- cbind(data, df.dummy)
  
  # No random effects
  if(is.null(id)){
    v_formulas <- lapply(paste(names(df.dummy),VarAdj,sep = "~"),
                         formula)
    glm.p <- partial(glm, data = df, family = binomial(link="logit"))
    l_mods <- lapply(v_formulas, glm.p)
  }else {
    # With random effects
    v_formulas <- lapply(paste(names(df.dummy), 
                               paste(VarAdj, "+ (1|", id, ")", collapse = ""),
                               sep = "~"),
                         formula)
    glmer.p <- partial(glmer, data = df, family = binomial(link="logit"))
    l_mods <- mclapply(v_formulas, glmer.p, mc.cores = detectCores()-1)
  }
  
  # Output dataframe
  l_preds <- predict(l_mods, newdata = df, type = "response")
  t_pred <- as.data.frame(Reduce(cbind,l_preds))
  names(t_pred) <- paste("e", c(1:length(t_pred)), sep = "")
  t_pred$e.min <- do.call(pmin, t_pred)
  denominator <- rowSums(t_pred[,c(1:4)]*df.dummy)
  data$MW <- t_pred$e.min/denominator
  
  return(data)
}

