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
  lapply(c("lme4","varhandle","purrr", "parallel"),
         require, 
         character.only = T)
  
  # Settings
  data <- as.data.frame(data)
  VarTrt <- paste(formula[[2]])
  VarAdj <- paste(formula[[3]][2])
  
  # Dummy variables
  df.dummy <- as.data.frame(to.dummy(data[,VarTrt], "tt"))
  dummy.names <- names(df.dummy)
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
  
  NIPls <- lapply(c(1:length(l_mods)),FUN = function(i){
    unique(as.numeric(row.names(ranef(l_mods[[i]])$NIP)))
  })
  
  length_NIPls <- lapply(c(1:length(NIPls)), FUN = function(i){
    length(NIPls[[i]])
  })
  
  indice <- min(which(min(unique(unlist(length_NIPls)))==unlist(length_NIPls)))
  newdata <-  df %>%
      filter(NIP %in% NIPls[[indice]])
  l_preds <- predict(l_mods, newdata, type = "response")
  t_pred <- as.data.frame(Reduce(cbind,l_preds))
  names(t_pred) <- paste("e", c(1:length(t_pred)), sep = "")
  t_pred$e.min <- do.call(pmin, t_pred)
  denominator <- rowSums(t_pred[,c(1:4)]*newdata[,dummy.names])
  newdata$MW <- t_pred$e.min/denominator
  
  return(newdata)
}

