# titre: pava.r
# auteur: Pierre-O Goffard
# date: 06/08/2023
# Functions to implement the PAVA algorithm and make plots


#' pava: Fit an isotonic regression model using the Pool adjacent violator algorithm
#'
#' @param x vector of predictor
#' @param y variable to predict
#' @param w weight of each observation
#'
#' @return a dataframe with the isotonic fit (iso_fit) and another dataframe with the original data (data)
#' @export
#'
#' @examples
pava <-  function(x, y, w){
  df_init <- data.frame(x= x, y = y, w = w)
  df <- dplyr::summarize( dplyr::group_by(df_init, x), w_sum = sum(w), y = stats::weighted.mean(y, w))
  while(any(c(dplyr::lead(df$y)[-length(df$y)], Inf) + .00001 < df$y)){
    v <- which(c(dplyr::lead(df$y)[-length(df$y)], Inf) < df$y)
    new_y <- ((df$y * df$w_sum)[v + 1]  + (df$y * df$w_sum)[v]) / (df$w_sum[v] + df$w_sum[v+1])
    df$y[v] <- new_y; df$y[v +1] <- new_y
  }
  res <- list(df, df_init); names(res) <- c("iso_fit", "data")
  return(res)
}




