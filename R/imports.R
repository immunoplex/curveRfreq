utils::globalVariables(c(
  "converged", "value",
  "yhat", "model", "estimate", "conf.low", "conf.high"
))

#' @importFrom stats AIC BIC as.formula coef df.residual formula median nls nls.control residuals sd setNames vcov qt predict fitted
#' @importFrom utils setTxtProgressBar txtProgressBar globalVariables
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @importFrom data.table :=
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_hline geom_linerange facet_wrap coord_flip scale_y_continuous expansion labs theme_bw theme element_text xlab ylab
#' @importFrom patchwork plot_annotation
#' @importFrom reshape2 melt
"_PACKAGE"
