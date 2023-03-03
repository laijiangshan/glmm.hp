#' Plot for a \code{\link{glmm.hp}} object
#'
#' @param x A \code{\link{glmm.hp}} object.
#' @param  plot.perc Logical;if TRUE, the bar plot (based on ggplot2 package) of the percentage to individual effects of variables or groups towards total explained variation, the default is FALSE to show plot with original individual effects.
#' @param  n Integer; which marginal R2 in output of r.squaredGLMM to plot. 
#' @param ... unused
#' @return a ggplot object
#' @author {Jiangshan Lai} \email{lai@njfu.edu.cn}
#' @export
#' @examples
#' library(MuMIn)
#' library(lme4)
#' mod1 <- lmer(Sepal.Length ~ Petal.Length + Petal.Width +(1 | Species), data = iris)
#' plot(glmm.hp(mod1))


plot.glmmhp <- function(x, plot.perc = FALSE, n = 1,...){
  if (!inherits(x, "glmmhp")){
    stop("x should be the output of glmm.hp()")
  }

  if (plot.perc){
    tips3 = data.frame(variable = rownames(x[[n+1]]),
                       value = as.numeric(x[[n+1]][, "I.perc(%)"]))
    gg = ggplot2::ggplot(tips3, ggplot2::aes(x = stats::reorder(variable,-value), y = value)) + ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Variables", y = "% Individual effect to Rsquare (%I)") + ggplot2::theme(axis.text = element_text(size = 10)) + ggplot2::theme(axis.title = element_text(size = 13))+ggplot2::labs(title =paste(names(x)[n+1]))
  } else {
    tips2 = data.frame(variable = rownames(x[[n+1]]),
                       value = as.numeric(x[[n+1]][, "Individual"]))
    gg = ggplot2::ggplot(tips2, ggplot2::aes(x = stats::reorder(variable,-value), y = value)) + ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Variables", y = "Individual effect") + ggplot2::theme(axis.text = element_text(size = 10)) + ggplot2::theme(axis.title = element_text(size = 13))+ggplot2::labs(title =paste(names(x)[n+1]))
  }

  return(gg)

}
