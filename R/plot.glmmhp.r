#' Plot for a \code{\link{glmm.hp}} object
#'
#' @param x A \code{\link{glmm.hp}} object.
#' @param  plot.perc Logical;if TRUE, the bar plot (based on ggplot2 package) of the percentage to individual effects of variables or groups towards total explained variation, the default is FALSE to show plot with original individual effects.
#' @param ... unused
#' @return a ggplot object
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @export
#' @examples
#'mod <- lmer(Biomass ~ Year + Temperature + Precipitation +SpeciesDiversity+(1 | Population),data = biomass)
#'plot(glmm.hp(mod))


plot.glmmhp <- function(x, plot.perc = FALSE, ...){
  if (class(x) != "glmmhp"){
    stop("x should be the output of glmm.hp()")
  }

  if (plot.perc){
    tips3 = data.frame(variable = rownames(x$Hier.part),
                       value = as.numeric(x$Hier.part[, "I.perc(%)"]))
    gg = ggplot2::ggplot(tips3, ggplot2::aes(x = stats::reorder(variable,-value), y = value)) + ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Variables", y = "% Individual effect to Rsquare (%I)") + ggplot2::theme(axis.text = element_text(size = 10)) + ggplot2::theme(axis.title = element_text(size = 13))
  } else {
    tips2 = data.frame(variable = rownames(x$Hier.part),
                       value = as.numeric(x$Hier.part[, "Individual"]))
    gg = ggplot2::ggplot(tips2, ggplot2::aes(x = stats::reorder(variable,-value), y = value)) + ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Variables", y = "Individual effect") + ggplot2::theme(axis.text = element_text(size = 10)) + ggplot2::theme(axis.title = element_text(size = 13))
  }

  return(gg)

}
