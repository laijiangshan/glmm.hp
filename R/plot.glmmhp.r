#' Plot for a \code{\link{glmm.hp}} object
#'
#' @param x A \code{\link{glmm.hp}} object.
#' @param  plot.perc Logical;if TRUE, the bar plot (based on ggplot2 package) of the percentage to individual effects of variables or groups towards total explained variation, the default is FALSE to show plot with original individual effects.
#' @param color Color of variables.
#' @param  n Integer; which marginal R2 in output of r.squaredGLMM to plot. 
#' @param  dig Integer; number of decimal places in Venn diagram. 
#' @param ... unused
#' @return a ggplot object
#' @author {Jiangshan Lai} \email{lai@njfu.edu.cn}
#' @export
#' @examples
#' library(MuMIn)
#' library(lme4)
#' mod1 <- lmer(Sepal.Length ~ Petal.Length + Petal.Width +(1 | Species), data = iris)
#' a <- glmm.hp(mod1)
#' plot(a)
#' mod3 <- lm(Sepal.Length ~ Petal.Length+Petal.Width,data = iris)
#' plot(glmm.hp(mod3,type="R2"))
#' plot(glmm.hp(mod3,commonality=TRUE),color = c("#8DD3C7", "#FFFFB3"))

plot.glmmhp <- function(x, plot.perc = FALSE, color = NULL,n = 1,dig = 4,...){
  if (!inherits(x, "glmmhp")){
    stop("x should be the output of glmm.hp()")
  }
if(x$type=="hierarchical.partitioning")
  {if (plot.perc){
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
if(x$type=="commonality.analysis")
{ 
  Var.part <- as.data.frame(x[[n+1]])
  #Var.part <- Var.part[which(Var.part$Fractions >=cutoff), ]
  Var.part$Fractions <- round(Var.part$Fractions,dig)
  nvar <- length(x$variables)
  Constrained <- Var.part$Fractions[2^nvar]
  if (!nvar%in% 2:4)
    stop("Venn diagram supports only 2-4 variables")
  else if (nvar == 2)
	Var <- Var.part$Fractions[1:3]
  else if (nvar == 3)
    Var <- Var.part$Fractions[c(1:4, 6, 5, 7)]
  else if (nvar == 4)
    Var <- Var.part$Fractions[c(1:5, 7, 6, 8:10, 12, 11, 14, 13, 15)]
    vegan::showvarparts(part = nvar, bg = color, Xnames = x$variables, labels = as.character(c(Var, 1-Constrained)))
}

}
