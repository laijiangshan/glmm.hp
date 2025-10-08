#' Automatically Recommend an Appropriate GLM Family
#'
#' This function assists users in selecting an appropriate \code{family} and link function
#' for generalized linear models (GLM) based on the distributional properties of the response variable.
#' It provides a quick diagnostic summary, distribution plots, and an optional AIC comparison among candidate models.
#'
#' @param y A numeric vector representing the response variable.
#' @param plot Logical; if \code{TRUE}, a histogram and boxplot of \code{y} will be drawn. Default is \code{TRUE}.
#' @param aic_test Logical; if \code{TRUE}, a simple AIC comparison across candidate GLM families will be performed. Default is \code{FALSE}.
#'
#' @details
#' The function inspects the basic characteristics of the response variable, including its range, mean,
#' variance, and whether it contains only integers or proportions. Based on these diagnostics, it suggests
#' one or more candidate GLM families among:
#' \itemize{
#'   \item \code{"gaussian"} – continuous response, may include negative values
#'   \item \code{"poisson"} or \code{"quasipoisson"} – integer count data (possibly overdispersed)
#'   \item \code{"Gamma"} or \code{"inverse.gaussian"} – strictly positive continuous data
#'   \item \code{"binomial"} – proportion or binary data (0–1 range)
#' }
#'
#' The suggested link functions are: \code{"identity"} for Gaussian, \code{"log"} for Poisson/Gamma,
#' and \code{"logit"} for Binomial.
#'
#' @return
#' A list containing:
#' \item{family}{Suggested GLM family}
#' \item{link}{Suggested link function}
#'
#' @examples
#' # Example 1: Continuous response (can be negative)
#' set.seed(123)
#' y1 <- rnorm(100)
#' detect_glm_family(y1)
#'
#' # Example 2: Count data
#' y2 <- rpois(100, lambda = 5)
#' detect_glm_family(y2)
#'
#' # Example 3: Proportion data
#' y3 <- rbeta(100, 2, 5)
#' detect_glm_family(y3)
#'
#' @export
detect_glm_family <- function(y, plot = TRUE, aic_test = FALSE) {
  if (!is.numeric(y)) stop("y must be a numeric vector.")
  y <- y[!is.na(y)]
  
  cat("Basic information:\n")
  cat("Min:", min(y), "\n")
  cat("Max:", max(y), "\n")
  cat("Mean:", mean(y), "\n")
  cat("Variance:", var(y), "\n")
  cat("Contains negative values:", any(y < 0), "\n")
  cat("All integers:", all(abs(y - round(y)) < 1e-6), "\n")
  
  # Plot distribution
  if (plot) {
    par(mfrow = c(1, 2))
    hist(y, main = "Distribution of response variable", xlab = "y", 
         col = "skyblue", border = "white")
    boxplot(y, main = "Boxplot", horizontal = TRUE, col = "lightgreen")
    par(mfrow = c(1, 1))
  }
  
  # Determine type
  n_int <- all(abs(y - round(y)) < 1e-6)
  has_neg <- any(y < 0)
  is_prop <- all(y >= 0 & y <= 1)
  
  suggestion <- ""
  link <- ""
  
  if (is_prop) {
    suggestion <- "binomial (or beta if continuous proportions)"
    link <- "logit"
  } else if (n_int && all(y >= 0)) {
    mean_y <- mean(y)
    var_y <- var(y)
    if (var_y > 1.5 * mean_y) {
      suggestion <- "quasipoisson or negative binomial"
    } else {
      suggestion <- "poisson"
    }
    link <- "log"
  } else if (!has_neg && min(y) > 0) {
    mean_y <- mean(y)
    var_y <- var(y)
    if (abs(var_y - mean_y^2) < abs(var_y - mean_y^3)) {
      suggestion <- "Gamma"
      link <- "log"
    } else {
      suggestion <- "inverse.gaussian"
      link <- "1/mu^2"
    }
  } else {
    suggestion <- "gaussian"
    link <- "identity"
  }
  
  cat("\nRecommended GLM family:", suggestion, "\n")
  cat("Suggested link function:", link, "\n")
  
  # Optional AIC comparison
  if (aic_test && length(y) > 10) {
    x <- seq_along(y)
    test_models <- list(
      gaussian = try(glm(y ~ x, family = gaussian), silent = TRUE),
      poisson = try(glm(y ~ x, family = poisson), silent = TRUE),
      quasipoisson = try(glm(y ~ x, family = quasipoisson), silent = TRUE),
      Gamma = try(glm(y ~ x, family = Gamma(link = "log")), silent = TRUE),
      binomial = try(glm(y ~ x, family = binomial), silent = TRUE)
    )
    valid_models <- test_models[!sapply(test_models, inherits, "try-error")]
    if (length(valid_models) > 1) {
      cat("\nAIC comparison among candidate families:\n")
      print(sapply(valid_models, AIC))
    }
  }
  
  cat("\nExample formula:\n")
  cat("glm(response ~ predictor1 + predictor2, family =", suggestion, ", data = your_data)\n")
  
  invisible(list(family = suggestion, link = link))
}
