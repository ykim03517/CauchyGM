#'
#' A Robust Genetic Model Based SNP-set Test using CauchyGM
#'
#' @param formula an object of class 'formula': a symbolic description of the NULL model to be fitted. No need to specify outcome (y) in the formula.
#' @param data a list of data frames which contain the outcome (Y) and the covariates (Z) in the model.
#' @param out_type an indicator of the outcome type. 'C' for the continuous outcome and 'D' for the dichotomous outcome.
#' @return an object that includes model parameters and residuals of the null model (there are no genetic effects on the outcome phenotypes).
#' @author Yeonil Kim
#' @import SKAT
#' @import Rcpp
#' @export

Cauchy_Null_model <- function(formula, data, out_type = "D") {

    y <- data$Y
    Z <- data.frame(data$Z)

    if (!is.null(Z)) {

        design.Z <- model.matrix(formula, Z)[, -1]
        Y_design.Z <- data.frame(y, design.Z)
        null_model <- as.formula(paste("y ~ ", paste(colnames(design.Z), collapse = "+")))

        if (out_type == "C") {

            lm_null <- lm(null_model, data = Y_design.Z)

        } else if (out_type == "D") {

            lm_null <- glm(null_model, data = Y_design.Z, family = "binomial")
        }

    } else {

        null_model <- as.formula(paste("y ~ ", 1))
        Y_design.Z <- data.frame(y, matrix(1, nrow = length(y)))

        if (out_type == "C") {

            lm_null <- lm(null_model, data = Y_design.Z)

        } else if (out_type == "D") {

            lm_null <- glm(null_model, data = Y_design.Z, family = "binomial")
        }
    }

    null <- list(fitted = lm_null$fitted, residuals = lm_null$residuals, null.model = lm_null$model)

    return(null)

}


