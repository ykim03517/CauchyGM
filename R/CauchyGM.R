#'
#' A Robust Genetic Model Based SNP-set Test using CauchyGM
#'
#' @param data a list of data frames which contain the outcome (Y), the genotype (X_add) and the covariates (Z) in the model.
#' @param method a character value to specify a test procedure (default='CauchyGM'). 'CauchyGM' represents a robust genetic model based SNP-set test. 'CauchyGM-O' represents a integreated test procedure that combines CauchyGM with SKAT and burdent test.
#' @param out_type an indicator of the outcome type. 'C' for the continuous outcome and 'D' for the dichotomous outcome.
#' @param int a numeric value to set an interval of a equal sized grid between 0 and 1. This grid indicates genotypes (0=AA, Aa = 0 < c < 1, 1=aa)
#' @param weights.beta a vector of parameters of beta weights (default= c(1,1))
#' @param Cauchy_Null an output object of the Cauchy_Null_model function
#' @param SKAT_Null an output object of the SKAT_Null (continuous) or SKAT_Null_Model_MomentAdjust (binary) function
#' @return p-value of CauchyGM or CauchyGM-O
#' @author Yeonil Kim
#' @useDynLib CauchyGM
#' @importFrom Rcpp sourceCpp
#' @import MASS
#' @import SKAT
#' @export



CauchyGM <- function(data, method = "CauchyGM", out_type = "D", int = 0.5, weights.beta = c(1, 1), Cauchy_Null,
    SKAT_Null) {

    seqC <- seq(0, 1, by = int)
    y <- data$Y
    g <- data$X_add

    covar <- Cauchy_Null$model[, -1]

    n <- length(y)
    m <- length(seqC)

    ## Dimension of genotype data and creating grid for genotype data
    dim = dim(g)[2]

    if (is.null(dim)) {
        p <- 1
        grid.g = array(0, c(n, m, p))

        grid.g[, , 1] = gridG(g, seqC)

    } else {
        p <- dim
        grid.g = array(0, c(n, m, p))

        for (j in 1:p) {
            grid.g[, , j] = gridG(g[, j], seqC)
        }
    }

    L <- c()
    ind_list <- list()
    for (c in 1:m) {

        MAF_ind <- MAF(2 * grid.g[, c, 1:p])
        ind_rm <- which((MAF_ind < 0.01) == TRUE)
        ind_list[[c]] <- ind_rm
        L <- c(L, c(length(ind_rm)))

    }

    ind <- array(NA, c(p, m))
    pval <- c()
    for (c in 1:m) {

        if (L[c] == p) {

            burden_pval = SKAT(as.matrix(2 * grid.g[, c, ]), SKAT_Null, weights.beta = c(1, 25), r.corr = 1)$p.value
            pval = c(pval, burden_pval)

            ind[, c] <- 0

        } else if (L[c] < p) {

            if (L[c] == 0) {

                U_V = eff_Score(y, grid.g[, c, ], Z = covar, Cauchy_Null, out_type = out_type)
                score = colSums(U_V[[1]])^2/U_V[[2]]
                pval = c(pval, 1 - pchisq(score, 1))

                ind[, c] <- 1:p

            } else {

                grid.g_common = grid.g[, c, -ind_list[[c]]]
                grid.g_rare = grid.g[, c, ind_list[[c]]]

                U_V = eff_Score(y, grid.g_common, Z = covar, Cauchy_Null, out_type = out_type)
                score = colSums(U_V[[1]])^2/U_V[[2]]
                pval = c(pval, 1 - pchisq(score, 1))

                burden_pval = SKAT(as.matrix(2 * grid.g_rare), SKAT_Null, weights.beta = c(1, 25), r.corr = 1)$p.value
                pval = c(pval, burden_pval)
            }

            ind[ind_list[[c]], c] <- 0
            ind[-ind_list[[c]], c] <- c(1:p)[-ind_list[[c]]]

        }

    }

    MAF_out <- c()
    MAF_add <- MAF(g)

    for (c in 1:m) {

        if (L[c] != 0) {
            MAF_temp <- c(MAF_add[which(ind[, c] != 0)], mean(MAF_add[which(ind[, c] == 0)]))
            MAF_out <- c(MAF_out, MAF_temp)

        } else {
            MAF_temp <- c(MAF_add[which(ind[, c] != 0)])
            MAF_out <- c(MAF_out, MAF_temp)
        }
    }

    MAF_out = MAF_out[pval != 1]
    pval = pval[pval != 1]

    w = (dbeta(MAF_out, weights.beta[1], weights.beta[2]) * sqrt(MAF_out * (1 - MAF_out)))^2
    w = w/sum(w)
    T = sum((w) * tan((0.5 - pval) * pi))
    pval_proposed = 1/2 - atan(T)/pi


    if (method == "CauchyGM") {

        return(list(pval = pval_proposed))

    } else if (method == "CauchyGM-O") {

        if (outcome == "D") {
            pval_SKAT = SKATBinary(as.matrix(g), SKAT_Null, weights.beta = c(1, 1))$p.value
            pval_burden = SKATBinary(as.matrix(g), SKAT_Null, weights.beta = c(1, 1), r.corr = 1)$p.value
        } else if (outcome == "C") {
            pval_SKAT = SKAT(as.matrix(g), SKAT_Null, weights.beta = c(1, 1))$p.value
            pval_burden = SKAaT(as.matrix(g), SKAT_Null, weights.beta = c(1, 1), r.corr = 1)$p.value
        }

        pval = c(pval_proposed, pval_SKAT, pval_burden)

        T = sum((1/length(pval)) * tan((0.5 - pval) * pi))
        pval_comb = 1/2 - atan(T)/pi

        return(list(pval = pval_comb))

    }

}



# design.Z <- model.matrix(null_model, data$Z)[, -1] if (outcome == 'D') { SKAT_Null <-
# SKAT_Null_Model_MomentAdjust(data$Y ~ design.Z, type.Resampling = 'bootstrap.fast') } else if
# (outcome == 'C') { # SKAT_Null <- SKAT_Null_Model(data$Y ~ design.Z, out_type = outcome) SKAT_Null <-
# SKAT_Null_Model(data$Y ~ design.Z, out_type = outcome, type.Resampling = 'bootstrap.fast') }








