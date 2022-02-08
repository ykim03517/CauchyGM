#'
#' A Robust Genetic Model Based SNP-set Test: CauchyGM
#'
#' @param y a numeric vector of the outcome phenotype (continuous or binary phenotype). 
#' @param G a numeric genotype vector/matrix. 
#' @param Z a numeric vector/matrix for covariate(s). Default=NULL. 
#' @param null_obj an object that includes model parameters and residuals of the null model (i.e., there are no genetic effects on the outcome phenotypes). Please run CauchyGM_Null to get this object. 
#' @param out_type an indicator of the outcome type. 'C' for the continuous outcome and 'D' for the dichotomous outcome. 
#' @return A list of numeric vectors/matrices that contains score statistics (U) for each subject and SNP and that include the variance of the score statistics.
#' @author Yeonil Kim


eff_Score <- function(y, G, Z = NULL, null_obj, out_type = "D") {

    G = as.matrix(G)
    m = dim(G)[2]
    n = length(y)
    Z = as.matrix(Z)

    if (out_type == "D") {

        p_hat = null_obj$fitted
        Z = cbind(rep(1, n), Z)  # adding one column in covariate matrix for intercept (n x q+1) 
        p = dim(Z)[2]
        resid = (y - p_hat)

        y = array(rep(y, m), c(n, m))
        p_hat = array(rep(p_hat, m), c(n, m))
        pp_hat = p_hat * (1 - p_hat)

        if (!is.null(Z)) {
            # score function for beta1
            s_b1_i = array(rep(resid, m), c(n, m)) * G  # n x m
            s_b1 = colSums(s_b1_i)

            # score function for sigma = (b0 and delta)
            s_sigma_i = Z * array(rep(resid, p), c(n, p))  # n x q+1

            # I(beta1, sigma = c(beta0, delta))
            I.b1.sig = -t(t(Z) %*% (pp_hat * G))  # m x (q+1)

            # I(sigma, sigma)
            s_sigma = Z * array(rep(pp_hat[, 1], p), c(n, p))  # n x q+1

            I.sig.sig = -t(Z) %*% s_sigma  # q+1 x q+1

            # individual score function for beta1
            U_nus = t(I.b1.sig %*% ginv(I.sig.sig) %*% t(s_sigma_i))  # n x m
            U = s_b1_i - U_nus  # n x m
            V = colSums(U^2)

        } else {
            ## no covariates
            Z = array(rep(1, n), c(n, 1))  # adding one column for intercept (n x 1) 
            dim_Z = 1  # dim Z = 1

            # score function for beta1
            s_b1_i = array(rep(resid, m), c(n, m)) * G  # n x m
            s_b1 = colSums(s_b1_i)

            # score function for sigma = (b0 and delta)
            s_sigma_i = array(resid, c(n, 1))  # n x 1

            # I(beta1, sigma = c(beta0, delta))
            I.b1.sig = -colSums(pp_hat * G)  # m x 1

            # I(sigma, sigma)
            I.sig.sig = -sum(pp_hat[, 1])  # 1 x 1 scalar 

            # individual score function for beta1
            U_nus = t((I.b1.sig * I.sig.sig^(-1)) %*% t(s_sigma_i))  # n x m
            U = s_b1_i - U_nus  # n x m
            V = colSums(U^2)

        }



    } else if (out_type == "C") {

        Z = cbind(rep(1, n), Z)  # adding one column in covariate matrix for intercept (n x q+1) 
        p = dim(Z)[2]
        resid = null_obj$residuals
        var = (1/n) * sum(resid^2)
        y = array(rep(y, m), c(n, m))

        if (!is.null(Z)) {

            U = array(rep(resid, m), c(n, m)) * G  #n x m
            cov_Z = t(Z) %*% Z

            var2 = numeric(m)
            for (i in 1:m) {
                G_temp = array(rep(G[, i], p), c(n, p))
                GZ = G_temp * Z
                GZ = colSums(GZ)
                var2[i] = t(GZ) %*% solve(cov_Z) %*% GZ
            }

            V = var * (colSums(G^2) - var2)

        } else {

            # score for beta1
            s_b1_i = array(rep(resid, m), c(n, m)) * G  # n x m
            s_b1 = colSums(s_b1_i)  # 1 x m 

            # score for sigma
            s_b0_i = resid
            s_b0 = sum(resid)
            s_b0_sq_i = resid^2
            s_b0_sq = sum(resid^2)

            s_sigma_i = cbind((1/var) * s_b0_i, (-1/2) * (1/var) + (1/(2 * var^2)) * s_b0_sq_i)  # n x 2

            # I(beta1, sigma=c(beta0, var))
            I.1.1 = -(1/var) * colSums(G)
            I.1.2 = -(1/var^2) * s_b1

            I.b1.sig = cbind(I.1.1, I.1.2)  # m x 2

            # I(sigma, sigma)
            I.2.11 = -n/(var)
            I.2.12 = -1/(var^2) * s_b0
            I.2.22 = n/(2 * var^2) - 1/(var^3) * s_b0_sq

            I.sig.sig = matrix(c(I.2.11, I.2.12, I.2.12, I.2.22), nrow = 2)

            U_nus = t(I.b1.sig %*% ginv(I.sig.sig) %*% t(s_sigma_i))  # n x m
            U = s_b1_i - U_nus
            V = colSums(U^2)

        }

    }

    return(list(U = U, V = V))

}
