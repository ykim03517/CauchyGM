#'
#' A Robust Genetic Model Based SNP-set Test using CauchyGM
#'
#' @param g a numeric genotype vector/matrix. Each genotype should be coded as 0, 1, and 2 for AA, Aa, and aa, where A is a major allele and a is a minor allele.  
#' @return MAF for each SNP
#' @author Yeonil Kim 
#' @export


MAF <- function(g) {

    n <- dim(g)[1]
    p <- dim(g)[2]

    if (p == 1) {

        MAF_out <- sum(g)/2 * n

        if (MAF_out >= 0.5) {
            MAF_out = 1 - MAF_out
        }

    } else {

        MAF_out <- apply(g, 2, function(g) sum(g)/(2 * n))
        ind <- which(MAF_out >= 0.5)

        if (sum(ind) != 0) {

            MAF_out[ind] = 1 - MAF_out[ind]
        }

    }
    return(MAF_out)

}


