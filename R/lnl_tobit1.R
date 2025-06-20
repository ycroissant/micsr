lnl_tp <- function(param, X, y, wt, Z = NULL, scedas = c("exp", "pnorm"), sum = TRUE, gradient = FALSE,
                   hessian = FALSE, left = 0, right = Inf, sample = "censored"){

    CENS <- sample == "censored"
    TRUNC <- ! CENS
    LEFT <- ! is.infinite(left) & is.infinite(right)
    RIGHT <- is.infinite(left) & ! is.infinite(right)
    TWO <- ! is.infinite(left) & ! is.infinite(right)
    if (TWO & TRUNC) stop("the two-limits truncated model is not supported")
    Ia <- y <= left
    Ib <- y >= right
    Io <- (y > left & y < right)

    K <- ncol(X)
    N <- length(y)
    beta <- param[1:K]
    bX <- drop(X %*% beta)
    sigo <- param[K + 1]

    .scedas <- match.arg(scedas)
    if (.scedas == "exp"){
        fh <- function(x) exp(x)
        gh <- function(x) exp(x)
        hh <- function(x) exp(x)
    }
    if (.scedas == "pnorm"){
        fh <- function(x) pnorm(x)
        gh <- function(x) dnorm(x)
        hh <- function(x) -x * dnorm(x)
    }

    if (is.null(Z)){
        heter <- FALSE
        sig <- sigo
    }
    else{
        heter <- TRUE
        J <- ncol(Z)
        gamma <- param[(K + 2) : (K + 1 + J)]
        gZ <- as.numeric(Z %*% gamma)
        sig <- sigo * fh(gZ)
    }
    
    za <- (left  - bX) / sig
    zb <- (right - bX) / sig
    ze <- (    y - bX) / sig
    lnl_com <-  Io * (- log(sig) - 0.5 * log(2 * pi) - 1 / 2 * ze ^ 2)

    if (CENS){
        lnl_spec <- 0
        if (LEFT + TWO)  lnl_spec <-            Ia * pnorm(za, log.p = TRUE)
        if (RIGHT + TWO) lnl_spec <- lnl_spec + Ib * pnorm(- zb, log.p = TRUE)
    }
    if (TRUNC){
        if (LEFT)  lnl_spec <- - Io * pnorm(- za, log.p = TRUE)
        if (RIGHT) lnl_spec <- - Io * pnorm(  zb, log.p = TRUE)
        if (TWO)   lnl_spec <- -  Io * log(pnorm(zb) - pnorm(za))
    }
    lnl <- lnl_com + lnl_spec
    if (sum) lnl <- sum(lnl * wt)
    if (gradient){
        g_com_beta <- Io * ze
        g_com_sig <-  Io * (ze ^ 2 - 1)
        if (CENS){
            g_spec_beta <- g_spec_sig <- 0
            if (LEFT | TWO){
                g_spec_beta <-  - Ia * mills(za)
                g_spec_sig <-   - Ia * mills(za) * za
            }
            if(RIGHT | TWO){
                g_spec_beta <-   g_spec_beta +  Ib * mills(- zb)
                g_spec_sig <-    g_spec_sig +  Ib * mills(-zb) * zb
            }
        }
        if (TRUNC){
            if (LEFT){
                g_spec_beta <- -  Io * mills(- za)
                g_spec_sig <-  -  Io * mills(- za) * za
            }
            if (RIGHT){
                g_spec_beta <-  Io * mills(zb)
                g_spec_sig <-   Io * mills(zb) * zb
            }
            if (TWO){
                Dphi <- pnorm(zb) - pnorm(za)
                g_spec_beta <-  Io * (mills(zb) * pnorm(zb) - mills(za) * pnorm(za)) / Dphi
                g_spec_sig <-  Io * (mills(zb) * pnorm(zb) * zb - mills(za) * pnorm(za) * za) / Dphi
            }
        }
        if (heter){
            grad <- cbind((g_com_beta + g_spec_beta) * X,
            (g_com_sig  + g_spec_sig) * cbind(fh(gZ), sigo * gh(gZ) * Z)) / sig
        }
        else{
            grad <- cbind((g_com_beta + g_spec_beta) * X,
            (g_com_sig  + g_spec_sig) ) / sig
        }
        if (sum) grad <- apply(wt * grad, 2, sum)
        attr(lnl, "gradient") <- grad
    }
    
    if (hessian){
        h_com_beta_beta <- - Io
        h_com_beta_sig <- - 2 * Io * ze
        h_com_sig_sig <- Io * (- 3 * ze ^ 2 + 1)
        Pa <- pnorm(za)
        Pb <- pnorm(- zb)
        Pe <- 1 - Pa - Pb

        E_h_com_beta_beta <- - Pe
        E_h_com_beta_sig <- - 2 * (dnorm(za) - dnorm(zb))
        E_h_com_sig_sig <- - 2 * Pe
        if (! LEFT) E_h_com_sig_sig <- E_h_com_sig_sig + (zb * dnorm(zb))
        if (! RIGHT) E_h_com_sig_sig <- E_h_com_sig_sig - (za * dnorm(za))
            
        if (CENS){
            h_spec_beta_beta <- h_spec_beta_sig <- h_spec_sig_sig <- 0
            if (LEFT | TWO){
                h_spec_beta_beta <-  Ia * mills(za, 1)
                h_spec_beta_sig <-    Ia * (mills(za) + mills(za, 1) * za)
                h_spec_sig_sig <-  Ia * (za ^ 2 * mills(za, 1) + 2 * za * mills(za))
            }
            if (RIGHT | TWO){
                h_spec_beta_beta <- h_spec_beta_beta + Ib * mills(-zb, 1)
                h_spec_beta_sig <-  h_spec_beta_sig  + Ib * ( - mills(-zb) + mills(- zb, 1) * zb)
                h_spec_sig_sig <-   h_spec_sig_sig   + Ib * (zb ^ 2 * mills(- zb, 1) - 2 * zb * mills(- zb))
            }
        }
        if (TRUNC){
            if (LEFT){
                h_spec_beta_beta <- - Io *  mills( - za, 1)
                h_spec_beta_sig <-    Io * (mills(- za) - za * mills(- za, 1))
                h_spec_sig_sig <-     Io * (2 * mills(- za) * za - mills(- za, 1) * za ^ 2)
            }
            if (RIGHT){
                h_spec_beta_beta <- - Io * mills(zb, 1)
                h_spec_beta_sig <-    Io * (- mills(zb) - zb * mills(zb, 1))
                h_spec_sig_sig <-     Io * (- 2 * mills(zb) * zb - mills(zb, 1) * zb ^ 2)
            }
            if (TWO){
                DELTA <- pnorm(zb) - pnorm(za)
                A <- mills(zb) * pnorm(zb) - mills(za) * pnorm(za)
                B <- mills(zb) * pnorm(zb) * zb - mills(za) * pnorm(za) * za
                MILLS <- A / DELTA
                A_sig <- - pnorm(zb) * (mills(zb) ^ 2 + mills(zb, 1)) * zb + pnorm(za) * (mills(za) ^ 2 + mills(za, 1)) * za
                B_sig <- - pnorm(zb) * zb * (mills(zb, 1) * zb + mills(zb) ^ 2 * zb + mills(zb)) +
                     pnorm(za) * za * (mills(za, 1) * za + mills(za) ^ 2 * za + mills(za))
                D_sig <- - mills(zb) * pnorm(zb) * zb + mills(za) * pnorm(za) * za
                g_spec_beta <- - Io * MILLS
                g_spec_sig <- - Io * (mills(zb) * pnorm(zb) * zb - mills(za) * pnorm(za) * za) / Dphi
                DD <-  (pnorm(zb) * (mills(zb, 1) + mills(zb) ^ 2) - pnorm(za) * (mills(za, 1) + mills(za) ^ 2)) / Dphi -
                    (mills(zb) * pnorm(zb) - mills(za) * pnorm(za)) ^ 2 / (pnorm(zb) - pnorm(za)) ^ 2
                h_spec_beta_beta <- -Io * DD
                h_spec_beta_sig <-   -Io * (MILLS - (A_sig - MILLS * D_sig) / DELTA)
                h_spec_sig_sig <- -Io * (B / DELTA - (B_sig * DELTA - B * D_sig) / DELTA ^ 2)
            }
        }
        if (heter){
            H_bb <- crossprod(wt * (h_com_beta_beta + h_spec_beta_beta) * X / sig, X / sig)
            H_bs <- crossprod(wt * (h_com_beta_sig  + h_spec_beta_sig ) * X / sig, cbind(fh(gZ), sigo * gh(gZ) * Z) / sig)
            h_ss <- (h_com_sig_sig + h_spec_sig_sig)
            A <- sum(wt * h_ss * fh(gZ) ^ 2 / sig ^ 2)
            g_s <- (g_com_sig + g_spec_sig) / sig
            B <- apply(wt * gh(gZ) * ( h_ss  / sig ^ 2 * sigo* fh(gZ) + g_s) * Z, 2, sum)
            D <- crossprod(wt * (h_ss * sigo ^ 2 * gh(gZ) ^ 2 / sig ^ 2 + sigo * g_s * hh(gZ)) * Z, Z)
            H_ss <- rbind(c(A, B), cbind(B, D))
            attr(lnl, "hessian") <- rbind(cbind(H_bb, H_bs),
                                          cbind(t(H_bs), H_ss))
        }
        else{
            H_bb <- crossprod(wt * (h_com_beta_beta + h_spec_beta_beta) * X, X)
            H_bs <- apply(wt * (h_com_beta_sig + h_spec_beta_sig) * X, 2, sum)
            H_ss <- sum(wt * (h_com_sig_sig + h_spec_sig_sig))
            attr(lnl, "hessian") <- rbind(cbind(H_bb, H_bs),
                                          c(H_bs, H_ss)) / sig ^ 2
            
            I_bb <- crossprod(wt * (E_h_com_beta_beta + h_spec_beta_beta) * X, X)
            I_bs <- apply(wt * (E_h_com_beta_sig + h_spec_beta_sig) * X, 2, sum)
            I_ss <- sum(wt * (E_h_com_sig_sig + h_spec_sig_sig))

            Im1 <- - rbind(cbind(I_bb, I_bs), c(I_bs, I_ss)) / sig ^ 2
            if (! is_definite_positive(Im1)){
                attr(lnl, "info") <- NA
            } else {
                attr(lnl, "info") <- solve(Im1)
            }
            ## attr(lnl, "info") <- solve(- rbind(cbind(I_bb, I_bs),
            ##                                    c(I_bs, I_ss)) / sig ^ 2)
        }
    }
    lnl
}

# Olsen reparametrization

lnl_tp_olsen <- function(param, X, y, wt, sum = TRUE, gradient = FALSE, hessian = FALSE,
                         left = 0, right = Inf, sample = "censored"){

    CENS <- sample == "censored"
    TRUNC <- ! CENS
    LEFT <- ! is.infinite(left) & is.infinite(right)
    RIGHT <- is.infinite(left) & ! is.infinite(right)
    TWO <- ! is.infinite(left) & ! is.infinite(right)
    if (TWO & TRUNC) stop("the two-limits truncated model is not supported")
    Ia <- y <= left
    Ib <- y >= right
    Io <- (y > left & y < right)

    K <- length(param) - 1
    N <- length(y)
    beta <- param[1:K]
    sig <- param[K + 1]

    bX <- as.numeric(X %*% beta)
    za <- (sig * left - bX)
    zb <- (sig * right - bX)
    ze <- (sig * y - bX)

    lnl_com <-  Io * ( log(sig) - 0.5 * log(2 * pi) - 1 / 2 * ze ^ 2)

    if (CENS){
        lnl_spec <- 0
        if (LEFT + TWO)  lnl_spec <-            Ia * pnorm(za, log.p = TRUE)
        if (RIGHT + TWO) lnl_spec <- lnl_spec + Ib * pnorm(- zb, log.p = TRUE)
    }
    if (TRUNC){
        if (LEFT)  lnl_spec <- - Io * pnorm(- za, log.p = TRUE)
        if (RIGHT) lnl_spec <- - Io * pnorm(  zb, log.p = TRUE)
        if (TWO)   lnl_spec <- - Io * log(pnorm(zb) - pnorm(za))
    }
    lnl <- lnl_com + lnl_spec
    if (sum) lnl <- sum(wt * lnl)
    
    if (gradient){
        g_com_beta <- Io * ze
        g_com_sig <-  Io * (1 / sig  - ze * y)
        if (CENS){
            g_spec_beta <- g_spec_sig <- 0
            if (LEFT | TWO){
                g_spec_beta <-  - Ia * mills(za)
                g_spec_sig <-     Ia * mills(za) * left
            }
            if(RIGHT | TWO){
                g_spec_beta <-   g_spec_beta +  Ib * mills(- zb)
                g_spec_sig <-    g_spec_sig  -  Ib * mills(-zb) * right
            }
        }
        if (TRUNC){
            if (LEFT){
                g_spec_beta <- -  Io * mills(- za)
                g_spec_sig <-     Io * mills(- za) * left
            }
            if (RIGHT){
                g_spec_beta <-    Io * mills(zb)
                g_spec_sig <-   - Io * mills(zb) * right
            }
            if (TWO){
                Dphi <- pnorm(zb) - pnorm(za)
                g_spec_beta <-   Io * (mills(zb) * pnorm(zb) - mills(za) * pnorm(za)) / Dphi
                g_spec_sig <-  - Io * (mills(zb) * pnorm(zb) * right - mills(za) * pnorm(za) * left) / Dphi
            }
        }
        grad <- cbind((g_com_beta +  g_spec_beta) * X,
                      g_com_sig + g_spec_sig)
        if (sum) grad <- apply(wt * grad, 2, sum)
        attr(lnl, "gradient") <- grad
    }
    
    if (hessian){
        h_com_beta_beta <- - Io
        h_com_beta_sig <-  Io * y
        h_com_sig_sig <- - Io * (1 / sig ^ 2 + y ^ 2)
        if (CENS){
            h_spec_beta_beta <- h_spec_beta_sig <- h_spec_sig_sig <- 0
            if (LEFT | TWO){
                h_spec_beta_beta <-  Ia * mills(za, 1)
                h_spec_beta_sig <-   - Ia * mills(za, 1) * left
                h_spec_sig_sig <-  Ia * mills(za, 1) * left ^ 2
            }
            if (RIGHT | TWO){
                h_spec_beta_beta <- h_spec_beta_beta + Ib * mills(-zb, 1)
                h_spec_beta_sig <-  h_spec_beta_sig  - Ib * mills(- zb, 1) * right
                h_spec_sig_sig <-   h_spec_sig_sig   + Ib * mills(- zb, 1) * right ^ 2
            }
        }
        if (TRUNC){
            if (LEFT){
                h_spec_beta_beta <- - Io * mills(- za, 1) 
                h_spec_beta_sig <-    Io * left * mills(- za, 1)
                h_spec_sig_sig <-   - Io * mills(- za, 1) * left ^ 2
            }
            if (RIGHT){
                h_spec_beta_beta <- - Io * mills(zb, 1)
                h_spec_beta_sig <-    Io * mills(zb, 1) * right
                h_spec_sig_sig <-   - Io * mills(zb, 1) * right ^ 2
            }
            if (TWO){
                DELTA <- pnorm(zb) - pnorm(za)
                A <- mills(zb) * pnorm(zb) - mills(za) * pnorm(za)
                B <- mills(zb) * pnorm(zb) * right - mills(za) * pnorm(za) * left
                MILLS <- A / DELTA
                A_sig <- - pnorm(zb) * (mills(zb) ^ 2 + mills(zb, 1)) * zb + pnorm(za) * (mills(za) ^ 2 + mills(za, 1)) * za
                B_sig <- - pnorm(zb) * zb * (mills(zb, 1) * zb + mills(zb) ^ 2 * zb + mills(zb)) +
                     pnorm(za) * za * (mills(za, 1) * za + mills(za) ^ 2 * za + mills(za))
                D_sig <- - mills(zb) * pnorm(zb) * zb + mills(za) * pnorm(za) * za
                g_spec_beta <- - Io * MILLS
                g_spec_sig <- - Io * (mills(zb) * pnorm(zb) * zb - mills(za) * pnorm(za) * za) / Dphi
                DD <-  (pnorm(zb) * (mills(zb, 1) + mills(zb) ^ 2) - pnorm(za) * (mills(za, 1) + mills(za) ^ 2)) / Dphi -
                    (mills(zb) * pnorm(zb) - mills(za) * pnorm(za)) ^ 2 / (pnorm(zb) - pnorm(za)) ^ 2
                h_spec_beta_beta <- Io * DD
                h_spec_beta_sig <-   Io * (MILLS - (A_sig - MILLS * D_sig) / DELTA)
                h_spec_sig_sig <- Io * (B / DELTA - (B_sig * DELTA - B * D_sig) / DELTA ^ 2)
            }
        }
        H_bb <- crossprod(wt * (h_com_beta_beta + h_spec_beta_beta) * X, X)
        H_bs <- apply(wt * (h_com_beta_sig + h_spec_beta_sig) * X, 2, sum)
        H_ss <- sum(wt * (h_com_sig_sig + h_spec_sig_sig))
        attr(lnl, "hessian") <- rbind(cbind(H_bb, H_bs),
                                      c(H_bs, H_ss))
      }
    lnl
}
    


