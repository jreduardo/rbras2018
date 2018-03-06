# =====================================================================
# Fitting double COM-Poisson models
# Joint modelling mean and variance for count data
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2018-03-05
# =====================================================================

# Negative of log-likelihood function
lldcmp <- function(params, beta = NULL,
                   X, Z, y, sumto = 500L) {
    n <- length(y)
    npar1 <- ncol(Z)
    # Obtem os parametros separados
    if (is.null(beta)) beta <- params[-(1:npar1)]
    gama <- params[1:npar1]
    nu <- exp(Z %*% gama)
    xb <- X %*% beta
    aux <- exp(xb) + (nu - 1) / (2 * nu)
    aux[aux < 0] <- 1e-10
    tl <- log(aux)
    # Obtem as constantes
    j <- 1:sumto
    zi <- vapply(1:n, FUN = function(i) {
        1 + sum(exp(nu[i] * j * tl[i] - nu[i] * lfactorial(j)))
    }, FUN.VALUE = numeric(1))
    # Calcula o negativo do log da função de verossimilhança
    ll <- sum(y * nu * tl) - sum(nu * lfactorial(y)) - sum(log(zi))
    return(-ll)
}

# Framework for fitting double COM-Poisson models
cmp <- function(formula, data, start = NULL, sumto = 500L,
                strategy = c("joint", "fixed")) {
    strategy <- match.arg(strategy)
    #-------------------------------------------
    # Separando as formulas
    if(length(formula[[3]]) > 1L &&
       formula[[3]][[1]] == as.name("|")) {
        ff <- formula
        formula[[3]][1] <- call("+")
        ffx <- . ~ .
        ffz <- ~ .
        ffx[[3]] <- ff[[3]][[2]]
        ffx[[2]] <- ff[[2]]
        ffz[[3]] <- ff[[3]][[3]]
        ffz[[2]] <- NULL
    } else {
        ffz <- ffc <- ff <- formula
        ffz[[2]] <- NULL
    }
    #-------------------------------------------
    # Obtendo as matrizes
    frame_beta <- model.frame(ffx, data)
    frame_gama <- model.frame(ffz, data)
    terms_beta <- attr(frame_beta, "terms")
    terms_gama <- attr(frame_gama, "terms")
    X <- model.matrix(terms_beta, frame_beta)
    Z <- model.matrix(terms_gama, frame_gama)
    y <- model.response(frame_beta)
    #-------------------------------------------
    # Parametros iniciais
    if (is.null(start)) {
        m0 <- glm.fit(x = X, y = y, family = poisson())
        start <- c(rep(0, ncol(Z)), "beta" = m0$coefficients)
        names(start)[1:ncol(Z)] <- paste0("gama.", colnames(Z))
    }
    #-------------------------------------------
    # Ajuste do modelo
    if (strategy == "joint") {
        bbmle::parnames(lldcmp) <- names(start)
        fixed <- list(X = X, Z = Z, y = y, sumto = sumto)
        fit <- bbmle::mle2(lldcmp,
                           start = start,
                           data = fixed,
                           method = "BFGS",
                           vecpar = TRUE)
    }
    if (strategy == "fixed") {
        np <- ncol(Z)
        bbmle::parnames(lldcmp) <- names(start[1:np])
        fixed <- list(beta = start[-(1:np)], X = X,
                      Z = Z, y = y, sumto = sumto)
        fit <- bbmle::mle2(lldcmp,
                           start = start[1:np],
                           data = fixed,
                           method = "BFGS",
                           vecpar = TRUE)
    }
    return(fit)
}

# Probability mass function for COM-Poisson
dcmp <- function(y, mu, phi, sumto = 500L) {
    vapply(y, function(yi) {
        exp(-lldcmp(params = c(phi, log(mu)),
                    X = cbind(1L),
                    Z = cbind(1L),
                    y = yi))
    }, numeric(1))
}

# Compute variance for a COM-Poisson variable
compute_variance <- function(mu, phi, sumto = 500L, tol = 1e-6) {
    # (poor) approximate variance
    nu <- exp(phi)
    ap_stderr <- sqrt(mu / nu)
    ymax <- ceiling(mu + 5 * ap_stderr)
    pmax <- dcmp(ymax, mu, phi, sumto = sumto)
    # Verifica se prob(ymax) é pequena o suficiente
    while (pmax > tol) {
        ymax <- ymax + 1L
        pmax <- dcmp(ymax, mu, phi, sumto = sumto)
    }
    yrange <- 1:ymax
    expec_y2 <- sum(yrange^2 * dcmp(yrange, mu, phi))
    variance <- expec_y2 - mu^2
    return(variance)
}
