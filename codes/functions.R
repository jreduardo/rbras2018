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
        fixed <- list(X = X, Z = Z, y = y, sumto = sumto,
                      formula_beta = ffx, formula_gama = ffz)
        fit <- suppressWarnings(
            bbmle::mle2(lldcmp,
                        start = start,
                        data = fixed,
                        method = "BFGS",
                        vecpar = TRUE)
        )
    }
    if (strategy == "fixed") {
        np <- ncol(Z)
        bbmle::parnames(lldcmp) <- names(start[1:np])
        fixed <- list(beta = start[-(1:np)], X = X,
                      Z = Z, y = y, sumto = sumto,
                      formula_beta = ffx, formula_gama = ffz)
        gama_fit <- suppressWarnings(
            bbmle::mle2(lldcmp,
                        start = start[1:np],
                        data = fixed,
                        method = "BFGS",
                        vecpar = TRUE)
        )
        start <- c(gama_fit@coef, start[-(1:np)])
        bbmle::parnames(lldcmp) <- names(start)
        fit <- suppressWarnings(
            bbmle::mle2(lldcmp,
                        start = start,
                        data = fixed[-1],
                        method = "BFGS",
                        vecpar = TRUE,
                        control = list(maxit = 0L))
        )
        fit@details$convergence <- 0L
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
    expec_y2 <- sum(yrange^2 * dcmp(yrange, mu, phi, sumto = sumto))
    variance <- expec_y2 - mu^2
    return(variance)
}
compute_variance <- Vectorize(
    compute_variance,
    vectorize.args = c("mu", "phi"))


# To organize the output of profile log-likelihood
myprofile <- function(model, which = 1:length(coef(model))) {
    sel <- c("param", "z", "focal")
    prof <- profile(model, which = which)
    as.data.frame(prof)[, sel]
}

# Get estimates and standard error for CMP outputs
get_coef <- function(object,
                     type = c("both", "gama", "beta"),
                     digits = getOption("digits"),
                     compact = TRUE,
                     add_letter = TRUE,
                     math_style = TRUE,
                     keep_rownames = FALSE) {
    type <- match.arg(type)
    x <- bbmle::summary(object)@coef
    x <- switch(
        type,
        "both" = x,
        "gama" = x[grepl("gama", rownames(x)), ,drop = FALSE],
        "beta" = x[grepl("beta", rownames(x)), ,drop = FALSE]
    )
    est <- x[, "Estimate"]
    std <- x[, "Std. Error"]
    ind <- abs(est / std) > qnorm(0.975)
    if (compact) {
        fest <- format(round(est, digits))
        fstd <- format(round(std, digits))
        if (math_style) fest <- gsub("-", "$-$", fest)
        res <- sprintf("%s (%s)", fest, fstd)
        if (add_letter) {
            tex <- if (!math_style) {
                       c(" a", "  ")
                   } else c("$^\\text{a}$",
                            "\\textcolor{white}{$^\\text{a}$}")
            let <- ifelse(ind, tex[1], tex[2])
            res <- paste0(res, let)
        }
        out <- data.frame("est (std)" = res, check.names = FALSE)
    } else {
        out <- data.frame(est = est, std = std)
        if (add_letter) out$sig <- ifelse(ind, "a", "")
    }
    if (keep_rownames) rownames(out) <- rownames(x)
    return(out)
}

# Complete estimates tables for multiples nested models
complete_coef <- function(object) {
    nmax <- max(vapply(object, nrow, integer(1)))
    lapply(object, function(x) {
        if (nrow(x) == nmax) return(x)
        x[(nrow(x) + 1):nmax, ] <- NA
        return(x)
    })
}

# Analysis of deviance table (TRV for nested models)
get_anova <- function(object, ..., print = TRUE) {
    if (is.list(object)) {
        mlist <- object
    } else {
        mlist <- c(object, list(...))
    }
    n   <- vapply(mlist, function(x) nrow(x@data$X), 0L)
    nps <- vapply(mlist, function(x) attr(logLik(x), "df"), 0)
    aic <- vapply(mlist, function(x) AIC(x), 0)
    dev <- vapply(mlist, function(x) -2*logLik(x), 0)
    cst <- -diff(dev)
    pvs <- pchisq(cst, df = diff(nps), lower.tail = FALSE)
    tab <- data.frame("gl" = n - nps,
                      "dev" = dev,
                      "AIC" = aic,
                      "chisq" = c(NA, cst),
                      "Pr(>Chisq)" = c(NA, pvs),
                      check.names = FALSE)
    rownames(tab) <- names(mlist)
    if (print) printCoefmat(tab, na.print = "", cs.ind = 1)
    invisible(tab)
}

# Prediction mean and dispersion for double COM-Poisson models
predict_cmp <- function(model,
                        newdata,
                        predict = c("mean", "dispersion"),
                        type = c("response", "link"),
                        interval = c("confidence", "none"),
                        level = 0.95,
                        augment_data = TRUE) {
    type <- match.arg(type)
    predict <- match.arg(predict)
    interval <- match.arg(interval)
    #-------------------------------------------
    Vcov <- vcov(model)
    indb <- grep("beta", rownames(Vcov))
    indg <- grep("gama", rownames(Vcov))
    #-------------------------------------------
    Vbeta <- Vcov[indb, indb, drop = FALSE]
    Vgama <- Vcov[indg, indg, drop = FALSE]
    Vbega <- Vcov[indb, indg]
    #-------------------------------------------
    if (predict == "mean") {
        Vcond <- Vbeta - tcrossprod(Vbega %*% solve(Vgama), Vbega)
        formula <- model@data$formula_beta
        formula[[2]] <- NULL
        coefs <- model@coef[indb]
    }
    if (predict == "dispersion") {
        Vcond <- Vgama - crossprod(Vbega, solve(Vbeta) %*% Vbega)
        formula <- model@data$formula_gama
        coefs <- model@coef[indg]
    }
    #-------------------------------------------
    Mat <- model.matrix(formula, newdata)
    est <- Mat %*% coefs
    if (interval == "none") {
        out <- data.frame("fit" = est)
    }
    if (interval == "confidence") {
        qn <- -qnorm((1 - level[1])/2)
        std <- sqrt(diag(tcrossprod(Mat %*% Vcond, Mat)))
        out <- data.frame("fit" = est,
                          "lwr" = est - qn * std,
                          "upr" = est + qn * std)
    }
    if (type == "response") out <- data.frame(apply(out, 2L, exp))
    if (augment_data) out <- cbind(newdata, out)
    return(out)
}
