#-----------------------------------------------------------------------
# Setup

#-----------------------------------------------------------------------
# Settings for outputs and load packages
## ---- setup

library(knitr)
library(xtable)
options(digits = 3, OutDec = ".",
        xtable.caption.placement = "top",
        xtable.booktabs = TRUE,
        xtable.sanitize.text.function = identity,
        xtable.sanitize.colnames.function = identity)

opts_chunk$set(
    warning = FALSE,
    message = FALSE,
    echo = FALSE,
    results = "hide",
    ## dev = "tikz",
    fig.width = 7,
    fig.height = 5,
    out.width = "0.99\\textwidth",
    fig.align = "center",
    fig.pos = "h",
    dev.args = list(family = "Palatino"))

library(dplyr)
library(tidyr)
library(purrr)
library(lattice)
library(latticeExtra)
library(gridExtra)
source("../codes/functions.R")
source("../codes/lattice-panels.R")

#-----------------------------------------------------------------------
# Read and descriptive of data

## ---- load-data
data(Paula, package = "labestData")
nitrofen <- PaulaEx4.6.20
nitrofen <- transform(nitrofen, dose = dose / 100)

orthopoly <- poly(nitrofen$dose, degree = 3)
colnames(orthopoly) <- paste0("dose", 1:3)
nitrofen <- cbind(nitrofen, orthopoly)

# Exploratory analysis
## ---- desc-nitrofen
xy1 <- xyplot(novos ~ dose,
              data = nitrofen,
              xlab = "Nitrofen concentration level",
              ylab = "Number of live offspring",
              type = c("p", "g", "smooth"),
              sub = "(a)",
              spread = 0.1,
              panel = panel.beeswarm)
mv <- nitrofen %>%
    group_by(dose) %>%
    summarise(mu = mean(novos), va = var(novos))
xlim <- ylim <- extendrange(c(mv$mu, mv$va), f = 0.05)

xy2 <- xyplot(va ~ mu,
              data = mv,
              type = c("p", "g"),
              xlim = xlim,
              ylim = ylim,
              xlab = expression("Sample mean"~(bar(y))),
              ylab = expression("Sample variance"~(s^2)),
              sub = "(b)",
              panel = function(x, y, ...) {
                  panel.xyplot(x, y, ...)
                  panel.text(mv$mu, mv$va,
                             format(round(mv$dose, 3)),
                             pos = 4, cex = 0.9)
                  panel.abline(a = 0, b = 1, lty = 2)
              })

print(xy1, split = c(1, 1, 2, 1), more = TRUE)
print(xy2, split = c(2, 1, 2, 1), more = FALSE)

#-------------------------------------------
# Fit models

## ---- fit-models
mnames <- c("Constant", "Linear", "Quadratic", "Cubic")

# CMP Models
library(bbmle)

# Joint estimate both beta and gama parameters
m1CP <- cmp(novos ~ dose1 + dose2 + dose3 | 1,
            data = nitrofen)
m2CP <- cmp(novos ~ dose1 + dose2 + dose3 | dose1,
            data = nitrofen)
m3CP <- cmp(novos ~ dose1 + dose2 + dose3 | dose1 + dose2,
            data = nitrofen)
m4CP <- cmp(novos ~ dose1 + dose2 + dose3 | dose1 + dose2 + dose3,
            data = nitrofen)
modelsCP <- list(m1CP, m2CP, m3CP, m4CP)
names(modelsCP) <- mnames

# Double Poisson models
library(gamlss)
fmu <- novos ~ dose1 + dose2 + dose3
m1DP <- gamlss(fmu, ~ 1,
               family = DPO(),
               data = nitrofen)
m2DP <- gamlss(fmu, ~ dose1,
               family = DPO(),
               data = nitrofen)
m3DP <- gamlss(fmu, ~ dose1 + dose2,
               family = DPO(),
               data = nitrofen)
m4DP <- gamlss(fmu, ~ dose1 + dose2 + dose3,
               family = DPO(),
               data = nitrofen)
modelsDP <- list(m1DP, m2DP, m3DP, m4DP)
names(modelsDP) <- mnames

#-------------------------------------------
# Show results

## ---- tab-coef
cname <- c("Parameter", mnames)
nbeta <- sprintf("$\\beta_%s$", 0:3)
ngama <- sprintf("$\\gamma_%s$", 0:3)

# Estimates for COM-Poisson
tabCP_coef <- do.call(
    rbind,
    lapply(c("beta" = 1, "gama" = 2), function(i) {
        lapply(modelsCP, get_coef,
               type = c("beta", "gama")[i], digits = 3) %>%
            complete_coef() %>%
            do.call(cbind, .) %>%
            cbind("par" = list(nbeta, ngama)[[i]], .)
    }))

# Show the latex code for table
colnames(tabCP_coef) <- sprintf("\\multicolumn{1}{c}{%s}", cname)
print.xtable(xtable(tabCP_coef),
             hline.after = 0,
             only.contents = TRUE,
             include.rownames = FALSE,
             NA.string = "\\multicolumn{1}{c}{$-$}",
             add.to.row = list(
                 list(4),
                 "\\specialrule{0px}{0.2em}{0.3em}\n"))

## ---- tab-anova
tab_anova <- get_anova(modelsCP, print = FALSE)
rownames(tab_anova) <- mnames
colnames(tab_anova) <- c("D.f", "Deviance", "AIC",
                         "$\\rchi^2$", "$\\Pr(>\\rchi^2$)")

print.xtable(xtable(tab_anova, digits = c(0, 0, 3, 3, 4, 4)),
             include.rownames = TRUE,
             only.contents = TRUE,
             NA.string = "\\multicolumn{1}{c}{$-$}")

#-------------------------------------------
# Fitted values

# ---- fit-values
get_preds <- function(model, ...) {
    pred_type <- c("mean", "dispersion")
    names(pred_type) <- pred_type
    purrr::map_dfr(pred_type,
                   predict_cmp,
                   model = model,
                   ...,
                   .id = "prediction")
}

# Fitted mean (\mu_i) and dispersion (\nu_i) values
aux <- data.frame(dose = seq(0, 3, length.out = 100))
orthopoly <- poly(aux$dose, degree = 3)
colnames(orthopoly) <- paste0("dose", 1:3)
aux <- cbind(aux, orthopoly)

pred <-
    modelsCP %>%
    purrr::map_dfr(get_preds, newdata = aux, .id = "model") %>%
    mutate(model = forcats::fct_inorder(factor(model)),
           dose1 = NULL,
           dose2 = NULL,
           dose3 = NULL)

daaux <- subset(pred, prediction == "dispersion")
xy1 <- xyplot(fit ~ dose | model,
              type = "l",
              layout = c(2, 2),
              scales = "free",
              as.table = TRUE,
              xlab = "Nitrofen concentration level",
              ylab = expression(nu),
              sub = "(a)",
              ly = daaux$lwr,
              uy = daaux$upr,
              cty = "bands",
              fill = "gray80",
              alpha = 0.3,
              lty_bands = 3,
              auto.key = TRUE,
              prepanel = prepanel.cbH,
              panel = function(x, y, ..., subscripts) {
                  panel.cbH(x, y, ..., subscripts)
                  panel.segments(min(x), 1, max(x), 1, lty =2)
              },
              data = daaux)

# Compute mean and variance by COM-Poisson distribution
fun <- purrr::possibly(compute_variance, otherwise = NA)
pred2 <-
    pred %>%
    mutate(lwr = NULL, upr = NULL) %>%
    tidyr::spread(key = prediction, value = fit) %>%
    mutate(variance = purrr::map2_dbl(mean, log(dispersion), fun)) %>%
    tidyr::gather(key = moment, value = value, mean:variance)

xy2 <- xyplot(value ~ dose | model,
              groups = moment,
              type = "l",
              layout = c(2, 2),
              as.table = TRUE,
              xlab = "Nitrofen concentration level",
              ylab = " ",
              sub = "(b)",
              auto.key = list(
                  columns = 2,
                  points = FALSE,
                  lines = TRUE,
                  text = c("Mean", "Variance"),
                  cex = 0.9
              ),
              par.settings = list(
                  superpose.line = list(lwd = 2)
              ),
              subset = value > 0L,
              data = pred2)

print(xy1, split = c(1, 1, 2, 1), more = TRUE)
print(xy2, split = c(2, 1, 2, 1), more = FALSE)

# Fit models by different strategies
# ---- fit-strategies
strats <- set_names(c("fixed", "joint"))
formulas <- list(
    novos ~ dose1 + dose2 + dose3 | 1,
    novos ~ dose1 + dose2 + dose3 | dose1,
    novos ~ dose1 + dose2 + dose3 | dose1 + dose2,
    novos ~ dose1 + dose2 + dose3 | dose1 + dose2 + dose3)
names(formulas) <- mnames

# Fit models under two strategies
get_std <- function(object) summary(object)@coef[, 3]
datamodels <- map_dfr(strats, .id = "strategy", function(i)
    tibble(preditor = mnames, formula = formulas)) %>%
    mutate(fit = pmap(., function(..., preditor) 
        cmp(..., data = nitrofen))) %>% 
    mutate(loglik = map_dbl(.$fit, logLik)) %>% 
    mutate(coefs = map(.$fit, coef)) %>% 
    mutate(stds = map(.$fit, get_std))
dacomp <- datamodels %>% 
    dplyr::select(-formula, -fit) %>% 
    unnest()

# Computational times (based on 100 times)
# times <- map_dfr(formulas, .id = "predictor", function(f) {
#    out <- microbenchmark::microbenchmark(
#        "Fixed" = cmp(f, data = nitrofen, strategy = "fixed"),
#        "Joint" = cmp(f, data = nitrofen, strategy = "joint"))
#    as.data.frame(out)
#})
times <- readRDS("times.rds")

# Show maximized likelihoods
# ---- bar-strategies
aux <- bind_cols(filter(dacomp, strategy == "fixed"),
                 filter(dacomp, strategy == "joint")) %>% 
    mutate(diff_co = coefs - coefs1,
           diff_se = stds - stds1,
           diff_ll = loglik - loglik1,
           preditor1 = NULL) %>% 
    dplyr::select(matches("diff|preditor")) %>% 
    gather("stat", "diff", -preditor)
bwplot(~diff | stat,
       scales = "free",
       data = aux)

# ---- compare-strategies
aux <- bind_cols(filter(dacomp, strategy == "fixed"),
                 filter(dacomp, strategy == "joint")) %>% 
    dplyr::select(-matches("preditor|strategy"))
aux <- bind_cols(gather(aux[, 1:3]), gather(aux[, 4:6])) %>% 
    mutate(key = forcats::fct_inorder(factor(key)))
xy1 <- xyplot(value1 ~ value | key, 
              scales = "free",
              grid = TRUE,
              layout = c(3, 1),
              xlab = "Obtained from fixed strategy",
              ylab = "Obtained from joint strategy",
              # sub = "(a)",
              panel = function(...) {
                  panel.xyplot(...)
                  panel.abline(0, 1)
              },
              strip = strip.custom(
                  factor.levels = c("Log-likelihood",
                                    "Estimates",
                                    "Std. Errors")
              ),
              data = aux)
xy2 <- tactile::bwplot2(log10(time) ~ predictor,
                        groups = expr,
                        horizontal = FALSE,
                        grid = TRUE,
                        ylab = expression(log[10]*"(time)"),
                        # sub = "(b)",
                        auto.key = list(
                            # columns = 2
                            space = "right"
                        ),
                        data = times)

print(xy1, position = c(0, 0, 1, 0.5), more = TRUE)
print(xy2, position = c(0.1, 0.5, .9, 1), more = FALSE)

# ---- no-include

xyplot(log10(time) ~ predictor,
       groups = expr,
       #horizontal = FALSE,
       data = times)



str(as.data.frame(microbenchmark::microbenchmark(
    rnorm(1000), rpois(10000,10)
)))





as_tibble(dacomp)

dacomp
dacomp %>% 
    gather(variable, value, loglik:stds) %>%
    unite(temp, strategy, variable) %>%
    spread(temp, value)

dacomp %>% spread(strategy, stds)

# Log-likelihood
barchart(loglik ~ preditor, 
         groups = strategy,
         horizontal = FALSE,
         ylab = "Maximized log-likelihood",
         auto.key = list(
             title = "Strategy",
             cex.title = 1,
             corner = c(0.1, 0.9),
             text = c("Fixed", "Joint")
         ),
         data = dacomp)


%>% 
    mutate(loglik = map_dbl(fit, logLik))

# Joint estimate both beta and gama parameters
mfixed <- lapply(formulas, cmp, data = nitrofen, strategy = "fixed")
mjoint <- lapply(formulas, cmp, data = nitrofen, strategy = "joint")

lapply(formulas, function(f) {
    mfixed <- cmp(data = nitrofen, strategy = "fixed")
    mjoint <- lapply(formulas, cmp, data = nitrofen, strategy = "joint")
    
})




lapply(mfixed, coef)






m1CP <- cmp(,
            data = nitrofen)
m2CP <- cmp(novos ~ dose1 + dose2 + dose3 | dose1,
            data = nitrofen)
m3CP <- cmp(novos ~ dose1 + dose2 + dose3 | dose1 + dose2,
            data = nitrofen)
m4CP <- cmp(novos ~ dose1 + dose2 + dose3 | dose1 + dose2 + dose3,
            data = nitrofen)
modelsCP <- list(m1CP, m2CP, m3CP, m4CP)
names(modelsCP) <- mnames












prof <- profile(m2CP)
plot(prof)

get_coef(m1DP)

str(m1DP)

mat <- summary(m1DP)
prefix <- c(rep("beta", length(m1DP$mu.coefficients)), 
            rep("gamma", length(m1DP$sigma.coefficients)))
rownames(mat) <- sprintf("%s.%s", prefix, rownames(mat))

mat

rownames(mat)

mat[c(1:4), ]

vcov(m1DP)

str(summary(m1DP))
