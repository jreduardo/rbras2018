#-------------------------------------------
# Setting session

library(xtable)
options(xtable.caption.placement = "top",
        xtable.booktabs = TRUE,
        xtable.sanitize.text.function = identity)

library(dplyr)
library(lattice)
library(latticeExtra)

source("functions.R")
source("lattice-panels.R")

#-----------------------------------------------------------------------
# Section: Case study

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
              xlab = "Nível de concentração de nitrofeno",
              ylab = "Número de ovos eclodidos",
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
              xlab = expression("Média amostral"~(bar(y))),
              ylab = expression("Variância amostral"~(s^2)),
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
library(bbmle)
mnames <- c("Escalar", "Linear", "Quadrático", "Cúbico")

# Joint estimate both beta and gama parameters
m1j <- cmp(novos ~ dose1 + dose2 + dose3 | 1,
           data = nitrofen)
m2j <- cmp(novos ~ dose1 + dose2 + dose3 | dose1,
           data = nitrofen)
m3j <- cmp(novos ~ dose1 + dose2 + dose3 | dose1 + dose2,
           data = nitrofen)
m4j <- cmp(novos ~ dose1 + dose2 + dose3 | dose1 + dose2 + dose3,
           data = nitrofen)
models_joint <- list(m1j, m2j, m3j, m4j)
names(models_joint) <- mnames

# Estimate gama with fixed beta according Poisson model
m1f <- cmp(novos ~ dose1 + dose2 + dose3 | 1,
           data = nitrofen, strategy = "fixed")
m2f <- cmp(novos ~ dose1 + dose2 + dose3 | dose1,
           data = nitrofen, strategy = "fixed")
m3f <- cmp(novos ~ dose1 + dose2 + dose3 | dose1 + dose2,
           data = nitrofen, strategy = "fixed")
m4f <- cmp(novos ~ dose1 + dose2 + dose3 | dose1 + dose2 + dose3,
          data = nitrofen)
models_fixed <- list(m1f, m2f, m3f, m4f)
names(models_fixed) <- mnames

#-------------------------------------------
# Show results

## ---- tab-coef
cname <- c("Parâmetro", mnames)
nbeta <- sprintf("$\\beta_%s$", 1:4)
ngama <- sprintf("$\\gamma_%s$", 1:4)

tab_coef <- do.call(
    rbind,
    lapply(c("beta" = 1, "gama" = 2), function(i) {
        lapply(models_joint, get_coef,
               type = c("beta", "gama")[i], digits = 3) %>%
                complete_coef() %>%
            do.call(cbind, .) %>%
            cbind("par" = list(nbeta, ngama)[[i]], .)
    }))

colnames(tab_coef) <- sprintf("\\multicolumn{1}{c}{%s}", cname)
print.xtable(xtable(tab_coef),
             hline.after = 0,
             only.contents = TRUE,
             include.rownames = FALSE,
             NA.string = "\\multicolumn{1}{c}{$-$}",
             add.to.row = list(
                 list(4),
                 "\\specialrule{0px}{0.2em}{0.3em}\n")
             )

## ---- tab-anova
tab_anova <- get_anova(models_joint, print = FALSE)
rownames(tab_anova) <- mnames
colnames(tab_anova) <- c("G.l", "Deviance", "AIC",
                         "$\\rchi^2$", "$\\Pr(>\\rchi^2$)")

print.xtable(xtable(tab_anova, digits = c(0, 0, 3, 3, 4, 4)),
             include.rownames = TRUE,
             only.contents = TRUE,
             NA.string = "\\multicolumn{1}{c}{$-$}")

## ---- noinclude
