### R code from vignette source 'spareg.Rnw'

###################################################
### code chunk number 1: r setup
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("ggplot2")
library("spareg")
is_paper <- FALSE
tmpdir <- tempdir()


###################################################
### code chunk number 2: condition_vignette
###################################################
if (is_paper) {
  cat("\\Positivetrue")
} else {
  cat("\\Positivefalse")
}


###################################################
### code chunk number 3: spareg.Rnw:558-559 (eval = FALSE)
###################################################
## install.packages("spareg")


###################################################
### code chunk number 4: spareg.Rnw:562-563
###################################################
library("spareg")


###################################################
### code chunk number 5: spareg.Rnw:570-573
###################################################
set.seed(1234)
example_data <- simulate_spareg_data(n = 200, p = 2000, ntest = 100)
str(example_data)


###################################################
### code chunk number 6: spareg.Rnw:789-790
###################################################
obj <- screen_marglik()


###################################################
### code chunk number 7: spareg.Rnw:794-795
###################################################
obj


###################################################
### code chunk number 8: spareg.Rnw:798-799
###################################################
unclass(obj)


###################################################
### code chunk number 9: spareg.Rnw:945-947
###################################################
obj <- rp_gaussian()
obj


###################################################
### code chunk number 10: spareg.Rnw:950-951
###################################################
unclass(obj)


###################################################
### code chunk number 11: first_example
###################################################
set.seed(12)
spar_res <- spar(example_data$x, example_data$y,
  xval = example_data$xtest, yval = example_data$ytest,
  nummods = c(5, 10, 15, 20, 25, 30))
spar_res


###################################################
### code chunk number 12: spareg.Rnw:1030-1033
###################################################
spar_cv <- spar.cv(example_data$x, example_data$y,
  nummods = c(5, 10, 15, 20, 25, 30))
spar_cv


###################################################
### code chunk number 13: spareg.Rnw:1047-1048
###################################################
str(coef(spar_res))


###################################################
### code chunk number 14: spareg.Rnw:1061-1062
###################################################
str(coef(spar_cv))


###################################################
### code chunk number 15: spareg.Rnw:1158-1163 (eval = FALSE)
###################################################
## plot(spar_res)
## plot(spar_res, plot_type = "Val_numAct")
## plot(spar_res,  plot_type = "coefs")
## plot(spar_res, plot_type = "res-vs-fitted", xfit = example_data$xtest,
##   yfit = example_data$ytest)


###################################################
### code chunk number 16: plotmethod
###################################################
library(ggplot2)
p1 <- plot(spar_res)+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
p2 <- plot(spar_res, plot_type = "Val_numAct") +
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
p4 <- plot(spar_res, plot_type = "res-vs-fitted", xfit = example_data$xtest,
  yfit = example_data$ytest)
p3 <- plot(spar_res,  plot_type = "coefs")
p <- ggpubr::ggarrange(p1, p2, p3, p4,
       ncol = 4, nrow = 1)
p


###################################################
### code chunk number 17: parallel_example (eval = FALSE)
###################################################
## set.seed(1234)
## example_data3 <- simulate_spareg_data(n = 1000, p = 2000, ntest = 1000)
## library(doParallel)
## cl <- makeCluster(2, type = "PSOCK")
## registerDoParallel(cl)
## spar_res_par <- spar(example_data3$x, example_data3$y,
##   screencoef = screen_cor(), rp = rp_gaussian(),
##   nummods = 50, parallel = TRUE)
## stopCluster(cl)


###################################################
### code chunk number 18: spareg.Rnw:1231-1238
###################################################
generate_scr_sirs <- function(y, x, object) {
  res_screen <- do.call(function(...)
    VariableScreening::screenIID(x, y, ...),
    object$control)
  coefs <- res_screen$measurement
  coefs
}


###################################################
### code chunk number 19: sirs_chunk
###################################################
screen_sirs <- constructor_screencoef(
  "screen_sirs",
  generate_fun = generate_scr_sirs)


###################################################
### code chunk number 20: screen_sirs
###################################################
set.seed(123)
spar_example <- spar(example_data$x, example_data$y,
  screencoef = screen_sirs(type = "fixed",
    control = list(method = "SIRS")),
  rp = rp_sparse(psi = 1/sqrt(ncol(example_data$x))), measure = "mse")
spar_example


###################################################
### code chunk number 21: spareg.Rnw:1283-1288
###################################################
simulate_haar <- function(m, p) {
  R0 <- matrix(1/sqrt(p) * rnorm(p * m), nrow = p, ncol = m)
  RM <- qr.Q(qr(R0), complete = FALSE)
  t(RM)
}


###################################################
### code chunk number 22: spareg.Rnw:1320-1344
###################################################
update_rp_cannings <- function(...) {
  args <- rlang::list2(...)
  if (is.null(args$rp$control$family)) {
    family_string <- paste0(args$family$family,
                            "(", args$family$link, ")")
    args$rp$control$family_string <- family_string
    family <- args$family
    fit_family <- switch(family$family,
      "gaussian" = if (family$link == "identity") "gaussian" else family,
      "binomial" = if (family$link == "logit") "binomial" else family,
      "poisson"  = if (family$link == "log") "poisson" else family,
      family)
    args$rp$control$fit_family <- fit_family
  }
  if (is.null(args$rp$control$alpha)) args$rp$control$alpha <- 1
  if (is.null(args$rp$control$B2))  args$rp$control$B2 <- 50
  if (is.null(args$rp$control$xi)) {
    args$rp$control$xi <- 0.25
  } else {
    stopifnot("xi must be between 0 and 1."=
      args$rp$control$xi >= 0 &  args$rp$control$xi<= 1)
  }
  args$rp
}


###################################################
### code chunk number 23: generate_cannings
###################################################
generate_cannings <- function(rp, m, included_vector, x, y) {
  xs <- x[, included_vector]
  n <- nrow(x);  p <- ncol(xs)
  B2 <- rp$control$B2;  xi <- rp$control$xi
  id_test <- sample(n, size = n * xi)
  xtrain <- xs[-id_test, ];  xtest <- xs[id_test,]
  ytrain <- y[-id_test];  ytest <- y[id_test]
  control_glmnet <-
    rp$control[names(rp$control) %in% names(formals(glmnet::glmnet))]
  best_val <- 1e5
  family <- eval(parse(text = rp$control$family_string))
  for (b in seq_len(B2)) {
    RM <- simulate_haar(m, p)
    xrp <- tcrossprod(xtrain, RM)
    mod <- do.call(function(...)
      glmnet::glmnet(x = xrp, y = ytrain,
        family = rp$control$fit_family, ...),
      control_glmnet)
    coefs <- coef(mod, s = min(mod$lambda))
    eta_test <- (cbind(1, tcrossprod(xtest, RM)) %*% coefs)
    pred <- family$linkinv(as.vector(eta_test))
    out_perf <-  ifelse(family$family == "binomial",
      mean(((pred > 0.5) + 0) != ytest),
      mean((pred - ytest)^2))
    if (out_perf < best_val) {
      best_val <- out_perf; best_RM <- RM
    }
    rm(RM)
  }
  return(best_RM)
}


###################################################
### code chunk number 24: define_rp_cannings
###################################################
rp_cannings <- constructor_randomprojection(
  "rp_cannings",
  generate_fun = generate_cannings,
  update_fun = update_rp_cannings
)


###################################################
### code chunk number 25: data_cannings
###################################################
set.seed(1234)
example_data2 <- simulate_spareg_data(n = 100, p = 1000, ntest = 100)
ystar <- (example_data2$y > 0) + 0
ystarval <- (example_data2$ytest > 0) + 0


###################################################
### code chunk number 26: run_cannings
###################################################
file <- "cannings_example.rda"
if (file.exists(file)) {
  load(file)
} else {
  set.seed(12345)
  Sys.time()
  spar_example_1 <- spar(
    x = example_data2$x, y = ystar,
    xval = example_data2$xtest, yval = ystarval,
    family = binomial(),
    nus = 0,
    nummods = 50,
    rp = rp_cannings(control = list(lambda.min.ratio = 0.01)),
    measure = "class"
  )
  Sys.time()
  set.seed(12345)
  spar_example_2 <- spar(
    x = example_data2$x, y = ystar,
    xval = example_data2$xtest, yval = ystarval,
    family = binomial(),
    rp = rp_cw(data = TRUE),
    nus = 0, nummods = 50,
    measure = "class"
  )
  save(spar_example_1, spar_example_2,
       file = file)
}


###################################################
### code chunk number 27: spareg.Rnw:1436-1443 (eval = FALSE)
###################################################
## set.seed(12345)
## spar_example_1 <- spar(x = example_data2$x, y = ystar,
##   family = binomial(),
##   rp = rp_cannings(control = list(lambda.min.ratio = 0.01)),
##   nus = 0, nummods = 50,
##   xval = example_data2$xtest, yval = ystarval,
##   measure = "class")


###################################################
### code chunk number 28: spareg.Rnw:1446-1451 (eval = FALSE)
###################################################
## set.seed(12345)
## spar_example_2 <- spar(x = example_data2$x, y = ystar,
##   family = binomial(), rp = rp_cw(data = TRUE),
##   nus = 0, nummods = 50, xval = example_data2$xtest, yval = ystarval,
##   measure = "class")


###################################################
### code chunk number 29: compare
###################################################
spar_example_1$val_res[which.min(spar_example_1$val_res$Meas), ]
spar_example_2$val_res[which.min(spar_example_2$val_res$Meas), ]


###################################################
### code chunk number 30: spareg.Rnw:1469-1485
###################################################
model_glmrob <- function(y, z, object) {
  requireNamespace("robustbase")
  fam <- object$control$family
  if (fam$family == "gaussian" & fam$link == "identity") {
    glmrob_res <- do.call(function(...)
      robustbase::lmrob(y ~ as.matrix(z), ...),
      object$control)
  } else {
    glmrob_res <- do.call(function(...)
      robustbase::glmrob(y ~ as.matrix(z), ...),
      object$control)
  }
  intercept <- coef(glmrob_res)[1]
  gammas <- coef(glmrob_res)[-1]
  list(gammas = gammas, intercept = intercept)
}


###################################################
### code chunk number 31: spareg.Rnw:1491-1493
###################################################
spar_glmrob <- constructor_sparmodel(name = "glmrob",
  model_fun = model_glmrob)


###################################################
### code chunk number 32: contaminate_x_y
###################################################
perc_cont <- 0.25
x <- example_data2$x;
set.seed(123)
np <- ncol(x) * nrow(x)
id_outliers_x <- sample(seq_len(np), perc_cont * np)
x[id_outliers_x] <- x[id_outliers_x] + 50


###################################################
### code chunk number 33: compare_glmrob_glm
###################################################
set.seed(1234)
spar_rob_res <- spar(x, ystar, family = binomial(),
  model = spar_glmrob(), rp = rp_gaussian(msup = 25),
  measure = "mae")
set.seed(1234)
spar_res <- spar(x, ystar, family = binomial(),
  model = spar_glm(), rp = rp_gaussian(msup = 25),
  measure = "mae")
spar_rob_res
spar_res


###################################################
### code chunk number 34: urls
###################################################
if (is_paper) {
    url1 <- "https://web.archive.org/web/20150922051706/"
    url2 <- "http://isomap.stanford.edu/face_data.mat.Z"
    ## On Ubuntu or MacOS use code below
    ## On Windows, download locally and unzip
    if (!file.exists(file.path(tmpdir, "face_data.mat"))) {
      # system('uncompress face_data.mat.Z')
      library("R.matlab")
      download.file(paste0(url1, url2),
                    file.path(tmpdir, "face_data.mat.Z"))
      system(sprintf('uncompress %s', file.path(tmpdir, "face_data.mat.Z")))
    }
}


###################################################
### code chunk number 35: unzip_file
###################################################
if (is_paper) {
  library("R.matlab")
  facedata <- readMat(file.path(tmpdir,"face_data.mat"))
  x <- t(facedata$images)
  y <- facedata$poses[1,]
}


###################################################
### code chunk number 36: read_matlab (eval = FALSE)
###################################################
## library("R.matlab")
## facedata <- readMat("face_data.mat")
## x <- t(facedata$images)
## y <- facedata$poses[1,]


###################################################
### code chunk number 37: std_x
###################################################
if (is_paper) x[, apply(x, 2, sd) < 0.01] <- 0


###################################################
### code chunk number 38: std_x (eval = FALSE)
###################################################
## x[, apply(x, 2, sd) < 0.01] <- 0


###################################################
### code chunk number 39: std_x
###################################################
if (is_paper) {
  set.seed(123)
  ntot <- length(y); ntest <- ntot * 0.25
  testind <- sample(ntot, ntest, replace = FALSE)
  xtrain <- as.matrix(x[-testind, ]); ytrain <- y[-testind]
  xtest <- as.matrix(x[testind, ]); ytest <- y[testind]
}


###################################################
### code chunk number 40: spareg.Rnw:1602-1607 (eval = FALSE)
###################################################
## set.seed(123)
## ntot <- length(y); ntest <- ntot * 0.25
## testind <- sample(ntot, ntest, replace = FALSE)
## xtrain <- as.matrix(x[-testind, ]); ytrain <- y[-testind]
## xtest <- as.matrix(x[testind, ]); ytest <- y[testind]


###################################################
### code chunk number 41: spareg.Rnw:1622-1643
###################################################
if (is_paper) {
  file <- file.path(tmpdir, "faces_res.rda")
  if (!file.exists(file)) {
    set.seed(123)
    control_glmnet <- list(lambda.min.ratio = 0.001)
    library("spareg")
    Sys.time()
    spar_faces <- spar.cv(
      xtrain, ytrain,
      model = spar_glm(),
      screencoef = screen_glmnet(control = control_glmnet,
                                 reuse_in_rp = TRUE),
      rp = rp_cw(data = TRUE, ),
      nummods = c(5, 10, 20, 50),
      measure = "mse")
    Sys.time()
    save(spar_faces, file = file)
  } else {
    load(file)
  }
}


###################################################
### code chunk number 42: spareg.Rnw:1645-1653 (eval = FALSE)
###################################################
## library("spareg")
## set.seed(123)
## control_glmnet <- list(lambda.min.ratio = 0.001)
## spar_faces <- spar.cv(xtrain, ytrain,
##   model = spar_glm(),
##   screencoef = screen_glmnet(control = control_glmnet, reuse_in_rp = TRUE),
##   rp = rp_cw(data = TRUE),
##   nummods = c(5, 10, 20, 50), measure = "mse")


###################################################
### code chunk number 43: spareg.Rnw:1657-1658 (eval = FALSE)
###################################################
## spar_faces


###################################################
### code chunk number 44: spareg.Rnw:1660-1661
###################################################
if (is_paper) spar_faces


###################################################
### code chunk number 45: plotfacesmeasure (eval = FALSE)
###################################################
## plot(spar_faces)


###################################################
### code chunk number 46: plotfacesmeasurenummod5 (eval = FALSE)
###################################################
## plot(spar_faces, nummod = 5)


###################################################
### code chunk number 47: plotfacesmeasurereal_best
###################################################
if (is_paper) {
  p <- plot(spar_faces, digits = 3L) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     hjust = 1))
  p
}


###################################################
### code chunk number 48: plotfacesmeasurereal_1se
###################################################
if (is_paper) {
  p2 <- plot(spar_faces, nummod=5, digits = 3L) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     hjust = 1))
  p2
}


###################################################
### code chunk number 49: facecoef (eval = FALSE)
###################################################
## face_coef <- coef(spar_faces, opt_par = "best")
## str(face_coef)


###################################################
### code chunk number 50: facecoef
###################################################
if (is_paper) {
  face_coef <- coef(spar_faces, opt_par = "best")
  str(face_coef)
}


###################################################
### code chunk number 51: facecoef1se (eval = FALSE)
###################################################
## face_coef_1se <- coef(spar_faces, opt_par = "1se")
## str(face_coef_1se)


###################################################
### code chunk number 52: facecoef1se
###################################################
if (is_paper) {
  face_coef_1se <- coef(spar_faces, opt_par = "1se")
  str(face_coef_1se)
}


###################################################
### code chunk number 53: predictfaces
###################################################
if (is_paper) {
  ynew <- predict(spar_faces, xnew = xtest, coef = face_coef_1se)
}


###################################################
### code chunk number 54: predictfaces (eval = FALSE)
###################################################
## ynew <- predict(spar_faces, xnew = xtest, coef = face_coef_1se)


###################################################
### code chunk number 55: tarp_comparison_true (eval = FALSE)
###################################################
## tarp_faces <-  spar.cv(xtrain, ytrain,
##   model = spar_glm(),
##   screencoef = screen_cor(),
##   rp = rp_sparse(psi = 1/3),
##   nummods = c(5, 10, 20, 50), measure = "mse")
## ynew_tarp <- predict(tarp_faces, xnew = xtest,
##   coef = coef(tarp_faces, opt_par = "1se"))
## c("SPAR" = mean((ytest-ynew)^2),
##   "TARP" = mean((ytest-ynew_tarp)^2))


###################################################
### code chunk number 56: tarp_comparison_true
###################################################
if(is_paper) {
  tarp_faces <-  spar.cv(xtrain, ytrain,
                         model = spar_glm(),
                         screencoef = screen_cor(), rp = rp_sparse(psi = 1/3),
                         nummods = c(5, 10, 20, 50), measure = "mse")
  ynew_tarp <- predict(tarp_faces, xnew = xtest,
                       coef = coef(tarp_faces, opt_par = "1se"))
  c("SPAR" = mean((ytest-ynew)^2),  "TARP" = mean((ytest-ynew_tarp)^2))
}


###################################################
### code chunk number 57: plotfacesobs179
###################################################
if (is_paper) {
  i <- 179
  p <- ggplot(data.frame(X = rep(1:64,each=64),
                         Y = rep(64:1,64),
                         Z = facedata$images[,i]),
              aes(X, Y, fill = Z)) +
    geom_tile() +
    theme_void() +
    ggtitle(paste0("y = ", round(facedata$poses[1, i],1))) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  p
}


###################################################
### code chunk number 58: plotobsvspred179
###################################################
if (is_paper) {
  id <- 3
  p4 <- ggplot(data.frame(X = rep(1:64, each = 64),
                             Y = rep(64:1, 64),
                             effect = xtest[id,] * face_coef_1se$beta),
                  aes(X, Y, fill = effect)) +
    geom_tile() +
    theme_void() +
    scale_fill_gradient2() +
    ggtitle(bquote(hat(y) == .(round(ynew[id]))), ) +
    theme(plot.title = element_text(hjust = 0.5))
  p4
}


###################################################
### code chunk number 59: readdarwin (eval = FALSE)
###################################################
## tmpdir <- tempdir()
## download.file("https://archive.ics.uci.edu/static/public/732/darwin.zip",
##               file.path(tmpdir, "darwin.zip"))


###################################################
### code chunk number 60: readdarwinreal
###################################################
if (is_paper) {
  # if (file.exists(file.path(tmpdir, "data.zip"))) {
  #  unzip(file.path(tmpdir, "data.zip")exdir = )
  #} else {
  if (!file.exists(file.path(tmpdir, "darwin.zip"))) {
      download.file("https://archive.ics.uci.edu/static/public/732/darwin.zip",
                    file.path(tmpdir, "darwin.zip"))
      unzip(zipfile = file.path(tmpdir, "darwin.zip"),
            files = "data.csv",
            exdir = tmpdir)
  }
  darwin_tmp <- read.csv(file.path(tmpdir, "data.csv"),
                       stringsAsFactors = TRUE)
}


###################################################
### code chunk number 61: readdarwinzip (eval = FALSE)
###################################################
## unzip(zipfile = file.path(tmpdir, "darwin.zip"),
##   files = "data.csv",
##   exdir = file.path(tmpdir))
## darwin_tmp <- read.csv(file.path(tmpdir, "data.csv"),
##   stringsAsFactors = TRUE)


###################################################
### code chunk number 62: darwinimpute (eval = FALSE)
###################################################
## darwin_orig <- list(
##   x = darwin_tmp[, !(colnames(darwin_tmp) %in% c("ID", "class"))],
##   y = as.numeric(darwin_tmp$class) - 1)
## tmp <- cellWise::DDC(
##   as.matrix(darwin_orig$x),
##   list(returnBigXimp = TRUE, tolProb = 0.999, silent = TRUE))
## darwin <- list(x = tmp$Ximp, y = darwin_orig$y)


###################################################
### code chunk number 63: spareg.Rnw:1861-1870
###################################################
if (is_paper) {
  darwin_orig <- list(
    x = darwin_tmp[, !(colnames(darwin_tmp) %in% c("ID", "class"))],
    y = as.numeric(darwin_tmp$class) - 1)
  tmp <- cellWise::DDC(
    as.matrix(darwin_orig$x),
    list(returnBigXimp = TRUE, tolProb = 0.999, silent = TRUE))
  darwin <- list(x = tmp$Ximp, y = darwin_orig$y)
}


###################################################
### code chunk number 64: darwinstr (eval = FALSE)
###################################################
## str(darwin)


###################################################
### code chunk number 65: spareg.Rnw:1876-1877
###################################################
if (is_paper) str(darwin)


###################################################
### code chunk number 66: spardarwin (eval = FALSE)
###################################################
## set.seed(1234)
## spar_darwin <- spar.cv(darwin$x, darwin$y,
##   family = binomial(logit), screencoef = screen_glmnet(reuse_in_rp = TRUE),
##   nummods = c(10, 20, 30, 50), measure = "1-auc")
## spar_darwin


###################################################
### code chunk number 67: spardarwinreal
###################################################
if (is_paper) {
  file <- file.path(tmpdir, "darwin_res.rda")
  if (!file.exists(file)) {
    set.seed(1234)
    spar_darwin <- spareg::spar.cv(darwin$x, darwin$y,
                           family = binomial(logit),
                           screencoef = screen_glmnet(reuse_in_rp = TRUE),
                           nummods = c(10, 20, 30, 50),
                           measure = "1-auc")
    save(spar_darwin, file = file)
  } else {
    load(file = file)
  }
  spar_darwin
}


###################################################
### code chunk number 68: plotdarwinact (eval = FALSE)
###################################################
## plot(spar_darwin, plot_type = "Val_numAct")


###################################################
### code chunk number 69: plotdarwinactnumAct
###################################################
if (is_paper) {
  p <- plot(spar_darwin, plot_type = "Val_numAct", digits=3L)
  p
}


###################################################
### code chunk number 70: plotdarwincoef
###################################################
if (is_paper) {
  ntasks <- 25
  nfeat <- 18
  reorder_ind <- c(outer(
    (seq_len(ntasks) - 1) * nfeat,
    seq_len(nfeat), "+"))
  feat_names <- sapply(colnames(darwin$x)[seq_len(nfeat)],
                       function(name) substr(name, 1, nchar(name) - 1))
  p <- plot(spar_darwin,"coefs",coef_order = reorder_ind) +
    geom_vline(xintercept = 0.5 + seq_len(ntasks - 1) * ntasks,
               alpha = 0.2, linetype = 2) +
    annotate("text",x = (seq_len(nfeat) - 1) * ntasks + 12,
             y = 40,
             label=feat_names, angle = 90,
             size = 3)
  p
}


###################################################
### code chunk number 71: session_info
###################################################
sessionInfo()


