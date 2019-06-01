
#check if a value indicates there is something to do.
#Description
#Similar to a 'true-ish' type value. Everything but FALSE and NA should return TRUE.
#preprocess {growr}	R Documentation
# preprocess raw growth-curve measures
# Description preprocess raw growth-curve measures
#
# Usage
# preprocess(y, log_base = exp(1), runmed_k = 3, runmav_n = 5, force_inc = TRUE, bg_subtract = function(y) y - min(y), calibrate_fxn = function(y) y)
# Arguments
# y	numeric vector of values to preprocess

# log_base	<num>|<lgl> base of log for log(y / min(y)) transformation. Defaults to base [exp(1)]. set to NA or FALSE to skip.
#
# runmed_k
# integer width of median window, must be odd. Passed to stats::runmed(). Defaults to 3. set to NA or FALSE to skip.
#
# runmav_n
# integer width of moving average filter. Defaults to 5. set to NA or FALSE to skip.
#
# force_inc
# force y to be monotonically increasing?
#
#   bg_subtract
# numeric value to use for background substraction. Defaults to min(y). set to NA or FALSE to skip.
#
# calibrate_fxn
# function to calibrate y values. Usually for correcting for non-linear scattering at higher concentrations. Defaults to the identity function I, that simply returns y un-transformed (as is).


preprocess
function (y, log_base = exp(1), runmed_k = 3, runmav_n = 5, force_inc = TRUE,
  bg_subtract = function(y) y - min(y), calibrate_fxn = function(y) y)
{
  if (is_todo(log_base))
    y <- log(y/min(y), base = log_base)
  if (is_todo(runmed_k))
    y <- runmed(y, k = runmed_k)
  if (is_todo(runmav_n))
    y <- runmav(y, n = runmav_n)
  if (is_todo(force_inc))
    y <- enforce_mono_inc(y)
  if (is_todo(bg_subtract))
    y <- bg_subtract(y)
  if (is_todo(calibrate_fxn))
    y <- calibrate_fxn(y)
  y
}

gf_spline
function (time, n)
{
  smooth.spline(x = time, y = n)
}

# Description
# running mean smoothing of length n

# Usage
# runmav(y, n)
# Arguments
# y
# values to process
#
# n
# length of filter, width of window to apply mean smoothing to
#
#
summarise_fit2 = function (x, y, method = "smooth.spline"){
  data <- tibble::tibble(time = x, y = y)
  f <- switch(method,smooth.spline  = growr:::fit_smooth_spline, gompertz = growr:::fit_gompertz,
    richards = growr:::fit_richards, logistic = growr:::fit_logistic, avgG = myavgG2)
  res <- f(data)
  tibble::as_tibble(res[c("A", "mu", "lambda")])
}
### stats
> runmed
function (x, k, endrule = c("median", "keep", "constant"), algorithm = NULL,
  print.level = 0)
{
  n <- as.integer(length(x))
  if (is.na(n))
    stop(gettextf("invalid value of %s", "length(x)"), domain = NA)
  k <- as.integer(k)
  if (is.na(k))
    stop(gettextf("invalid value of %s", "'k'"), domain = NA)
  if (k < 0L)
    stop("'k' must be positive")
  if (k%%2L == 0L)
    warning(gettextf("'k' must be odd!  Changing 'k' to %d",
      k <- as.integer(1 + 2 * (k%/%2))), domain = NA)
  if (n == 0L) {
    x <- double()
    attr(x, "k") <- k
    return(x)
  }
  if (k > n)
    warning(gettextf("'k' is bigger than 'n'!  Changing 'k' to %d",
      k <- as.integer(1 + 2 * ((n - 1)%/%2))), domain = NA)
  algorithm <- if (missing(algorithm)) {
    if (k < 20L || n < 300L)
      "Stuetzle"
    else "Turlach"
  }
  else match.arg(algorithm, c("Stuetzle", "Turlach"))
  endrule <- match.arg(endrule)
  iend <- switch(endrule, median = , keep = 0L, constant = 1L)
  if (print.level)
    cat("runmed(*, endrule=", endrule, ", algorithm=", algorithm,
      ", iend=", iend, ")\n")
  res <- switch(algorithm, Turlach = .Call(C_runmed, as.double(x),
    1, k, iend, print.level), Stuetzle = .Call(C_runmed,
      as.double(x), 0, k, iend, print.level))
  if (endrule == "median")
    res <- smoothEnds(res, k = k)
  attr(res, "k") <- k
  res
}

### diffs = 1 less e.g. a1[17:21] = diffs[16:20]
###
### >
###
SummarizeGrowth
function (data_t, data_n, t_trim = 0, bg_correct = "min", blank = NA)
{
  if (is.list(data_t) == TRUE) {
    tryCatch({
      data_t <- unlist(data_t)
    }, error = function(e) {
      stop("data_t is not a vector")
    })
  }
  if (is.list(data_n) == TRUE) {
    tryCatch({
      data_n <- unlist(data_n)
    }, error = function(e) {
      stop("data_n is not a vector")
    })
  }
  stopifnot(is.vector(data_t))
  stopifnot(is.vector(data_n))
  if (!is.numeric(data_t) | !is.numeric(data_n)) {
    stop("Error: The input data (data_t and data_n) must be numeric.")
  }
  if (length(data_t) != length(data_n)) {
    stop("Error: The input data (data_t and data_n) must have the same number of rows")
  }
  if (!bg_correct %in% c("none", "min", "blank")) {
    stop(paste0(bg_correct, "is not a valid option for bg_correct"))
  }
  if (bg_correct == "blank") {
    if (is.list(blank) == TRUE) {
      tryCatch({
        blank <- unlist(blank)
      }, error = function(e) {
        stop("blank is not a vector")
      })
      stopifnot(is.vector(blank))
      if (!is.numeric(blank)) {
        stop("Error: The blank data must be numeric.")
      }
      if (length(blank) != length(data_n)) {
        stop("Error: The input data (data_n) and the background correction data (blank) must have the same number of rows.")
      }
    }
  }
  if (t_trim > 0) {
    t_max <- t_trim
    data_n <- data_n[data_t < t_trim]
    data_t <- data_t[data_t < t_trim]
    if (bg_correct == "blank") {
      blank <- blank[data_t < t_trim]
    }
  }
  else {
    t_max <- max(data_t, na.rm = TRUE)
  }
  if (bg_correct == "blank") {
    data_n <- data_n - blank
    data_n[data_n < 0] <- 0
  }
  else if (bg_correct == "min") {
    data_n <- data_n - min(data_n)
  }
  tryCatch({
    log_mod = FitLogistic(data_t, data_n)
  }, error = function(e) {
  })
  if (exists("log_mod") == FALSE) {
    vals <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  }
  else {
    p <- summary(log_mod)$coefficients
    k <- p[1]
    k_se <- p[4]
    k_p <- p[10]
    n0 <- p[2]
    n0_se <- p[5]
    n0_p <- p[11]
    r <- p[3]
    r_se <- p[6]
    r_p <- p[12]
    t_inflection <- TAtInflection(k, n0, r)
    DT <- MaxDt(r)
    sigma <- summary(log_mod)$sigma
    df <- summary(log_mod)$df[2]
    auc_l <- AreaUnderCurve(k, n0, r, 0, t_max)$value
    auc_e <- EmpiricalAreaUnderCurve(data_t, data_n, t_max)
    vals <- c(k, k_se, k_p, n0, n0_se, n0_p, r, r_se, r_p,
      sigma, df, t_inflection, DT, auc_l, auc_e)
  }
  val_names <- c("k", "k_se", "k_p", "n0", "n0_se", "n0_p",
    "r", "r_se", "r_p", "sigma", "df", "t_mid", "t_gen",
    "auc_l", "auc_e")
  vals <- stats::setNames(as.list(vals), val_names)
  class(vals) <- "gcvals"
  if (exists("log_mod") == FALSE) {
    vals$note <- "cannot fit data"
    log_mod <- ""
  }
  else if (k < n0) {
    vals$note <- "questionable fit (k < n0)"
  }
  else if (t_inflection < 0) {
    vals$note <- "questionable fit"
  }
  else {
    vals$note <- ""
  }
  ret <- list(vals = vals, model = log_mod, data = list(t = data_t,
    N = data_n))
  class(ret) <- "gcfit"
  return(ret)
}
##########
> growthcurver:::FitLogistic
function (data_t, data_n)
{
  if (!is.vector(data_t) | !is.vector(data_n)) {
    stop("Error: The input data (data_t and data_n) must be vectors.")
  }
  if (!is.numeric(data_t) | !is.numeric(data_n)) {
    stop("Error: The input data (data_t and data_n) must be numeric.")
  }
  if (length(data_t) != length(data_n)) {
    stop("Error: The input data (data_t and data_n) must have the same length.")
  }
  d <- data.frame(cbind(data_t, data_n))
  names(d) <- c("t", "n")
  k_init <- max(data_n)
  n0_init <- min(data_n[data_n > 0])
  glm_mod <- stats::glm(n/k_init ~ t, family = stats::quasibinomial("logit"),
    data = d)
  r_init <- stats::coef(glm_mod)[[2]]
  if (r_init <= 0) {
    r_init <- 0.001
  }
  suppressWarnings(nls_mod <- tryCatch(minpack.lm::nlsLM(n ~
      k/(1 + ((k - n0)/n0) * exp(-r * t)), start = list(k = k_init,
        n0 = n0_init, r = r_init), control = list(maxiter = 500),
    lower = c(stats::median(data_n), 0, 0), upper = c(Inf,
      max(data_n), Inf), data = d), error = function(e) {
        stop("Error: Growthcurver FitLogistic cannot fit data.")
      }))
  return(nls_mod)
}

