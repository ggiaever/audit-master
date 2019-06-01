

#A	Maximum growth value.
#mu	Maximum slope.
#lambda	Lag-phase.
# myavgG2 = function(data){
#
#   mycalib = function(od){
#     nod = (od-0.0672)/0.2973
#     nod
#   }
#
#   mylag = function(vect){
#   diffs = diff(vect)
#   sign = sign(diffs)
#   rign = rle(sign)
#   wrign=which(rign$lengths == max(rign$lengths) & rign$values == 1)
#   if(length(wrign) > 0) {
#
#   lens = 1:(wrign-1)
#
#   sm = sum(rign$lengths[lens])
#   sm = sm + 1
#   }
#
#   sm
#
#   }
#
#   x2 = mycalib(data$y)
#
#   xlag = mylag(x2)
#
#   wx = which(x2 >= 2)
#
#   if(length(wx) > 0) wx = wx[1]
#
#   tx = data$time[wx] - data$time[xlag]
#
#   if(xlag == 0) G = 0 else G = tx/5
#
#   if(G > 10) G = G/3600
#
#   init_smooth_spline = function (x, y, ...)
#   {
#     fit <- smooth.spline(x, y, cv = TRUE)
#     deriv1 <- predict(fit, deriv = 1)
#     slope_max_pos <- which.max(deriv1$y)
#     x_at_slope_max <- x[slope_max_pos]
#     y_at_slope_max <- y[slope_max_pos]
#     dy_at_slope_max <- deriv1$y[slope_max_pos]
#     f <- function(x) {
#       p <- predict(fit, x)
#       p$y
#     }
#     mu = dy_at_slope_max
#     b = y_at_slope_max - (dy_at_slope_max * x_at_slope_max)
#     lambda = -b/mu
#     integral <- integrate(f, min(x), max(x))$value
#     list(mu = mu, b = b, lambda = lambda, A = max(y), integral = integral,
#       fit = fit)
#   }
#
#   init = init_smooth_spline(x = data$time,y = data$y)
#
#   init$mu = G
#
#
#
#   init$lambda = xlag
#
#   init
# }

AvgG = function(data){

  delta = 0.66

  time = data$time

  y = data$y

  y = y - min(y)

  wy = which(y >= delta)

  if(length(wy) > 0) wy = wy[1]

  ty = data$time[wy]

  if (length(wy) != 0) G = ty/5 else G = 0

  if(G > 10) G = round(G/3600,2)

  init_smooth_spline = function (x, y, ...)
  {
    fit <- smooth.spline(x, y, cv = TRUE)
    deriv1 <- predict(fit, deriv = 1)
    slope_max_pos <- which.max(deriv1$y)
    x_at_slope_max <- x[slope_max_pos]
    y_at_slope_max <- y[slope_max_pos]
    dy_at_slope_max <- deriv1$y[slope_max_pos]
    f <- function(x) {
      p <- predict(fit, x)
      p$y
    }
    mu = dy_at_slope_max
    b = y_at_slope_max - (dy_at_slope_max * x_at_slope_max)
    lambda = -b/mu
    integral <- integrate(f, min(x), max(x))$value
    list(mu = mu, b = b, lambda = lambda, A = max(y), integral = integral,
      fit = fit)
  }

  init = init_smooth_spline(x = data$time,y = data$y)

  init$mu = G

  #init$lambda = xlag

  init
}





#A	Maximum growth value.
#mu	Maximum slope.
#lambda	Lag-phase.
summarise_fit2 = function (x, y, method = "smooth.spline"){
  data <- tibble::tibble(time = x, y = y)
  f <- switch(method,smooth.spline  = growr:::fit_smooth_spline, gompertz = growr:::fit_gompertz,
    richards = growr:::fit_richards, logistic = growr:::fit_logistic, avgG = AvgG)
  res <- f(data)
  tibble::as_tibble(res[c("A", "mu", "lambda")])
}


###
fit_smooth_spline = function (data)
{
  init_smooth_spline(data$time, data$y)
}
#####
init_smooth_spline = function (x, y, ...)
{
  fit <- smooth.spline(x, y, cv = TRUE)
  deriv1 <- predict(fit, deriv = 1)
  slope_max_pos <- which.max(deriv1$y)
  x_at_slope_max <- x[slope_max_pos]
  y_at_slope_max <- y[slope_max_pos]
  dy_at_slope_max <- deriv1$y[slope_max_pos]
  f <- function(x) {
    p <- predict(fit, x)
    p$y
  }
  mu = dy_at_slope_max
  b = y_at_slope_max - (dy_at_slope_max * x_at_slope_max)
  lambda = -b/mu
  integral <- integrate(f, min(x), max(x))$value
  list(mu = mu, b = b, lambda = lambda, A = max(y), integral = integral,
    fit = fit)
}

####
preprocess = function (y, log_base = exp(1), runmed_k = 3, runmav_n = 5, force_inc = TRUE,
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
