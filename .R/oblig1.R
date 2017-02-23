## ---- Preliminary
gastricXelox <- read.csv(paste("https://www.math.ntnu.no/~jarlet/levetid/",
                         "gastricXelox.csv", sep=""))
library(survival)

## ---- Task 1
time <- gastricXelox$timeWeeks
delta <- gastricXelox$delta

data <- Surv(time, delta)
R_KM <- survfit(data ~ 1)
plot(R_KM)

## ---- Task 2
unique_time_indices <- !duplicated(time)
unique_time_weeks <- time[unique_time_indices]
unique_time_event_indices <- (R_KM$n.event != 0)
unique_time_event_weeks <- unique_time_weeks[unique_time_event_indices]

Z_KM <- -log(R_KM$surv[unique_time_event_indices])
# Z_KM <- -cumsum(log(1-R_KM$n.event[unique_time_event_indices]/
#                       R_KM$n.risk[unique_time_event_indices]))
plot(unique_time_event_weeks,Z_KM)

Z_NA <- cumsum(R_KM$n.event[unique_time_event_indices]/
                 R_KM$n.risk[unique_time_event_indices])
plot(unique_time_event_weeks, Z_NA)

n <- length(time)
r <- length(R_KM$n.event[R_KM$n.event!=0])  # no. of events
Tau <- vector(mode="double", length=r)
j <- 2
Tau[1] <- Tau[2] <- n * time[1]
for (i in 2:n) {
  if (time[i] == time[i-1]) next
  Tau[j] <- Tau[j] + (n-i+1)*(time[i] - time[i-1])
  if (delta[i]) {
    if (j == r) break
    j <- j+1
    Tau[j] <- Tau[j-1]
  }
}

plot((1:r)/r, Tau/tail(Tau, n=1))

## ---- Task 3
W <- sum(Tau[1:(r-1)]/tail(Tau, n=1))
Z <- (W - (r-1)/2) / sqrt((r-1)/12)
print(pnorm(Z))

W_ex_last <- sum(Tau[1:(r-2)]/Tau[r-1])
Z_ex_last <- (W_ex_last - (r-2)/2) / sqrt((r-2)/12)
F_Z_ex_last <- pnorm(Z_ex_last)

## ---- Task 4
print(R_KM, print.rmean=TRUE)
mean_KM <- 93.9
SE_KM <- 15.6
median_KM <- 44.5

s <- sum(time)
theta_hat <- s/r
median_theta_hat <- qexp(0.5, rate=1/theta_hat)
SD_theta_hat <- theta_hat/sqrt(r)
CI_95_ln <- theta_hat * c(exp(-1.96*SD_theta_hat/theta_hat), 
                      exp(1.96*SD_theta_hat/theta_hat))

plot(R_KM)
curve(1 - pexp(x, rate=1/theta_hat), from=0, to=253, add=TRUE)

## ---- Task 5
l <- function(theta, r, s) { return (-r*log(theta) - s/theta) }
curve(l(x, r, s), from=50, to=200)
abline(v=theta_hat)
abline(h = l(theta_hat, r, s) - pchisq(0.95, df=1)/2)
left <- uniroot( function(x) 2*(l(theta_hat, r, s) - l(x, r, s)) - pchisq(0.95, df=1)/2, 
                lower = 50, upper = 114)
right <- uniroot( function(x) 2*(l(theta_hat, r, s) - l(x, r, s)) - pchisq(0.95, df=1)/2, 
                 lower = 114, upper = 200)
CI_95_chisq = c(left$root, right$root)

## ---- Task 6
plot(R_KM)

weibfit <- survreg(data ~ 1, dist="weibull")
lls_scale = weibfit$scale
intercept = weibfit$coefficients["(Intercept)"]
curve(1-pweibull(x, shape=1/lls_scale, scale=exp(intercept)), add=TRUE)

lnfit <- survreg(data ~ 1, dist="lognormal")
intercept = lnfit$coefficients["(Intercept)"]
curve(1-plnorm(x, meanlog=intercept, sdlog=lnfit$scale), add=TRUE)

expfit = survreg(data ~ 1, dist="exponential")
intercept = expfit$coefficients["(Intercept)"]
curve(1-pexp(x, rate=exp(-intercept)), add=TRUE)

# library(FAdist)
# llogfit <- survreg(data ~ 1, dist="loglogistic")
# intercept = llogfit$coefficients["(Intercept)"]
# curve(1-pllog(x, shape=intercept, scale=llogfit$scale), add=TRUE)

