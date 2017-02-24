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
event_n <- length(R_KM$n.event[R_KM$n.event!=0])  # no. of events

Tau <- vector(mode="double", length=event_n)
j <- 2
Tau[1] <- Tau[2] <- n * time[1]

for (i in 2:n) {
  if (time[i]!= time[i-1]) {
    Tau[j] <- Tau[j] + (n-i+1)*(time[i] - time[i-1])
  }
  if (time[i] == unique_time_event_weeks[j]) {
    if (j == event_n) {
      print('bla')
      break
    }
    j <- j+1
    Tau[j] <- Tau[j-1]
  }
}

plot((1:event_n)/event_n, Tau/tail(Tau, n=1))

## ---- Task 3
W <- sum(Tau[1:(event_n-1)]/tail(Tau, n=1))
Z <- (W - (event_n-1)/2) / sqrt((event_n-1)/12)
print(pnorm(Z))

W_ex_last <- sum(Tau[1:(event_n-2)]/Tau[event_n-1])
Z_ex_last <- (W_ex_last - (event_n-2)/2) / sqrt((event_n-2)/12)
F_Z_ex_last <- pnorm(Z_ex_last)

## ---- Task 4
print(R_KM, print.rmean=TRUE)
mean_KM <- 93.9
SE_KM <- 15.6
median_KM <- 44.5

s <- sum(time)
r <- sum(delta)
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
left <- uniroot( function(x) 2*(l(theta_hat, r, s) - l(x, r, s)) - pchisq(0.95, df=1), 
                lower = 50, upper = 89)
right <- uniroot( function(x) 2*(l(theta_hat, r, s) - l(x, r, s)) - pchisq(0.95, df=1), 
                 lower = 89, upper = 150)
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

library(FAdist)
llogfit <- survreg(data ~ 1, dist="loglogistic")
intercept = llogfit$coefficients["(Intercept)"]
curve(1-pllog(x, shape=1.5, scale=4), add=TRUE)

## Task 7


## spm
# intercept, best fit
# tolkning av levetider
# CI median
# message feb 17
