## ---- Preliminary
gastricXelox <- read.csv(paste("https://www.math.ntnu.no/~jarlet/levetid/",
                         "gastricXelox.csv", sep=""))
library(survival)

## ---- Task 1
time = gastricXelox$timeWeeks
delta = gastricXelox$delta

data <- Surv(time, delta)
R_KM <- survfit(data ~ 1)
plot(R_KM)

## ---- Task 2
unique_time_indices <- !duplicated(time)
unique_time_weeks <- time[unique_time_indices]

Z_KM <- -log(R_KM$surv)
# Z_KM <- -cumsum(log(1-R_KM$n.event/R_KM$n.risk))
plot(unique_time_weeks,Z_KM)

Z_NA <- cumsum(R_KM$n.event/R_KM$n.risk)
plot(unique_time_weeks, Z_NA)

n <- length(time)
event_n <- length(R_KM$n.event[R_KM$n.event!=0])
Tau <- vector(mode="double", length=event_n)
j <- 2
Tau[1] <- Tau[2] <- n * time[1]
for (i in 2:n) {
  if (time[i] == time[i-1]) next
  Tau[j] <- Tau[j] + (n-i+1)*(time[i] - time[i-1])
  if (delta[i]) {
    if (j == event_n) break
    j <- j+1
    Tau[j] <- Tau[j-1]
  }
}

## ---- Task 3

