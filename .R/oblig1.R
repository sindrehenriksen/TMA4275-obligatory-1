## ---- Preliminary
gastricXelox <- read.csv(paste("https://www.math.ntnu.no/~jarlet/levetid/",
                         "gastricXelox.csv", sep=""))
library(survival)

## ---- Task 1
data <- Surv(gastricXelox$timeWeeks, gastricXelox$delta)
R_KM <- survfit(data ~ 1)
plot(R_KM)

## ---- Task 2
unique_time_indices = !duplicated(gastricXelox$timeWeeks)
unique_time_weeks = gastricXelox$timeWeeks[unique_time_indices]

Z_KM = -log(R_KM$surv)
plot(unique_time_weeks,Z_KM)

Z_NA <- cumsum(R_KM$n.event/R_KM$n.risk)
plot(unique_time_weeks, Z_NA)

