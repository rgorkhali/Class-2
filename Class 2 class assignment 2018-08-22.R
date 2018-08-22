#class excercise 1.

#Rakshya Gorkkhali, 08/22/2018, this program does the KM analysis and plot

#dim9survData) gives dimensions of your dataset, Intercept?
# Load the dataset, survival package which has function for KM analysis

require(datasets)
require(survival)

# Check your dataset
head(lung)

# Make an object survData
survData <-lung

# use package survival
survObj <- Surv(time=survData$time, event=survData$status==2, type='right')

# Fit your plot to the intercept
fit <- survfit(survObj~1)

# Plot the fit
plot(fit, main="KM plot for Lung cancer data", xlab="time in days", ylab="Percent survival")

#Fit your plot to the sexes
fit <- survfit(survObj ~ survData$sex==1)

# Plot the new fit
plot(fit, main="KM plot for Lung cancer data in males and females", xlab="time in days", ylab="Percent survival", col=c(1,2))

#assigning labels
legend('topright', c("Male","Female"), lty=1)