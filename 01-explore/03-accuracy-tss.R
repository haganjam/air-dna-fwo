
# simnulate the detection probability

# define parameters
n_sites <- 100               # total number of sites
n_species_A <- 50            # sites with species A present
false_negative_rate <- 0.10  # false negative rate
false_positive_rate <- 0.05  # false positive rate

# 1. generate true presence/absence data for species A
true_presence <- rep(0, n_sites)
true_presence[sample(1:n_sites, n_species_A)] <- 1

# 2. generate observed detection data based on true presence and error rates
observed_detection <- sapply(true_presence, function(presence) {
  if (presence == 1) {
    # species is actually present, consider false negative rate
    rbinom(1, 1, 1 - false_negative_rate)
  } else {
    # species is not present, consider false positive rate
    rbinom(1, 1, false_positive_rate)
  }
})

# combine into a data frame
data <- data.frame(
  site = 1:n_sites,
  true_presence = true_presence,
  observed_detection = observed_detection
)

# display the first few rows of the dataset
head(data)

# generate the confusion matrix

# calculate components of the confusion matrix
TP <- sum(data$true_presence == 1 & data$observed_detection == 1) # True Positives
FN <- sum(data$true_presence == 1 & data$observed_detection == 0) # False Negatives
FP <- sum(data$true_presence == 0 & data$observed_detection == 1) # False Positives
TN <- sum(data$true_presence == 0 & data$observed_detection == 0) # True Negatives

# calculate N
N <- (TP + FN + FP + TN)

# create the confusion matrix
confusion_matrix <- matrix(c(TP, FN, FP, TN), nrow = 2, byrow = TRUE,
                           dimnames = list("Observed" = c("Detected", "Not Detected"),
                                           "Actual" = c("Present", "Absent")))
print(confusion_matrix)

# calculate different metrics of accuracy, sensitivity, specificity, TSS

# accuracy
ac <- (TP + TN)/N

# sensitivity
se <- TP/(TP + FN)

# specificity
sp <- TN/(FP + TN)

# TSS
(se + sp) - 1

