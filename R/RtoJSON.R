install.packages("jsonlite")
library(jsonlite)

# Create sample data
values <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Create a list with keys and values
data <- list()
for (i in 1:length(values)) {
  key <- paste0("key", i)
  data[[key]] <- values[i]
}

# Convert list to JSON
json_data <- jsonlite::toJSON(data, auto_unbox = TRUE)

# Write JSON to file
write(json_data, file = "output.json")
