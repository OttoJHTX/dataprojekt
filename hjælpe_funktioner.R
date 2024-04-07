### FILL WITH 0's ###

fill_gaps <- function(df) {
  # Initialize an empty list to store the new rows
  new_rows <- list()
  for(i in 1:(nrow(df) - 1)) {
    current_position <- df$start[i]
    next_position <- df$start[i + 1]
    if(next_position - current_position > 1) {
      # Generate the missing positions
      missing_positions <- (current_position + 1):(next_position - 1)
      for(missing_pos in missing_positions) {
        new_row <- df[i, ]
        new_row$start <- missing_pos
        new_row$score <- 0 # Set score to 0 for the new row
        new_rows[[length(new_rows) + 1]] <- new_row
      }
    }
  }	
  
  if(length(new_rows) > 0) {
    new_rows_df <- do.call(rbind, new_rows)
    # Combine with the original dataframe and sort
    df <- rbind(df, new_rows_df)
    df <- df[order(df$start), ]
  }
  
  return(df)
}


# Ny

gap_fill <- function(df) {
  final_df <- data.frame(head(df, 0L)) %>% select(start, score, transcript_number)
  for (transcript in unique(df$transcript_number)) {
    current_df <- df[df$transcript_number == transcript,]
    rng <- min(current_df$start):max(current_df$start)
    zeroes <- data.frame(start = rng, score = rep(0, length(rng)), transcript_number = transcript)
    current_df <- zeroes %>% 
      left_join(current_df, by = "start", keep = FALSE) %>% 
      mutate(score = ifelse(is.na(score.y), score.x, score.y)) %>%
      mutate(transcript_number = ifelse(is.na(transcript_number.y), transcript_number.x, transcript_number.y)) %>% 
      select(start, transcript_number, score)
    final_df <- rbind(final_df, current_df)
  }
  
  return (final_df)
}


### AGGREGATE DF WITH BINSIZE ###

create_aggregated_df <- function(df, bin_size) {
  # Determine the number of chunks
  n <- nrow(df)
  num_chunks <- ceiling(n / bin_size)
  
  
  
  medians <- numeric(num_chunks)
  means <- numeric(num_chunks)
  
  for(i in 1:num_chunks) {
    # Calculate indices for slicing the dataframe
    start_index <- (i - 1) * 100 + 1
    end_index <- min(i * 100, n)
    
    chunk <- df[start_index:end_index, ]
    medians[i] <- median(chunk$start)
    means[i] <- mean(chunk$score)
  }	
  
  
  aggregated_df <- data.frame(start = medians, score = means)
  
  return(aggregated_df)
}

### Fitter model ###

mse <- function(model, y_observed, root = F, y_predicted = NULL) {
  if (is.null(y_predicted)) {
    y_predicted <- predict(model)
  }
  r <- ifelse(root, 1, 2)
  output <- mean((y_observed - y_predicted)^r)
  return (output)
}

r_sq <- function(model, y_observed, y_predicted = NULL) {
  if (is.null(y_predicted)) {
    y_predicted <- predict(model)
  }
  
  # Calculate mean values
  mean_observed <- mean(y_observed)
  mean_predicted <- mean(y_predicted)
  
  # Calculate sums for the formula
  numerator <- sum((y_observed - mean_observed) * (y_predicted - mean_predicted))
  denominator_observed <- sum((y_observed - mean_observed)^2)
  denominator_predicted <- sum((y_predicted - mean_predicted)^2)
  
  # Calculate correlation coefficient (r)
  r <- numerator / sqrt(denominator_observed * denominator_predicted)
  
  # Calculate R-squared
  return (r^2)
  
}

# DOUBLE SIGMOID

double_sigmoid <- function(x, A1, B1, C1, A2, B2, C2) {
  A1 / (1 + exp((B1 - x) / C1)) + A2 / (1 + exp((B2 - x) / C2))
}
