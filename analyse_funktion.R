library(rtracklayer)
library(GenomicRanges)
library(tidyverse)

getwd()
setwd("C:/Users/sejeo/Documents/Uni/4. semester/Dataprojekt")
# e.g. analyse(ctrl_path, test_path, annotation_path, "chr12", "GAPDH")
analyse <- function(ctrl_path, test_path, annotation_path, chr, gene_name) {
  
  cps <- import(test_path)
  ctrl <- import(ctrl_path)
  annot_gr <- import(annotation_path)
  
  ts_annot_gr <- subset(annot_gr, type == "transcript" & gene_name == gene_name)
  ctrl <- subset(ctrl, seqnames == chr)
  cps <- subset(cps, seqnames == chr)

  ### Create DF with transcript no. column, position (start & end), score etc. ###
  
  subsetOverlaps_ctrl <- as.list(rep(0,11))
  subsetOverlaps_cps <- as.list(rep(0,11))
  
  for (i in seq_along(ts_annot_gr)) {
    current_observation <- ts_annot_gr[i]
    
    start_coords <- start(current_observation) - 500
    end_coords <- end(current_observation) + 100000
    
    new_ranges <- GRanges(
      seqnames = chr,
      ranges = IRanges(start = start_coords, end = end_coords)
    )
    
    subsetOverlaps_ctrl[i] <- (subsetByOverlaps(ctrl, new_ranges))
    subsetOverlaps_cps[i] <- (subsetByOverlaps(cps, new_ranges))
  }
  
  ### Create an empty data frame to store the results ###
  ctrl_df <- data.frame(transcript_number = numeric(),
                        seqnames = character(),
                        start = numeric(),
                        end = numeric(),
                        score = numeric(),
                        stringsAsFactors = FALSE)
  
  cps_df <- data.frame(transcript_number = numeric(),
                        seqnames = character(),
                        start = numeric(),
                        end = numeric(),
                        score = numeric(),
                        stringsAsFactors = FALSE)
  
  ctrl_overlaps <- subsetOverlaps_ctrl
  cps_overlaps <- subsetOverlaps_cps
  
  ### Fill empty df ###
  
  for (i in seq_along(ctrl_overlaps)) {
    # Extract the current GRanges object
    current_gr <- ctrl_overlaps[[i]]
    
    # Extract range and score values from the current GRanges object
    range_values <- as.data.frame(current_gr)
    
    # Add transcript number column
    range_values$transcript_number <- i
    
    # Append range and score values to the result data frame
    ctrl_df <- rbind(ctrl_df, range_values)
  }
  
  ctrl_df$transcript_number <- as.factor(ctrl_df$transcript_number)
  
  for (i in seq_along(cps_overlaps)) {
    # Extract the current GRanges object
    current_gr <- cps_overlaps[[i]]
    
    # Extract range and score values from the current GRanges object
    range_values <- as.data.frame(current_gr)
    
    # Add transcript number column
    range_values$transcript_number <- i
    
    # Append range and score values to the result data frame
    cps_df <- rbind(cps_df, range_values)
  }
  
  cps_df$transcript_number <- as.factor(cps_df$transcript_number)
  
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
  
  # USAGE
  
  ctrl_filled <- gap_fill(ctrl_df)
  cps_filled <- gap_fill(cps_df)
  
  length(ctrl_filled$start[ctrl_filled$transcript_number == 10])
  
  max(ctrl_filled$start[ctrl_filled$transcript_number == 10]) - min(ctrl_filled$start[ctrl_filled$transcript_number == 10])
  
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
  
  # USAGE
  
  ctrl_agg <- create_aggregated_df(ctrl_df, 100)
  cps_agg <- create_aggregated_df(cps_df, 100)
  
  ### Find kroppen ###
  
  # Get diffs:
  diffs_ctrl <- c(0, diff(ctrl_agg$score))
  ctrl_agg$diff <- diffs_ctrl 
  
  diffs_cps <- c(0, diff(cps_agg$score))
  cps_agg$diff <- diffs_cps 
  
  # Getting start and termination position:
  krop_ctrl <- data.frame(start = subset(ctrl_agg, diff == max(diffs_ctrl))$start - 100, slut = subset(ctrl_agg, diff == min(diffs_ctrl))$start)
  krop_cps <- data.frame(start = subset(cps_agg, diff == max(diffs_cps))$start - 100, slut = subset(cps_agg, diff == min(diffs_cps))$start)
  
  ### Plots ###
  
  ctrl_mean <- mean(ctrl_filled$score[ctrl_filled$start > krop_ctrl$start & ctrl_filled$start < krop_ctrl$slut])
  cps_mean <- mean(cps_filled$score[cps_filled$start > krop_cps$start & cps_filled$start < krop_cps$slut])
  
  
  plot_ctrl <- subset(ctrl_filled, start < 6545000) %>% mutate(score = score / ctrl_mean)
  model <- lm(data = filter(plot_ctrl, start > krop_ctrl$slut), score ~ start)
  
  plot(data = plot_ctrl, score ~ start)
  abline(v = krop_ctrl$start, lty = "dashed")
  abline(v = krop_ctrl$slut, lty = "dashed")
  abline(model, col = "red")
  
  plot_cps <- subset(cps_filled, start < 6545000) %>% mutate(score = score / cps_mean)
  
  ggplot(plot_ctrl) + geom_line(aes(start, score)) + geom_vline(xintercept = krop_ctrl$start, linetype = "dashed") + geom_vline(xintercept = krop_ctrl$slut, linetype = "dashed") #+ geom_hline(yintercept = -ctrl_mean)
  ggplot(plot_cps) + geom_line(aes(start, score)) + geom_vline(xintercept = krop_cps$start, linetype = "dashed") + geom_vline(xintercept = krop_cps$slut, linetype = "dashed")#+ geom_hline(yintercept = -cps_mean)
  
  b <- left_join(plot_cps, plot_ctrl, by = "start") %>% pivot_longer(c("score.x", "score.y"), names_to = "Data", values_to = "Score")
  
  ggplot(b, aes(start, Score, color = Data)) + geom_line() + geom_line() + ggtitle("Coverage for test data and control data") + ylab("norm coverage") + xlab("position (bp)") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  ggplot() + geom_line(plot_ctrl, mapping = aes(start, score), color = "red") + geom_line(plot_cps, mapping = aes(start, score), color = "blue") + ggtitle("Coverage for test data and control data") + ylab("norm coverage") + xlab("position (bp)") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  ### Statistisk test for om kroppen er ens for ctrl og test ###
  
  # Ens varians
  
  a <- subset(ctrl_agg, start %in% krop_ctrl$start:krop_ctrl$slut) %>% 
    mutate(score = score - mean(.$score))
  plot(a$start, a$score)
  
  b <- subset(cps_agg, start %in% krop_cps$start:krop_cps$slut) %>% 
    mutate(score = score - mean(.$score))
  plot(b$start, b$score)
  
  ggplot() + geom_line(a, mapping = aes(start, score), color = "red") + geom_line(b, mapping = aes(start, score), color = "blue")
  
  var.test(a$score, b$score)
  
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
  
  library(minpack.lm)
  
  ## Only test
  
  temp_ctrl <- subset(plot_ctrl, start > krop_ctrl$slut)
  temp_ctrl <- subset(temp_ctrl, transcript_number == 10)
  
  # temp <- head(create_aggregated_df(temp, nrow(temp) / 100))
  
  x <- temp_ctrl$start
  y <- temp_ctrl$score
  x2 <- 1:length(x)
  
  # Create a data frame
  data <- data.frame(x = x, y = y, x2 = x2)
  # Adjust starting values
  
  # Polynomial 
  
  poly2_model <- lm(y ~ poly(x2, 2))
  
  mse(model = poly2_model, y_observed = y, root = F)
  r_sq(poly2_model, y)
  
  poly4_model <- lm(data = data, y ~ poly(x2, 4))
  
  mse(poly4_model, y, root = F)
  r_sq(poly4_model, y)
  
  # Exponential
  
  exp_model <- nlsLM(y ~ b * exp(-c * x2), data = data, start = list(b = 1, c = 0.01))
  
  mse(exp_model, y, root = F)
  r_sq(exp_model, y)
  
  # Gaussian
  
  gaus_model <- nlsLM(y ~ h * dnorm(x2, mean = m, sd = v), data = data, start = list(h = 10, m = 0, v = 1000))
  
  mse(gaus_model, y, root = F)
  r_sq(gaus_model, y)
  
  # Plot the data
  plot(data$x2, data$y, main = "Fitted Gaussian", xlab = "bp after body end", ylab = "norm cov")
  x_seq <- seq(min(data$x2), max(data$x2))
  y_pred <- predict(gaus_model, newdata = data.frame(x2 = x_seq))
  lines(x_seq, y_pred, col = "red", lwd = 2)
  
  # fitted model on test data
  temp_cps <- subset(plot_cps, start > krop_cps$slut)
  temp_cps <- subset(temp_cps, transcript_number == 10)
  
  data_cps <- data.frame(x = temp_cps$start, y = temp_cps$score, x2 = 1:length(temp_cps$start))
  
  plot(data_cps$x2, data_cps$y, main = "Exponential Model on Test Data", xlab = "bp after body end", ylab = "norm cov")
  x_seq <- seq(min(data_cps$x2), max(data_cps$x2))
  y_pred <- predict(exp_model, newdata = data.frame(x2 = x_seq))
  lines(x_seq, y_pred, col = "red", lwd = 2)
  
  mse(exp_model, data_cps$y, y_predicted = y_pred)
  r_sq(exp_model, data_cps$y, y_predicted = y_pred)
  
  ## Test and control
  
  # Full
  data_diff <- left_join(data, data_cps, by = "x2", suffix = c(".ctrl", ".cps")) %>% 
    mutate(y_diff = y.cps - y.ctrl) %>% 
    pivot_longer(c("y.ctrl", "y.cps"), names_to = "Data", values_to = "y")
  
  ggplot(data_diff, aes(x2, y, color = Data)) + geom_point()
  
  # Diff
  
  ggplot(data_diff, aes(x2, y_diff)) + geom_point() + geom_line(aes(x2, y, color = Data))
  
  double_sigmoid <- function(x, A1, B1, C1, A2, B2, C2) {
    A1 / (1 + exp((B1 - x) / C1)) + A2 / (1 + exp((B2 - x) / C2))
  }
  
  initial_guess <- list(A1 = 0.1, B1 = 0, C1 = 1, A2 = 0.1, B2 = 0, C2 = 1)  # Initial parameter estimates
  
  dbl_sigmoid_model <- nlsLM(data = data_diff, y_diff ~ double_sigmoid(x2, A1, B1, C1, A2, B2, C2), start = initial_guess)
  
  summary(dbl_sigmoid_model)
  
  mse(dbl_sigmoid_model, data_diff$y_diff)
  r_sq(dbl_sigmoid_model, data_diff$y_diff)
  
  plot(data_diff$x2, data_diff$y_diff, main = "Fitted Double Sigmoid", xlab = "bp after body end", ylab = "norm cov")
  x_seq <- seq(min(data_diff$x2), max(data_diff$x2))
  y_pred <- predict(dbl_sigmoid_model, newdata = data.frame(x2 = x_seq))
  lines(x_seq, y_pred, col = "red", lwd = 2)
}
