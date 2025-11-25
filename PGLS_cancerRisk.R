# Load necessary packages
library(ape)
library(nlme)
library(geiger)
library(ggplot2)
library(ggpmisc)
library(phytools)

# Read data from file
data <- read.table("cancerRisk.csv", header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Read phylogenetic tree
tree <- read.tree("tree2.nwk")

# Data preprocessing
# Remove rows with NA CMR values
data <- data[!is.na(data$CMR), ]

# Check if all species in data are in the tree
species_in_data <- data$Species
species_in_tree <- tree$tip.label

# Find species not in the tree
missing_species <- setdiff(species_in_data, species_in_tree)
if(length(missing_species) > 0) {
  print(paste("The following species are not in the phylogenetic tree:", paste(missing_species, collapse=", ")))
  # Remove species not in the tree from data
  data <- data[data$Species %in% species_in_tree, ]
}

# Ensure tree species order matches data
tree <- keep.tip(tree, data$Species)

# Check variable distributions to decide on transformations
print("Summary of variables:")
print(summary(data[, c("Species_body_mass_kg", "Adult_life_expectancy_days", "CMR")]))

# Create improved model comparison function
compare_models <- function(formula, data, tree) {
  # Ensure no NA values in data
  model_data <- model.frame(formula, data)
  if(any(is.na(model_data))) {
    print("Warning: NA values found in data, removing these rows")
    complete_cases <- complete.cases(model_data)
    data <- data[complete_cases, ]
    tree <- keep.tip(tree, data$Species)
  }
  
  # Prepare model list
  models <- list()
  model_names <- c()
  
  # Brownian Motion model
  model_bm <- tryCatch({
    gls(formula, 
        correlation = corBrownian(phy = tree, form = ~Species),
        data = data, method = "ML")
  }, error = function(e) {
    print(paste("Brownian Motion model error:", e$message))
    return(NULL)
  })
  
  if(!is.null(model_bm)) {
    models[["bm"]] <- model_bm
    model_names <- c(model_names, "bm")
  }
  
  # Ornstein-Uhlenbeck model - using more stable method
  # Try multiple initial values, including larger ones
  ou_initial_values <- c(0.1, 0.5, 1.0, 2.0, 5.0)
  model_ou <- NULL
  
  for(alpha in ou_initial_values) {
    model_ou <- tryCatch({
      gls(formula, 
          correlation = corMartins(alpha, phy = tree, form = ~Species, fixed = FALSE),
          data = data, method = "ML",
          control = glsControl(msMaxIter = 200, msVerbose = FALSE, opt = "optim"))
    }, error = function(e) {
      return(NULL)
    })
    
    if(!is.null(model_ou)) {
      print(paste("OU model successfully fitted, initial alpha =", alpha))
      break
    }
  }
  
  if(!is.null(model_ou)) {
    models[["ou"]] <- model_ou
    model_names <- c(model_names, "ou")
  } else {
    print("All OU model initial values failed")
  }
  
  # Pagel's lambda model
  model_lambda <- tryCatch({
    gls(formula, 
        correlation = corPagel(0.5, phy = tree, form = ~Species, fixed = FALSE),
        data = data, method = "ML",
        control = glsControl(msMaxIter = 200, msVerbose = FALSE))
  }, error = function(e) {
    print(paste("Lambda model error:", e$message))
    return(NULL)
  })
  
  if(!is.null(model_lambda)) {
    models[["lambda"]] <- model_lambda
    model_names <- c(model_names, "lambda")
  }
  
  if(length(models) == 0) {
    stop("All models failed")
  }
  
  # Calculate AIC values
  aic_values <- sapply(models, AIC)
  aic_table <- data.frame(
    Model = names(models),
    AIC = aic_values,
    stringsAsFactors = FALSE
  )
  aic_table$deltaAIC <- aic_table$AIC - min(aic_table$AIC)
  aic_table$weight <- exp(-0.5 * aic_table$deltaAIC) / sum(exp(-0.5 * aic_table$deltaAIC))
  
  # Select best model
  best_model_idx <- which.min(aic_table$AIC)
  best_model_name <- aic_table$Model[best_model_idx]
  best_model <- models[[best_model_name]]
  
  # Return results
  return(list(
    models = models,
    aic_table = aic_table,
    best_model = best_model,
    best_model_name = best_model_name,
    data = data,
    tree = tree
  ))
}

# Function to calculate adjusted R² for GLS models
calculate_adjusted_r2 <- function(model, data, y_var) {
  n <- nrow(data)
  p <- length(coef(model)) - 1  # number of predictors (excluding intercept)
  
  # Calculate R²
  y <- data[[y_var]]
  y_pred <- predict(model)
  ss_res <- sum((y - y_pred)^2)
  ss_tot <- sum((y - mean(y))^2)
  r2 <- 1 - (ss_res / ss_tot)
  
  # Calculate adjusted R²
  adj_r2 <- 1 - (1 - r2) * (n - 1) / (n - p - 1)
  
  return(list(r2 = r2, adj_r2 = adj_r2))
}

# Analysis 1: CMR vs Adult life expectancy
print("Analysis 1: CMR vs Adult Life Expectancy")

# Test different transformations for both variables
transformations <- list(
  list(x_trans = "none", y_trans = "none", 
       x_formula = CMR ~ Adult_life_expectancy_days,
       x_lab = "Adult Life Expectancy (days)"),
  list(x_trans = "log", y_trans = "none", 
       x_formula = CMR ~ log10(Adult_life_expectancy_days),
       x_lab = "Adult Life Expectancy (days, log scale)"),
  list(x_trans = "none", y_trans = "log", 
       x_formula = log10(CMR) ~ Adult_life_expectancy_days,
       x_lab = "Adult Life Expectancy (days)"),
  list(x_trans = "log", y_trans = "log", 
       x_formula = log10(CMR) ~ log10(Adult_life_expectancy_days),
       x_lab = "Adult Life Expectancy (days, log scale)")
)

best_aic1 <- Inf
best_model1 <- NULL
best_trans1 <- NULL

for(trans in transformations) {
  tryCatch({
    model_result <- compare_models(trans$x_formula, data, tree)
    current_aic <- min(model_result$aic_table$AIC)
    
    if(current_aic < best_aic1) {
      best_aic1 <- current_aic
      best_model1 <- model_result
      best_trans1 <- trans
    }
    
    print(paste("Transformation (x:", trans$x_trans, ", y:", trans$y_trans, 
                ") AIC:", current_aic))
  }, error = function(e) {
    print(paste("Failed for transformation (x:", trans$x_trans, ", y:", trans$y_trans, 
                "):", e$message))
  })
}

if(!is.null(best_model1)) {
  model_results1 <- best_model1
  print(paste("Best transformation for Analysis 1: x =", best_trans1$x_trans, 
              ", y =", best_trans1$y_trans))
  print(model_results1$aic_table)
  
  # Update data and tree
  data1 <- model_results1$data
  tree1 <- model_results1$tree
  
  # Get statistics for best model
  summary1 <- summary(model_results1$best_model)
  
  # Calculate R² and adjusted R²
  if(best_trans1$y_trans == "none") {
    r2_values1 <- calculate_adjusted_r2(model_results1$best_model, data1, "CMR")
  } else {
    # For log-transformed CMR, we need to calculate R² differently
    y_pred <- predict(model_results1$best_model)
    y_actual <- log10(data1$CMR)
    ss_res <- sum((y_actual - y_pred)^2)
    ss_tot <- sum((y_actual - mean(y_actual))^2)
    r2_values1 <- list(
      r2 = 1 - (ss_res / ss_tot),
      adj_r2 = 1 - (1 - (1 - (ss_res / ss_tot))) * (nrow(data1) - 1) / (nrow(data1) - 1 - 1)
    )
  }
  
  r2_1 <- r2_values1$r2
  adj_r2_1 <- r2_values1$adj_r2
  p_value1 <- summary1$tTable[2, 4]
} else {
  stop("No valid model found for Analysis 1")
}

# Analysis 2: CMR vs Species body mass
print("Analysis 2: CMR vs Species Body Mass")

# Test different transformations for both variables
transformations2 <- list(
  list(x_trans = "log", y_trans = "none", 
       x_formula = CMR ~ log10(Species_body_mass_kg),
       x_lab = "Species Body Mass (kg, log scale)"),
  list(x_trans = "none", y_trans = "none", 
       x_formula = CMR ~ Species_body_mass_kg,
       x_lab = "Species Body Mass (kg)"),
  list(x_trans = "log", y_trans = "log", 
       x_formula = log10(CMR) ~ log10(Species_body_mass_kg),
       x_lab = "Species Body Mass (kg, log scale)"),
  list(x_trans = "none", y_trans = "log", 
       x_formula = log10(CMR) ~ Species_body_mass_kg,
       x_lab = "Species Body Mass (kg)")
)

best_aic2 <- Inf
best_model2 <- NULL
best_trans2 <- NULL

for(trans in transformations2) {
  tryCatch({
    model_result <- compare_models(trans$x_formula, data, tree)
    current_aic <- min(model_result$aic_table$AIC)
    
    if(current_aic < best_aic2) {
      best_aic2 <- current_aic
      best_model2 <- model_result
      best_trans2 <- trans
    }
    
    print(paste("Transformation (x:", trans$x_trans, ", y:", trans$y_trans, 
                ") AIC:", current_aic))
  }, error = function(e) {
    print(paste("Failed for transformation (x:", trans$x_trans, ", y:", trans$y_trans, 
                "):", e$message))
  })
}

if(!is.null(best_model2)) {
  model_results2 <- best_model2
  print(paste("Best transformation for Analysis 2: x =", best_trans2$x_trans, 
              ", y =", best_trans2$y_trans))
  print(model_results2$aic_table)
  
  # Update data and tree
  data2 <- model_results2$data
  tree2 <- model_results2$tree
  
  # Get statistics for best model
  summary2 <- summary(model_results2$best_model)
  
  # Calculate R² and adjusted R²
  if(best_trans2$y_trans == "none") {
    r2_values2 <- calculate_adjusted_r2(model_results2$best_model, data2, "CMR")
  } else {
    # For log-transformed CMR, we need to calculate R² differently
    y_pred <- predict(model_results2$best_model)
    y_actual <- log10(data2$CMR)
    ss_res <- sum((y_actual - y_pred)^2)
    ss_tot <- sum((y_actual - mean(y_actual))^2)
    r2_values2 <- list(
      r2 = 1 - (ss_res / ss_tot),
      adj_r2 = 1 - (1 - (1 - (ss_res / ss_tot))) * (nrow(data2) - 1) / (nrow(data2) - 1 - 1)
    )
  }
  
  r2_2 <- r2_values2$r2
  adj_r2_2 <- r2_values2$adj_r2
  p_value2 <- summary2$tTable[2, 4]
} else {
  stop("No valid model found for Analysis 2")
}

# Create plotting function - Nature style
create_plot <- function(data, x_var, y_var, x_lab, y_lab, model_results, r2, adj_r2, p_value, best_model_name, y_trans) {
  # Create base plot - Nature style
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    # Add confidence interval first (so it's behind the points)
    { 
      # Create prediction data
      x_range <- range(data[[x_var]], na.rm = TRUE)
      new_data <- data.frame(x = seq(x_range[1], x_range[2], length.out = 100))
      colnames(new_data) <- x_var
      
      # Calculate confidence intervals using standard linear model
      formula_str <- paste(y_var, "~", x_var)
      lm_model <- lm(as.formula(formula_str), data = data)
      pred <- predict(lm_model, newdata = new_data, interval = "confidence")
      
      # Create prediction data frame
      pred_df <- data.frame(
        x = new_data[[x_var]],
        fit = pred[, "fit"],
        lwr = pred[, "lwr"],
        upr = pred[, "upr"]
      )
      
      # Add confidence interval
      geom_ribbon(data = pred_df, aes(x = x, ymin = lwr, ymax = upr), 
                  fill = "#E1E5EA", alpha = 0.7, inherit.aes = FALSE)
    } +
    # Add data points - Nature style dark blue
    geom_point(color = "#2E5A87", size = 3, alpha = 0.8, shape = 16) +
    # Add fitted line - Nature style blue
    {
      # Create prediction data
      x_range <- range(data[[x_var]], na.rm = TRUE)
      new_data <- data.frame(x = seq(x_range[1], x_range[2], length.out = 100))
      colnames(new_data) <- x_var
      
      # Calculate fitted line using standard linear model
      formula_str <- paste(y_var, "~", x_var)
      lm_model <- lm(as.formula(formula_str), data = data)
      pred <- predict(lm_model, newdata = new_data, interval = "confidence")
      
      # Create prediction data frame
      pred_df <- data.frame(
        x = new_data[[x_var]],
        fit = pred[, "fit"]
      )
      
      # Add fitted line
      geom_line(data = pred_df, aes(x = x, y = fit), color = "#1F77B4", size = 1.5, inherit.aes = FALSE)
    } +
    # Labels and theme
    labs(x = x_lab, y = y_lab) +
    # Nature-style theme
    theme(
      # Clear background and grid
      panel.background = element_rect(fill = "white", colour = "black", size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Axis lines
      axis.line = element_line(colour = "black", size = 0.5),
      
      # Text styling
      axis.title = element_text(size = 14, face = "bold", colour = "black"),
      axis.text = element_text(size = 12, colour = "black"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, colour = "black"),
      
      # No legend
      legend.position = "none"
    )
  
  # Add model information and statistics - Nature style
  model_label <- paste0("Best model: ", best_model_name, "\n",
                       "R² = ", round(r2, 3), ", Adj. R² = ", round(adj_r2, 3), 
                       ", p = ", round(p_value, 4))
  
  # Position annotation in top-left corner
  p <- p + 
    annotate("text", x = -Inf, y = Inf, 
             label = model_label, 
             hjust = -0.05, vjust = 1.5, size = 4.5,
             fontface = "bold", color = "black")
  
  return(p)
}

# Plot 1: CMR vs Adult life expectancy
if(best_trans1$y_trans == "none") {
  y_lab1 <- "CMR"
  y_var1 <- "CMR"
} else {
  y_lab1 <- "CMR (log scale)"
  y_var1 <- "log_CMR"
  data1$log_CMR <- log10(data1$CMR)
}

if(best_trans1$x_trans == "none") {
  x_var1 <- "Adult_life_expectancy_days"
} else {
  x_var1 <- "log_life_expectancy"
  data1$log_life_expectancy <- log10(data1$Adult_life_expectancy_days)
}

p1 <- create_plot(
  data = data1,
  x_var = x_var1,
  y_var = y_var1,
  x_lab = best_trans1$x_lab,
  y_lab = y_lab1,
  model_results = model_results1,
  r2 = r2_1,
  adj_r2 = adj_r2_1,
  p_value = p_value1,
  best_model_name = model_results1$best_model_name,
  y_trans = best_trans1$y_trans
) + 
  ggtitle("PGLS: CMR vs Adult Life Expectancy")

# Plot 2: CMR vs Species body mass
if(best_trans2$y_trans == "none") {
  y_lab2 <- "CMR"
  y_var2 <- "CMR"
} else {
  y_lab2 <- "CMR (log scale)"
  y_var2 <- "log_CMR"
  data2$log_CMR <- log10(data2$CMR)
}

if(best_trans2$x_trans == "none") {
  x_var2 <- "Species_body_mass_kg"
} else {
  x_var2 <- "log_body_mass"
  data2$log_body_mass <- log10(data2$Species_body_mass_kg)
}

p2 <- create_plot(
  data = data2,
  x_var = x_var2,
  y_var = y_var2,
  x_lab = best_trans2$x_lab,
  y_lab = y_lab2,
  model_results = model_results2,
  r2 = r2_2,
  adj_r2 = adj_r2_2,
  p_value = p_value2,
  best_model_name = model_results2$best_model_name,
  y_trans = best_trans2$y_trans
) + 
  ggtitle("PGLS: CMR vs Species Body Mass")

# Display plots
print(p1)
print(p2)

# Save as PDF
ggsave("PGLS_CMR_vs_LifeExpectancy.pdf", p1, width = 8, height = 6)
ggsave("PGLS_CMR_vs_BodyMass.pdf", p2, width = 8, height = 6)

# Output model summaries
print("PGLS Model 1: CMR vs Adult Life Expectancy")
print(summary(model_results1$best_model))
print(paste("Best model:", model_results1$best_model_name))
print(paste("Transformations: x =", best_trans1$x_trans, ", y =", best_trans1$y_trans))
print(paste("R² =", round(r2_1, 4), "Adj. R² =", round(adj_r2_1, 4)))

print("PGLS Model 2: CMR vs Species Body Mass")
print(summary(model_results2$best_model))
print(paste("Best model:", model_results2$best_model_name))
print(paste("Transformations: x =", best_trans2$x_trans, ", y =", best_trans2$y_trans))
print(paste("R² =", round(r2_2, 4), "Adj. R² =", round(adj_r2_2, 4)))

# Output AIC comparison tables
print("Model comparison results - CMR vs Adult Life Expectancy:")
print(model_results1$aic_table)

print("Model comparison results - CMR vs Species Body Mass:")
print(model_results2$aic_table)
