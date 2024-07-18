library(seminr)
library(tidyverse)
library(logging)
library(MASS)
seminr::mobi

# set up model for testing purposes
measurements <- constructs(
  composite("Image",        multi_items("IMAG", 1:5)),
  composite("Expectation",  multi_items("CUEX", 1:3)),
  composite("Value",        multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3))
)
structure <- relationships(
  paths(from = c("Image", "Expectation"), to = "Value"),
  paths(from = "Value", to = "Satisfaction")
)

# Estimate using PLS-PM from model parts defined earlier  
pls_model <- estimate_pls(data = mobi, 
                          measurement_model = measurements, 
                          structural_model = structure)
summary(pls_model)
plot(pls_model)
pls_model$path_coef
model <- specify_model(measurements, structure)


# actual functions for package


# get all latent variables from model
getLatentVariables <- function(model){
  sm <- model$structural_model
  # get all strings vom source and target column of sm
  src <- sm[,"source"]
  tgt <- sm[,"target"]
  all <- c(src, tgt) |> unique() 
  # remove all interaction term (e.g. Test*Trial)
  all <- all[!grepl("\\*", all)]
  all
}

# get an empty path matrix from model
getEmptyPathMatrix <- function(model){
  latent_variables <- getLatentVariables(model)
  path_matrix <- matrix(0, nrow = length(latent_variables), ncol = length(latent_variables)) 
  rownames(path_matrix) <- latent_variables
  colnames(path_matrix) <- latent_variables
  
  sm <- model$structural_model
  logdebug(sm)
  for(i in 1:nrow(sm)){
    # debug
    #i <- 1
    src <- sm[i, "source"]
    tgt <- sm[i, "target"]
    path_matrix[src, tgt] <- 1
  }
  
  path_matrix
}

# set a path coefficient in path matrix
setPathCoef <- function(path_matrix, from, to, value){
  path_matrix[from, to] <- value
  path_matrix
}


# convert path_matrix to covariance matrix sigma
path_matrix_to_sigma <- function(path_matrix){
  # take square roots of all entries
  #path_matrix <- path_matrix |> sqrt()
  
  for (i in 1:nrow(path_matrix)){
    path_matrix[i,i] <- 1
  }
  # copy upper triangle to lower triangle
  for (i in 1:nrow(path_matrix)){
    for (j in 1:ncol(path_matrix)){
      if (i > j){
        path_matrix[i,j] <- path_matrix[j,i]
      }
    }
  }
  path_matrix
}

# test if matrix is positive definite
ispositivedefinite<-function(x, tol=1e-8) {
  all(eigen(x, only.values=TRUE)$values > tol)
}

# simulate structural data from path matrix 
generateStructuralData <- function(path_matrix, n = 1000){
  
  sigma <- path_matrix_to_sigma(path_matrix)  
  raw_data <- mvrnorm(n, mu = rep(0, nrow(sigma)), Sigma = sigma, empirical = TRUE)
  
  res_data <- as.data.frame(raw_data)
  names(res_data) <- rownames(path_matrix)
  
  res_data
}

# test the code

getLatentVariables(model)
path_matrix <- getEmptyPathMatrix(model)

path_matrix <- path_matrix |> 
  setPathCoef("Image", "Value", runif(1,0.15,0.8)) |> 
  setPathCoef("Expectation", "Value", runif(1,0.15,0.8)) |>
  setPathCoef("Value", "Satisfaction", runif(1,0.15,0.8))

# generate structural data
sd_data <- path_matrix |> 
  generateStructuralData()

eigen(path_matrix_to_sigma(path_matrix), only.values=TRUE)

ispositivedefinite(path_matrix_to_sigma(path_matrix))
path_matrix_to_sigma(path_matrix)


# test linear model frrom image to expectation
lm(data = as.data.frame(sd_data), Image ~ Value)
lm(data = as.data.frame(sd_data), Expectation ~ Value)  
lm(data = as.data.frame(sd_data), Value ~ Satisfaction) 


complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

complement(sd_data$Value, 0.45)

# generate indicator data
generateMeasurementData <- function(model, sd_data, cronbach_alpha = 0.8, cronbach_alpha_sd = 0.1){
  measurements <- model$measurement_model
  res_data <- data.frame(id = 1:nrow(sd_data))
  for (i in 1:length(measurements)){
    # debug
    #i <- 1
    construct <- measurements[i]
    items <- construct$composite
    for (j in 1:(length(items)/3)){
      # debug
      #j <- 2
      # get the item name as the second element of the list
      construct_name <- items[(j-1)*3 + 1]
      itemname <- items[(j-1)*3 + 2]
      itemtype <- items[(j-1)*3 + 3]
      coef <- rnorm(1, cronbach_alpha, cronbach_alpha_sd)
      # generate data for item
      column <- as.vector(complement(sd_data[[construct_name]], coef) )
     
      column_frame <- data.frame(column)
      names(column_frame)[1] <- itemname
      
      res_data <- res_data |> cbind(column_frame)
     
    }
  }
  res_data
}

mm_data <- tibble( generateMeasurementData(model, sd_data) )

new_model <- estimate_pls(mm_data, model = model)


plot(new_model)
plot(pls_model)

model$measurement_model[1]$composite[1]

