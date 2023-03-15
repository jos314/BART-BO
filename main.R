library(tidyverse)
library(tidymodels)
library(dplyr)
library(ggplot2)
library(WSGeometry) # Create 2d grid
library(phonTools)  # zeros function
library(purrr)      # map
library(ramify)     # argmax, linspace

library(tidybayes)  # Revisar esta librería


#Rosenbrock function
rosen <- function(x) {
  return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
}

#Otra fun
gradient <- function(x) {
  return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]), 200 * (x[2] - x[1] * x[1]) ) )
}

#IN: Data: D = (y,x); y in R^k, x in R^(kxh)
#OUT: Candidate model of function f

#Begin with s = 1
#While iter < maxiter and ..
#3.- Train the probabilistic model based on D
#4.- Calculate the selected aquisition function u(x|D)
#5.- Choose the next experiment point x(s+1)
#6.- Evaluate x(s+1) on f and add it D
#7.- s++
#8.- Return candidate model of the function f


################  Sin  ######################

# Este script hace optimización bayesiana sobre la función sin() como función de
# caja negra, tomando el modelo bart de tidymodels como método para hacer 
# aproximaciones de la función. 

set.seed(1)
dimFx <- 1 # Dimensión de la entrada de la función 
dimFy <- 1 # Dimensión de la salida de la función
s <- 1 
maxiter <- 5 # máximo de iteraciones a correr con el modelo 
num_point_pred <- 200  #Num de puntos a predecir en cada x_est
init_points <- 8 #Init evals of function
n <- 100 #Tamaño del grid
minx <- -pi
maxx <- pi

# Datos iniciales para evaluar la función
x <- as.data.frame(ramify::linspace(a = minx, b = maxx, n = init_points)) %>% 
  rename(x1 = 1)
y <- sin(x) %>% rename(y1 = 1)
D <- cbind(x,y)

# Calculamos los datos iniciales 
D %>% ggplot(aes(x = x1, y = y1)) +
  #geom_line() +
  geom_point() +
  ggtitle("Initial function evaluations: sin(x)") +
  xlab("x") + ylab("sin(x)")

# Calculamos x que minimice f de los datos iniciales
x_est <- D[argmin(D['y1'], rows = FALSE),1] 


# Definimos la funion de utilidad
acq_function <- function(f_est_hat, grid_x) {
  # Compute the mean and standard deviation of the predicted values
  #x_est_hat minimices f
  #f_est_hat <- f_hat(x_est_hat) bart evaluated on the minimicer of f
  #f_est_hat 
  u <- as.data.frame(zeros(nrow = nrow(grid_x), ncol = 2))
  u[,1] <- grid_x[,1]
  u[,2] <- max(0, f_est_hat - grid_x[, 'lci'] )
  colnames(u) <-  c('x1', 'utility_exp')
  return(u)
}

# Data Frame para guardar los resultados obtenidos en cada iteración (D + x_s)
resultados <- tibble(x1 = numeric(), y = numeric())


# Mientras algunas restricciones se cumplan, seguir buscando el punto que optimice f 
while (s < maxiter){   # mejora < eps...

  # Definiendo el modelo BART
  sin_bart <- bart()%>% 
    set_engine("dbarts") %>% 
    set_mode("regression") 
  
  sin_bart_fit <- sin_bart %>% 
    fit_xy(x = x, y = y)
  
  # Se genera el grid en 1d
  x_grid <- as.data.frame(ramify::linspace(a = minx, b = maxx, n = n)) %>% 
    cbind(zeros(n,4)) 
  names(x_grid) <-  c('x1', 'interval_max_incert', 'uci', 'lci', 'mu')
  
  
  #pred <- predict(sin_bart_fit, new_data = x_grid[i, , drop = FALSE])
  # **Revisar los métodos para simular
  for (i in 1:nrow(x_grid)){
    pred_at_x <- map_dfr(1:num_point_pred, ~predict(sin_bart_fit, new_data = x_grid[i,1, drop = FALSE]), type = "prob") %>% 
      as.data.frame()
    mu <- pred_at_x$fit
    sigma <- pred_at_x$sd
    
    lci <- quantile(pred_at_x$.pred, probs = 0.025)
    mu <- mean(pred_at_x$.pred)
    uci <- quantile(pred_at_x$.pred, probs = 0.975)
    x_grid[i, 2:5] <- c(abs(uci-lci), uci, lci, mu)
  }
  
  # Identificar f_est_hat 
  f_est_hat <- D$y1 %>% min()
  
  # Evaluamos la función de utilidad
  u <- acq_function(f_est_hat, x_grid)
  
  # Añadir un punto a los datos D
  x_est_prov <- u
  
  #Actualizar x_est
  if (x_est > x_est_prov){
    x_est <- x_est_prov
  }
  
  #Graficar el paso 
  p <- D %>% ggplot(aes(x = x1, y = y1)) +
    #geom_line() +
    geom_point() +
    geom_line(data=x_grid, aes(x=x1, y=uci), color='red') +
    geom_line(data=x_grid, aes(x=x1, y=lci), color='red') +
    ggtitle(paste0("BO step ", s)) +
    xlab("x") + ylab("sin(x)")
  
  print(p)
  
  #Graficar la funcion de utilidad
  p <- u %>% ggplot(aes(x1, utility_exp)) +
    #geom_point() +
    geom_line() +
    ggtitle(paste0("Acquisition function in step ", s)) +
    ylab("Expected value of u") + xlab("x")
  
  print(p)
    
  
  s <- s + 1
}


#while (s < maxiter){   # mejora < eps...
  
  # Fitting BART model
  
  sin_bart <- bart() %>% 
    set_engine("dbarts") %>% 
    set_mode("regression") 
  
  sin_bart_fit <- sin_bart %>% 
    fit_xy(x = x, y = y)
  
  # Grid 1d
  x_grid <- ramify::linspace(a = minx, b = maxx, n = n) %>% 
    as.data.frame() %>% 
    cbind(zeros(n,4)) 
  names(x_grid) <-  c('x1', 'expected_ux', 'uci', 'lci')
  
  #pred <- predict(sin_bart_fit, new_data = x_grid[i, , drop = FALSE])
  
  for (i in 1:nrow(x_grid)){
    pred_at_x <- map_dfr(1:num_point_pred, ~predict(sin_bart_fit, new_data = x_grid[i,1, drop = FALSE])) %>% 
      as.data.frame()
    lci <- quantile(pred_at_x$.pred, probs = 0.025)
    uci <- quantile(pred_at_x$.pred, probs = 0.975)
    x_grid[i,3] <- uci
    x_grid[i,4] <- lci
    x_grid[i,2] <- max(0,x_est-lci)
  }
  x_grid$expected_ex <- max(0, x_est - x_grid$lci)
  x_est_prov <- x_grid[argmax(x_grid['interval_incert'], rows = FALSE),1] 
  xy <- rbind(xy, c(x_est, sin(x_est)))
  
  #Graficar el paso ... .. .

  s <- s + 1
#}


#format(round(xy, 3), nsmall = 2)

# Graphing the variance and update of each iteration

xy %>% ggplot(aes(x = x1, y = y1)) +
  geom_line()

