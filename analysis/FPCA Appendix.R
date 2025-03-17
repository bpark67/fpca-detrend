library(tidyverse)

rm(list = ls())

# SIMULATION PARAMETERS
B = 100

# DGP PARAMETERS
m = 20
n = 50
t = seq(0, 1, length.out = n)

beta = 2 # MEAN
sigma = 1 # MEASUREMENT ERROR

# COMPLICATED TIME TREND
#f = function(t){(t)^(1/25) - 1.25*(t)^(37/4)}

# sharpness of change in trend setting

#f = function(t){0.1+(sin(pi*t))^7}

# change point setting
f = function(t){
  par = 0.4
  if(t<=0.45){
    return(t*par/0.45+0.1)
  }
  else if(0.45<t && t<0.55){
    return(par + (1-2*par)/0.1*(t-0.45)+0.1)
  }
  else{
    return(1-par*(1-t)/0.45+0.1)
  }
}

# FILLERS
results = data.frame(
  replicate = rep(1:B, times = 3),
  type = rep(c("simple", "gee", "detrend"), each = B),
  estimate = rep(NA, times = 3*B),
  se = rep(NA, times = 3*B)
)
all_ft = list()

pb = txtProgressBar(min = 0, max = B, style = 3)
for(b in 1:B){
  set.seed(b)
  
  # GENERATE RANDOM EFFECTS
  re = rnorm(m, mean = 0, sd = 1)
  
  
  # GENERATE Y
  Y = matrix(NA, m, n)
  for(i in 1:m){
    Y[i,] = beta + re[i]*sapply(t,f) + rnorm(n, mean = 0, sd = sigma)
  }
  
  
  df = data.frame(
    id = rep(1:m, each = n),
    Y = c(t(Y)),
    X = 1
  )
  
  # STEP 1: SIMPLE REGRESSION
  simple = lm(Y ~ -1 + X,
              data = df)
  
  results[b, 3:4] = summary(simple)$coefficients[1:2]
  
  # STEP 2: GEE
  gee = geepack::geeglm(
    Y ~ -1 + X,
    data = df,
    id = id,
    corstr = "exchangeable",
    std.err = "san.se"
  )
  
  results[(B + b), 3:4] =summary(gee)$coefficients[1:2] %>% as.numeric
  
  # STEP 3: FPCA
  Yt = rbind(Y, t)
  
  fpca = Yt %>% FunOnFun::tibbleFormat(t) %>% FunOnFun::fpcaFormat(id_col = "id")
  
  ## BOOTSTRAP TO ESTIMATE SE
  R = 10
  boot_estimates = rep(NA, R)
  
  pbr = txtProgressBar(min = 0, max = R, style = 2)
  
  for(r in 1:R){
    boot_indices = sample(1:m, m, replace = T)
    boot_Y = Y[boot_indices,]
    
    boot_df = data.frame(
      id = rep(1:m, each = n),
      Y = c(t(boot_Y)),
      X = 1
    )
    
    boot_Yt = rbind(Y, t)
    
    boot_fpca = boot_Yt %>% FunOnFun::tibbleFormat(t) %>% FunOnFun::fpcaFormat(id_col = "id")
    
    boot_res = fdapace::FPCA(
      boot_fpca$Variable1,
      boot_fpca$Time,
      list(dataType = "Sparse",
           error = T,
           kernel = "epan",
           verbose = F,
           nRegGrid = length(t))
    )
    
    boot_trend = t(boot_res$phi[,1, drop = F] %*% t(boot_res$xiEst[, 1, drop = F]))
    boot_detrend = boot_Y - boot_trend
    
    boot_df_detrend = data.frame(
      id = rep(1:m, each = n),
      Y = c(t(boot_detrend)),
      X = 1
    )
    
    boot_gee_detrend = geepack::geeglm(
      Y ~ -1 + X,
      data = boot_df_detrend,
      id = id,
      corstr = "exchangeable",
      std.err = "san.se"
    )
    
   boot_estimates[r] = boot_gee_detrend$coefficients %>% as.numeric
    setTxtProgressBar(pbr, r)
  }
  
  close(pbr)
  

  res = fdapace::FPCA(fpca$Variable1,
                      fpca$Time,
                      list(dataType = "Sparse",
                           error = T,
                           kernel = "epan",
                           verbose = F,
                           nRegGrid = length(t)))
  
  all_ft[[b]] = res$phi[,1]
  trend = t(res$phi[,1, drop = F] %*% t(res$xiEst[, 1, drop = F]))
  detrend = Y[1:m,] - trend
  
  df_detrend = data.frame(
    id = rep(1:m, each = n),
    Y = c(t(detrend)),
    X = 1
  )
  
  gee_detrend = geepack::geeglm(
    Y ~ -1 + X,
    data = df_detrend,
    id = id,
    corstr = "exchangeable",
    std.err = "san.se"
  )
  
  results[(2*B + b),3] = summary(gee_detrend)$coefficients[1] %>% as.numeric
  results[(2*B + b),4] = sd(boot_estimates)
  
  setTxtProgressBar(pb, b)
}

close(pb)

####### VISUALIZE

results$ub = results$estimate + 2*results$se
results$lb = results$estimate - 2*results$se

results$covered = (2 < results$ub) & (2 > results$lb)
results$bias = results$estimate - 2

# load("../results/results.RData")

results %>%
  ggplot(aes(x = replicate, color = covered)) +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0) +
  facet_wrap(~type, ncol = 1, labeller = labeller(
    type = c(
    "detrend" = "Proposed Approach",
    "gee" = "GEE",
    "simple" = "Naive Model"
  ))) +
  scale_color_manual(
    values = c("FALSE" = "#D55E00",  # Dark red
               "TRUE" = "#0072B2"),   # Navy blue
    labels = c("FALSE" = "Not Covered", "TRUE" = "Covered")
  ) +
  theme_bw() +
  labs(x = "Iteration",
       color = "Coverage") +
  theme(    
    strip.background = element_rect(fill = "skyblue", color = "black"),
    strip.text = element_text(face = "bold")
    )

results %>%
  ggplot(aes(x = type, y = abs(bias), fill = type)) +
  geom_boxplot() +
  theme_bw()

# save(results, file = "../results/results.RData")


trend = t(res$phi[,1, drop =F] %*% t(res$xiEst[, 1, drop = F]))  

df_trend = data.frame(
  time = rep(t, 101),
  Type = c(rep("Truth", each = length(t)), rep(rep("Estimated", each = length(t)), 100)),
  Iteration = rep(0:100, each = length(t)),
  f_t = c(sapply(t,f), unlist(all_ft))
)

# save(df_trend, file = "../results/df_trend.RData")

df_summary = df_trend %>%
  filter(Type == "Estimated") %>%
  group_by(time) %>%
  summarize(min_f = min(f_t),
            max_f = max(f_t), .group = "drop")

ggplot()+
  geom_ribbon(data = df_summary, aes(x = time, ymin = min_f, ymax = max_f),
              fill = "#a6cee3") +
  geom_line(data = df_trend %>% filter(Type == "Truth"), aes(x = time, y = f_t),
            color = "#1b4f72", linewidth = 1) +
  annotate("text", x = 0.5, y = 1, label = "Truth", 
           color = "#1b4f72", fontface = "bold", vjust = -1) +
  annotate("text", x = 0.85, y = 1.3, label = "Estimation Range", 
           color = "#a6cee3", fontface = "bold", vjust = 1) +
  theme_bw() +
  labs(x = "Time", y = expression(sapply(t,f))) +
  theme(legend.position = "none")

df_trend %>%
  ggplot(aes(x = time, y = f_t, group = Type, color = Type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Truth" = "black", "Estimated" = "red")) +
  theme_bw() + 
  labs(
    y = expression(sapply(t,f)),
    x = "Time"
  )


### COVARIANCE PLOT

# # Create grid coordinates
# x <- seq(1, 50, length.out = 50)
# y <- seq(1, 50, length.out = 50)
# 
# # Create a 3D surface plot
# fig <- plotly::plot_ly(
#   x = x, y = y, z = cov(Y), 
#   type = "surface",
#   colorscale = "Viridis"
# )
# 
# # Customize layout for publication quality
# fig <- fig %>%
#   plotly::layout(
#     title = "3D Visualization of Covariance Matrix",
#     scene = list(
#       xaxis = list(title = "Row Index"),
#       yaxis = list(title = "Column Index"),
#       zaxis = list(title = "Covariance"),
#       camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))  # Adjust viewing angle
#     )
#   )
# 
# # Show the plot
# fig
cv_mat = cov(Y)
diag(cv_mat) = NA
cv = reshape2::melt(cv_mat)

# Create the heatmap plot
ggplot(cv, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", mid = "lightgray", high = "red", midpoint = 0) +  # Use a colorblind-friendly scale
  labs(fill = "Covariance",
       y = element_blank(),
       x = element_blank()) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),  # Hide axis text if too cluttered
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
