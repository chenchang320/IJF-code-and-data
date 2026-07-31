library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)


hs_300 <- read_xlsx('/Users/mac/Desktop/GitHub-Englishi version/data/HS300 dataset.xlsx')


hs_300_dates <- as.Date(hs_300[[1]])
hs_300_returns <- as.numeric(hs_300[[2]]) / 100


Date <- hs_300_dates[1:974]


ES = as.matrix(read.csv('/Users/mac/Desktop/GitHub-Englishi version/data/simulated CVaR with vine.csv'))

CVaR_99 <- ES[, c(1, 4, 5, 6)] * (-1)
CVaR_95 <- ES[, c(8, 10, 11, 12)] * (-1)
CVaR_90 <- ES[, c(15, 16, 17, 18)] * (-1)


strategy_names_0.01 <- c("SCAD_0.01", "EWP", "MVP", "MVO")
strategy_names_0.05 <- c("SCAD_0.05", "EWP", "MVP", "MVO")
strategy_names_0.1 <- c("SCAD_0.1", "EWP", "MVP", "MVO")


df_99 <- data.frame(Date = Date, CVaR_99)
colnames(df_99)[-1] <- strategy_names_0.01

df_95 <- data.frame(Date = Date, CVaR_95)
colnames(df_95)[-1] <- strategy_names_0.05

df_90 <- data.frame(Date = Date, CVaR_90)
colnames(df_90)[-1] <- strategy_names_0.1


df_99_long <- melt(df_99, id.vars = "Date") %>% mutate(level = "99%")
df_95_long <- melt(df_95, id.vars = "Date") %>% mutate(level = "95%")
df_90_long <- melt(df_90, id.vars = "Date") %>% mutate(level = "90%")


p1 <- ggplot(df_99_long, aes(x = Date, y = value, color = variable)) +
  geom_line() +
  labs(x = "Date", y = "CVaR") +
  theme_minimal() +
  theme(legend.title = element_blank())

p2 <- ggplot(df_95_long, aes(x = Date, y = value, color = variable)) +
  geom_line() +
  labs(x = "Date", y = "CVaR") +
  theme_minimal() +
  theme(legend.title = element_blank())

p3 <- ggplot(df_90_long, aes(x = Date, y = value, color = variable)) +
  geom_line() +
  labs(x = "Date", y = "CVaR") +
  theme_minimal() +
  theme(legend.title = element_blank())


print(p1)
print(p2)
print(p3)






'================================'
library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)


hs_300 <- read_xlsx('/Users/mac/Desktop/GitHub-Englishi version/data/HS300 dataset.xlsx')
hs_300_dates <- as.Date(hs_300[[1]])
hs_300_returns <- as.numeric(hs_300[[2]]) / 100


Date <- hs_300_dates[1:974]


VaR = as.matrix(read.csv('/Users/mac/Desktop/GitHub-Englishi version/data/simulated VaR with vine.csv'))

VaR_99 <- VaR[, c(1, 4, 5, 6)] * (-1)
VaR_95 <- VaR[, c(8, 10, 11, 12)] * (-1)
VaR_90 <- VaR[, c(15, 16, 17, 18)] * (-1)


strategy_names_0.01 <- c("SCAD_0.01", "EWP", "MVP", "MVO")
strategy_names_0.05 <- c("SCAD_0.05", "EWP", "MVP", "MVO")
strategy_names_0.1  <- c("SCAD_0.1",  "EWP", "MVP", "MVO")


df_99_var <- data.frame(Date = Date, VaR_99)
colnames(df_99_var)[-1] <- strategy_names_0.01

df_95_var <- data.frame(Date = Date, VaR_95)
colnames(df_95_var)[-1] <- strategy_names_0.05

df_90_var <- data.frame(Date = Date, VaR_90)
colnames(df_90_var)[-1] <- strategy_names_0.1


df_99_long_var <- melt(df_99_var, id.vars = "Date") %>% mutate(level = "99%")
df_95_long_var <- melt(df_95_var, id.vars = "Date") %>% mutate(level = "95%")
df_90_long_var <- melt(df_90_var, id.vars = "Date") %>% mutate(level = "90%")


p1_var <- ggplot(df_99_long_var, aes(x = Date, y = value, color = variable)) +
  geom_line() +
  labs(x = "Date", y = "VaR") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p2_var <- ggplot(df_95_long_var, aes(x = Date, y = value, color = variable)) +
  geom_line() +
  labs(x = "Date", y = "VaR") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p3_var <- ggplot(df_90_long_var, aes(x = Date, y = value, color = variable)) +
  geom_line() +
  labs(x = "Date", y = "VaR") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


print(p1_var)
print(p2_var)
print(p3_var)





'================================'
library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)


hs_300 <- read_xlsx('/Users/mac/Desktop/GitHub-Englishi version/data/HS300 dataset.xlsx')
hs_300_dates <- as.Date(hs_300[[1]])
hs_300_returns <- as.numeric(hs_300[[2]]) / 100
Date <- hs_300_dates[1:974]


VAR = as.matrix(read.csv('/Users/mac/Desktop/GitHub-Englishi version/data/simulated VaR without vine.csv'))

with_copula_0.01 <- cbind(VaR[, 1], VAR[, 1]) * (-1)
with_copula_0.05 <- cbind(VaR[, 8], VAR[, 8]) * (-1)
with_copula_0.1  <- cbind(VaR[, 15], VAR[, 15]) * (-1)


strategy_names_0.01 <- c("SCAD_0.01", "SCAD_0.01_without_Copula")
strategy_names_0.05 <- c("SCAD_0.05", "SCAD_0.05_without_Copula")
strategy_names_0.1  <- c("SCAD_0.1",  "SCAD_0.1_without_Copula")


df_01 <- data.frame(Date = Date, with_copula_0.01)
colnames(df_01)[-1] <- strategy_names_0.01

df_05 <- data.frame(Date = Date, with_copula_0.05)
colnames(df_05)[-1] <- strategy_names_0.05

df_10 <- data.frame(Date = Date, with_copula_0.1)
colnames(df_10)[-1] <- strategy_names_0.1


df_01_long <- melt(df_01, id.vars = "Date")
df_05_long <- melt(df_05, id.vars = "Date")
df_10_long <- melt(df_10, id.vars = "Date")


plot_var_compare <- function(df_long, level_text) {
  ggplot(df_long, aes(x = Date, y = value, color = variable)) +
    geom_line(size = 1.1) +
    labs(x = "Date", y = "VaR") +
    theme_minimal(base_size = 13) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


p01 <- plot_var_compare(df_01_long, "VaR at 99%")
p05 <- plot_var_compare(df_05_long, "VaR at 95%")
p10 <- plot_var_compare(df_10_long, "VaR at 90%")


print(p01)
print(p05)
print(p10)



'================================'
library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)


hs_300 <- read_xlsx('/Users/mac/Desktop/GitHub-Englishi version/data/HS300 dataset.xlsx')
hs_300_dates <- as.Date(hs_300[[1]])
Date <- hs_300_dates[1:974]


CVAR = as.matrix(read.csv('/Users/mac/Desktop/GitHub-Englishi version/data/simulated CVaR without vine.csv'))

with_copula_0.01 <- cbind(ES[, 1], CVAR[, 1]) * (-1)
with_copula_0.05 <- cbind(ES[, 8], CVAR[, 8]) * (-1)
with_copula_0.1  <- cbind(ES[, 15], CVAR[, 15]) * (-1)


strategy_names_0.01 <- c("SCAD_0.01", "SCAD_0.01_without_Copula")
strategy_names_0.05 <- c("SCAD_0.05", "SCAD_0.05_without_Copula")
strategy_names_0.1  <- c("SCAD_0.1",  "SCAD_0.1_without_Copula")


df_01 <- data.frame(Date = Date, with_copula_0.01)
colnames(df_01)[-1] <- strategy_names_0.01

df_05 <- data.frame(Date = Date, with_copula_0.05)
colnames(df_05)[-1] <- strategy_names_0.05

df_10 <- data.frame(Date = Date, with_copula_0.1)
colnames(df_10)[-1] <- strategy_names_0.1


df_01_long <- melt(df_01, id.vars = "Date")
df_05_long <- melt(df_05, id.vars = "Date")
df_10_long <- melt(df_10, id.vars = "Date")


plot_cvar_compare <- function(df_long) {
  ggplot(df_long, aes(x = Date, y = value, color = variable)) +
    geom_line(size = 1.1) +
    labs(x = "Date", y = "CVaR") +
    theme_minimal(base_size = 13) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


p01 <- plot_cvar_compare(df_01_long)
p05 <- plot_cvar_compare(df_05_long)
p10 <- plot_cvar_compare(df_10_long)


print(p01)
print(p05)
print(p10)




