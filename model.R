# Konstruujemy drzewo niezależnie od ścieżki (PATH-INDEPENDENT) jako macierz, 
# gdzie t-ta kolumna to S_t moment wykonania akcji

# Dane
d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))    # d = 1 / u = u**(-1)
S_0 <- 50
r <- 0.02
K <- 48
T <- 2

matrix_size <- function(d_t, T){    # rozmiar macierzy kwadratowej N, gdzie uwzględniamy moment S_0
  return(T / d_t + 1)
}

N <- matrix_size(d_t, T)

binomial_tree <- function(S_0, u, d_t, T){
  
  N <- matrix_size(d_t, T)    # liczba kroków w drzewie, rozmiar macierzy
  A <- matrix(NA, nrow = N, ncol = N)
  
  for (i in 1:N){
    A[N:(N - i + 1), i] <- u**(seq(-(i - 1), (i - 1), 2))
  }

  return(S_0 * A)
}

S_T <- binomial_tree(S_0, u, d_t, T)

# Path-independent tree
y <- S_T[, 1:25]
x <- seq(1, 25, length = length(y))
colors <- ifelse(S_T > K, "#DEA435", "#40DE8E")
plot(x, y, log = "y", type = "p", col = colors, pch = 16,
     xlab = "Moments", ylab = "", main = "Path-independent tree", yaxt = "n")
legend("topleft", legend = c("S_T > K", "S_T < K"), col = c("#DEA435", "#40DE8E"), pch = 16, cex = 1.1, bty = "n")

wycena <- function(u, d, r, d_t, V_u, V_d){
  p <- (exp(r * d_t) - d) / (u - d)
  return(exp(-r * d_t) * (p * V_u + (1 - p) * V_d))
}
wycena <- Vectorize(wycena, c("V_u", "V_d"))

european_option <- function(S_0, u, d, r, K, d_t, T, type = "put"){
  
  N <- matrix_size(d_t, T)    # liczba kroków w drzewie, rozmiar macierzy
  S_T <- binomial_tree(S_0, u, d_t, T)
  B <- matrix(0, nrow = N, ncol = N)    # macierz payoff
  B[is.na(S_T)] <- NA
  ifelse(type == "put",
         B[, N] <- pmax(K - S_T[, N], 0),    # opcja put
         B[, N] <- pmax(S_T[, N] - K, 0))    # opcja call

  for (i in (N - 1):1){
    B[(N - i + 1):N, i] <- wycena(u, d, r, d_t, B[(N - i):(N - 1), i + 1], B[(N - i + 1):N, i + 1])
  }

  return(B)
}

american_option <- function(S_0, u, d, r, K, d_t, T, type = "put"){
  
  N <- matrix_size(d_t, T)    # liczba kroków w drzewie, rozmiar macierzy
  S_T <- binomial_tree(S_0, u, d_t, T)
  B <- matrix(0, nrow = N, ncol = N)    # macierz payoff
  B[is.na(S_T)] <- NA
  if(type == "put"){
    B[, N] <- pmax(K - S_T[, N], 0)
    for (i in (N - 1):1){
      for (j in (N - i + 1):N){
        a <- wycena(u, d, r, d_t, B[j - 1, i + 1], B[j, i + 1])
        b <- max(K - S_T[j, i], 0)
        B[j, i] <- max(a, b)
      }
    }
  }
  if(type == "call"){
    B[, N] <- pmax(S_T[, N] - K, 0)
    for (i in (N - 1):1){
      for (j in (N - i + 1):N){
        a <- wycena(u, d, r, d_t, B[j - 1, i + 1], B[j, i + 1])
        b <- max(S_T[j, i] - K, 0)
        B[j, i] <- max(a, b)
      }
    }
  }
  return(B)
}

EU_put <- round(european_option(S_0, u, d, r, K, d_t, T, type = "put")[N, 1], 2)
EU_call <- round(european_option(S_0, u, d, r, K, d_t, T, type = "call")[N, 1], 2)
AM_put <- round(american_option(S_0, u, d, r, K, d_t, T, type = "put")[N, 1], 2)
AM_call <- round(american_option(S_0, u, d, r, K, d_t, T, type = "call")[N, 1], 2)

#--------------------------------

# Momenty wykonania opcji amerykańskich put
AM_put <- american_option(S_0, u, d, r, K, d_t, T, type = "put")
y <- S_T[, 1:25]
x <- seq(1, 25, length = length(y))
colors <- ifelse(AM_put == K - S_T, "#3DC435", "#DE341F")
plot(x, y, log = "y", type = "p", col = colors, pch = 16,
     xlab = "Moments", ylab = "", main = "Momenty wykonania amerykańskich opcji put", yaxt = "n")

# Momenty wykonania opcji amerykańskich call
AM_call <- american_option(S_0, u, d, r, K, d_t, T, type = "call")
y <- S_T[, 1:25]
x <- seq(1, 25, length = length(y))
colors <- ifelse(AM_call == S_T - K, "#3DC435", "#DE341F")
plot(x, y, log = "y", type = "p", col = colors, pch = 16,
     xlab = "Moments", ylab = "", main = "Momenty wykonania amerykańskich opcji call", yaxt = "n")



execution_time_put <-  which(AM_put == K - S_T, arr.ind = TRUE)
colnames(execution_time_put) <- c("wariant", "moment")
execution_time_put <- execution_time_put[, 2:1]

# Momenty wykonania opcji amerykańskich call
execution_time_call <-  which(american_option(S_0, u, d, r, K, d_t, T, type = 'call') == S_T - K, arr.ind = TRUE)
colnames(execution_time_call) <- c("wariant", "moment")
execution_time_call <- execution_time_call[, 2:1]
execution_time_put <- data.frame(wariant = execution_time_put, moment = 1:length(execution_time_put))
execution_time_put <- execution_time_put[order(execution_time_put$moment), ]
execution_time_put
