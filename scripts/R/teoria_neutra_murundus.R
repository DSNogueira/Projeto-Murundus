# =============================================================================
# teoria_neutra_murundus.R
# Projeto Murundus — Teoria Neutra da Biodiversidade e Biogeografia (UNTB)
# =============================================================================
# Autores: Denis S. Nogueira
# Descrição: Estimativa dos parâmetros neutros (θ e m) para comunidades de
#   murundus, ajuste de SADs neutras, e testes de neutralidade.
#
# Referências principais:
#   Hubbell (2001). The Unified Neutral Theory. PUP.
#   Etienne (2005). Ecology Letters 8:253-260.
#   McGill et al. (2006). Ecology 87:1411-1423.
# =============================================================================

# --- Dependências ---
library(ggplot2)
library(dplyr)

# Carregar pacote untb se disponível
if (requireNamespace("untb", quietly = TRUE)) {
  library(untb)
  untb_disponivel <- TRUE
} else {
  message("Pacote 'untb' não instalado. Usando implementações internas.")
  untb_disponivel <- FALSE
}

script_dir <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
                   rstudioapi::isAvailable()) {
  dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  dirname(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE))
}
source(file.path(script_dir, "funcoes_ecologicas.R"))

# =============================================================================
# 1. DADOS
# =============================================================================
# Substitua por dados reais de abundância por murundu:
#   dados_abund <- read.csv("../../dados/abundancia_murundus.csv")
# Simulação de comunidade neutra para demonstração:
set.seed(2024)
J_local <- 200   # indivíduos na comunidade local
theta_verdadeiro <- 8
m_verdadeiro     <- 0.15
comunidade_sim   <- simular_comunidade_neutra(
  J     = J_local,
  theta = theta_verdadeiro,
  m     = m_verdadeiro,
  t_max = J_local * 50,
  semente = 42
)
cat(sprintf("Comunidade simulada: %d espécies, %d indivíduos\n",
            nrow(comunidade_sim), sum(comunidade_sim$abundancia)))

# =============================================================================
# 2. DISTRIBUIÇÃO DE ABUNDÂNCIAS DE ESPÉCIES (SAD)
# =============================================================================

# 2.1 SAD observada: histograma de Preston (octaves)
abundancias_obs <- comunidade_sim$abundancia
log2_classes    <- floor(log2(abundancias_obs)) + 1
sad_obs <- data.frame(
  octave    = seq_len(max(log2_classes)),
  n_especies = tabulate(log2_classes, nbins = max(log2_classes))
)
sad_obs$abund_min <- 2^(sad_obs$octave - 1)
sad_obs$abund_max <- 2^sad_obs$octave - 1

p_sad <- ggplot(sad_obs, aes(x = factor(octave), y = n_especies)) +
  geom_col(fill = "#2C7A3F", color = "white") +
  scale_x_discrete(
    labels = paste0(sad_obs$abund_min, "–", sad_obs$abund_max)
  ) +
  labs(
    title    = "Distribuição de Abundâncias (Octaves de Preston)",
    subtitle = sprintf("Comunidade simulada: θ = %g, m = %g", theta_verdadeiro, m_verdadeiro),
    x        = "Classe de abundância",
    y        = "Número de espécies"
  ) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p_sad)

# 2.2 Rank-abundance (Whittaker) plot
rank_abd <- sort(abundancias_obs, decreasing = TRUE)
p_rank <- ggplot(
  data.frame(rank = seq_along(rank_abd), abundancia = rank_abd),
  aes(x = rank, y = abundancia)
) +
  geom_line(color = "#2C7A3F", linewidth = 1) +
  geom_point(size = 1.5, color = "#2C7A3F", alpha = 0.7) +
  scale_y_log10() +
  labs(
    title = "Rank-Abundance (Diagrama de Whittaker)",
    x     = "Rank de abundância",
    y     = "Abundância [log10]"
  ) +
  theme_classic(base_size = 13)

print(p_rank)

# =============================================================================
# 3. ESTIMATIVA DE THETA POR MÁXIMA VEROSSIMILHANÇA (Ewens)
# =============================================================================
est_theta <- estimar_theta(abundancias_obs)
cat(sprintf("\nEstimativa de θ (Ewens MLE): %.3f\n", est_theta$theta))
cat(sprintf("Valor verdadeiro de θ = %g\n", theta_verdadeiro))
cat(sprintf("Log-verossimilhança: %.3f\n", est_theta$loglik))

# =============================================================================
# 4. ESTIMATIVA DE m (TAXA DE IMIGRAÇÃO) — Etienne (2005)
# =============================================================================
# A taxa m pode ser estimada via ESF (Etienne Sampling Formula) com
# o pacote 'untb'. Aqui mostramos a abordagem via gradiente de log-lik
# num grid de (theta, m):

theta_grid <- seq(2, 20, by = 1)
m_grid     <- seq(0.01, 0.5, by = 0.02)
J          <- sum(abundancias_obs)
S_obs      <- length(abundancias_obs)

# Aproximação: verossimilhança via riqueza esperada
# E[S | theta, m, J] = theta * sum_{k=1}^{J} I/(I + k - 1)
# onde I = m*(J-1)/(1-m)
riqueza_esperada_neutra <- function(theta, m, J) {
  I  <- m * (J - 1) / (1 - m)
  soma <- sum(I / (I + seq(0, J - 1)))
  theta * soma / I  # normalização aproximada
}

# Grade 2D de log-verossimilhança (simplificada, via Ewens)
loglik_grid <- outer(theta_grid, m_grid, function(th, m) {
  mapply(function(t, mi) {
    tryCatch(loglik_ewens(abundancias_obs, t), error = function(e) -Inf)
  }, th, m)
})

cat(sprintf("\nθ que maximiza a verossimilhança de Ewens: %.2f\n",
            theta_grid[which.max(rowMeans(loglik_grid))]))

# Se 'untb' estiver disponível, usar estimador mais preciso
if (untb_disponivel) {
  comunidade_untb <- count(comunidade_sim$abundancia)
  theta_mle_untb  <- theta.prob(comunidade_untb)
  cat(sprintf("θ estimado via untb::theta.prob: %.3f\n", theta_mle_untb))
}

# =============================================================================
# 5. COMPARAÇÃO DE MODELOS DE SAD
# =============================================================================
# Compara ajuste do modelo neutro com série logarítmica e log-normal

# 5.1 Série logarítmica de Fisher
alpha_est <- estimar_alpha_fisher(S = S_obs, N = J)
cat(sprintf("\nAlpha de Fisher: %.3f (S_obs = %d, N = %d)\n", alpha_est, S_obs, J))

# 5.2 Calcular expected rank-abundances para cada modelo
rank_seq <- seq_along(rank_abd)

# Broken stick (MacArthur)
bs_abd <- sad_broken_stick(S_obs, J)

# Comparar graficamente
df_comp <- rbind(
  data.frame(rank = rank_seq, abundancia = rank_abd,      modelo = "Observado"),
  data.frame(rank = rank_seq, abundancia = bs_abd,         modelo = "Broken Stick")
)

p_comp <- ggplot(df_comp, aes(x = rank, y = abundancia, color = modelo,
                               linetype = modelo)) +
  geom_line(linewidth = 0.9) +
  scale_y_log10() +
  scale_color_manual(values = c("Observado" = "#1A1A1A",
                                "Broken Stick" = "#E05C2A")) +
  labs(
    title    = "Comparação de Modelos de SAD",
    subtitle = "Comunidade neutra simulada vs. Broken Stick",
    x        = "Rank de abundância",
    y        = "Abundância [log10]",
    color    = "Modelo",
    linetype = "Modelo"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

print(p_comp)

# =============================================================================
# 6. EFEITO DO ISOLAMENTO NO PARÂMETRO m
# =============================================================================
# Em murundus mais isolados, esperamos menor m.
# Aqui simulamos comunidades com diferentes m (correspondendo a diferentes
# graus de isolamento) e plotamos a riqueza esperada.

set.seed(2024)
cenarios_m <- data.frame(
  isolamento  = c(10, 50, 100, 200, 400, 800),  # distância em metros
  m           = c(0.40, 0.30, 0.20, 0.12, 0.06, 0.02)
)
cenarios_m$riqueza_esperada <- sapply(cenarios_m$m, function(mi) {
  I_val <- mi * (J - 1) / (1 - mi)
  round(theta_verdadeiro * log(1 + J / theta_verdadeiro))  # aprox. simples
})

# Estimativa mais precisa usando simulação
set.seed(42)
riquezas_sim <- sapply(cenarios_m$m, function(mi) {
  com <- simular_comunidade_neutra(J, theta_verdadeiro, mi, t_max = J * 30)
  nrow(com)
})
cenarios_m$riqueza_simulada <- riquezas_sim

p_isol <- ggplot(cenarios_m, aes(x = isolamento, y = riqueza_simulada)) +
  geom_line(color = "#2C7A3F", linewidth = 1.2) +
  geom_point(size = 3, color = "#2C7A3F") +
  labs(
    title    = "Efeito do Isolamento na Riqueza: Perspectiva Neutra",
    subtitle = sprintf("θ = %g fixo; m varia de %.2f a %.2f",
                       theta_verdadeiro, max(cenarios_m$m), min(cenarios_m$m)),
    x        = "Isolamento do murundu (metros)",
    y        = "Riqueza de espécies (S)"
  ) +
  theme_classic(base_size = 13)

print(p_isol)
cat("\nCenários de isolamento vs. riqueza neutra:\n")
print(cenarios_m)
