# =============================================================================
# macroecologia_murundus.R
# Projeto Murundus — Macroecologia
# =============================================================================
# Autores: Denis S. Nogueira
# Descrição: Análise de padrões macroecológicos em murundus:
#   1. Distribuições de Abundância de Espécies (SADs): comparação de modelos
#   2. Relações de escalonamento: densidade × massa corporal (Lei de Damuth)
#   3. Relação riqueza-área e riqueza-produtividade
#   4. Gradientes de diversidade e Teoria Metabólica da Ecologia (MTE)
#
# Referências principais:
#   Brown & Maurer (1989). Science 243:1145-1150.
#   Brown et al. (2004). Ecology 85:1771-1789 (MTE).
#   Damuth (1987). Biol. J. Linn. Soc. 31:193-246.
#   McGill et al. (2007). Ecology Letters 10:995-1015. (SAD review)
# =============================================================================

# --- Dependências ---
library(ggplot2)
library(dplyr)

# iNEXT: rarefação e extrapolação (opcional)
if (requireNamespace("iNEXT", quietly = TRUE)) {
  library(iNEXT)
  inext_disp <- TRUE
} else {
  message("Instale 'iNEXT' para curvas de rarefação e extrapolação.")
  inext_disp <- FALSE
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
set.seed(2024)
n_spp <- 40

# Simular comunidade com estrutura realista (log-normal)
log_abd   <- rnorm(n_spp, mean = 3, sd = 1.2)
abundancias <- round(exp(log_abd))
J           <- sum(abundancias)
S           <- length(abundancias)

cat(sprintf("Pool regional simulado: %d espécies, %d indivíduos\n", S, J))

# Simular traços funcionais
massa_g <- exp(runif(n_spp, log(0.1), log(500)))  # 0.1 g a 500 g
area_corporal_cm2 <- 0.1 * massa_g^0.67  # relação empírica aproximada

# Simular densidade populacional (Lei de Damuth: D ~ M^{-3/4})
dens_teorica   <- 100 * massa_g^(-0.75)   # ind/ha (normalização arbitrária)
dens_observada <- dens_teorica * rlnorm(n_spp, 0, 0.4)  # + ruído log-normal

spp_df <- data.frame(
  especie    = paste0("sp", seq_len(n_spp)),
  abundancia = abundancias,
  massa_g    = massa_g,
  dens_obs   = dens_observada
)

# =============================================================================
# 2. DISTRIBUIÇÕES DE ABUNDÂNCIA (SAD) — COMPARAÇÃO DE MODELOS
# =============================================================================

# 2.1 Série Logarítmica vs. Log-Normal vs. Broken Stick
alpha_est <- estimar_alpha_fisher(S, J)
bs_abd    <- sad_broken_stick(S, J)

cat(sprintf("\nAlpha de Fisher: %.3f\n", alpha_est))

# Plot comparativo em rank-abundance
rank_seq    <- seq_len(S)
rank_obs    <- sort(abundancias, decreasing = TRUE)

df_sad <- rbind(
  data.frame(rank = rank_seq, abundancia = rank_obs,
             modelo = "Observado (log-normal simulado)"),
  data.frame(rank = rank_seq, abundancia = sort(bs_abd, decreasing = TRUE),
             modelo = "Broken Stick (MacArthur)")
)

p_rank <- ggplot(df_sad, aes(x = rank, y = abundancia,
                              color = modelo, linetype = modelo)) +
  geom_line(linewidth = 0.9) +
  scale_y_log10() +
  scale_color_manual(values = c("Observado (log-normal simulado)" = "#1A1A1A",
                                "Broken Stick (MacArthur)"        = "#E05C2A")) +
  labs(
    title    = "Distribuição de Abundâncias (Rank-Abundance)",
    subtitle = sprintf("S = %d; N = %d; α_Fisher = %.2f", S, J, alpha_est),
    x        = "Rank de Abundância",
    y        = "Abundância [log10]",
    color    = NULL, linetype = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")
print(p_rank)

# 2.2 Octaves de Preston
log2_rank  <- floor(log2(rank_obs)) + 1
sad_press  <- data.frame(
  octave     = seq_len(max(log2_rank)),
  n_especies = tabulate(log2_rank, nbins = max(log2_rank))
)

p_preston <- ggplot(sad_press, aes(x = octave, y = n_especies)) +
  geom_col(fill = "#2C7A3F", color = "white") +
  labs(
    title = "Octaves de Preston (Log-Normal)",
    x     = expression("Classe de abundância (log"[2]*")"),
    y     = "Número de espécies"
  ) +
  theme_classic(base_size = 13)
print(p_preston)

# =============================================================================
# 3. LEI DE DAMUTH: DENSIDADE × MASSA CORPORAL
# =============================================================================

# Ajuste por mínimos quadrados em log-log
lm_damuth <- lm(log(dens_obs) ~ log(massa_g), data = spp_df)
coef_dam  <- coef(lm_damuth)
exp_dam   <- coef_dam[2]
r2_dam    <- summary(lm_damuth)$r.squared

cat(sprintf("\nLei de Damuth (OLS log-log):\n"))
cat(sprintf("  Expoente empírico = %.3f  (esperado teórico = -0.750)\n", exp_dam))
cat(sprintf("  R² = %.3f\n", r2_dam))

# Linha teórica de Damuth
massa_seq <- exp(seq(log(min(spp_df$massa_g)), log(max(spp_df$massa_g)),
                     length.out = 200))
dens_teorica_seq <- exp(coef_dam[1]) * massa_seq^(-0.75)
dens_empir_seq   <- exp(coef_dam[1]) * massa_seq^exp_dam

p_damuth <- ggplot(spp_df, aes(x = massa_g, y = dens_obs)) +
  geom_point(size = 2, alpha = 0.6, color = "#2C7A3F") +
  geom_line(data = data.frame(massa_g = massa_seq, dens = dens_teorica_seq),
            aes(x = massa_g, y = dens), color = "#C0392B",
            linetype = "dashed", linewidth = 1) +
  geom_line(data = data.frame(massa_g = massa_seq, dens = dens_empir_seq),
            aes(x = massa_g, y = dens), color = "#1A6B3C", linewidth = 1) +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks() +
  labs(
    title    = "Lei de Damuth: Densidade × Massa Corporal",
    subtitle = sprintf("Expoente empírico = %.3f  |  Teórico (MTE) = -0.750",
                       exp_dam),
    x        = "Massa corporal (g) [log10]",
    y        = "Densidade populacional (ind/ha) [log10]",
    caption  = "Vermelho tracejado = Damuth teórico (-3/4)  |  Verde = ajuste empírico"
  ) +
  theme_classic(base_size = 13)
print(p_damuth)

# =============================================================================
# 4. TEORIA METABÓLICA DA ECOLOGIA (MTE) — TEMPERATURA E DIVERSIDADE
# =============================================================================
# Predição da MTE: ln(S) ~ -E/(k*T) onde E ≈ 0.65 eV
# Nos murundus: temperatura do solo varia entre a superfície dos murundus
# e a matriz; este gradiente de temperatura deve correlacionar com riqueza.

k_boltzmann <- 8.617333e-5  # eV/K
E_ativacao  <- 0.65          # eV (energia de ativação metabólica)
T_celsius   <- seq(25, 42, by = 1)
T_kelvin    <- T_celsius + 273.15
inv_kT      <- 1 / (k_boltzmann * T_kelvin)

# Predição MTE: S ∝ exp(-E/kT) ∝ exp(E/kT) para crescente com temperatura
S_mte <- exp(E_ativacao * inv_kT - min(E_ativacao * inv_kT))
S_mte_norm <- S_mte / max(S_mte) * 20  # normalizado para escala de 0-20 spp

df_mte <- data.frame(
  T_celsius = T_celsius,
  S_mte     = S_mte_norm
)

p_mte <- ggplot(df_mte, aes(x = T_celsius, y = S_mte)) +
  geom_line(color = "#E05C2A", linewidth = 1.3) +
  geom_point(size = 2, color = "#E05C2A") +
  labs(
    title    = "Predição da MTE: Riqueza × Temperatura do Solo",
    subtitle = sprintf("E = %.2f eV (Teoria Metabólica — Brown et al. 2004)", E_ativacao),
    x        = "Temperatura do solo (°C)",
    y        = "Riqueza de espécies predita (normalizada)"
  ) +
  theme_classic(base_size = 13)
print(p_mte)

# =============================================================================
# 5. RAREFAÇÃO E CURVAS DE ACUMULAÇÃO (iNEXT)
# =============================================================================
if (inext_disp) {
  cat("\nCalculando curvas de rarefação com iNEXT...\n")
  # Usar lista de comunidades simuladas (múltiplos murundus)
  set.seed(42)
  comunidades_iNEXT <- lapply(seq_len(5), function(i) {
    n_i <- sample(5:15, 1)
    abd <- round(rlnorm(n_i, meanlog = 2, sdlog = 1))
    sort(abd, decreasing = TRUE)
  })
  names(comunidades_iNEXT) <- paste0("Murundu_", seq_along(comunidades_iNEXT))

  out_iNEXT <- iNEXT(comunidades_iNEXT, q = 0, datatype = "abundance",
                     endpoint = 3 * max(sapply(comunidades_iNEXT, sum)))
  p_inext <- ggiNEXT(out_iNEXT, type = 1) +
    labs(
      title    = "Curvas de Rarefação e Extrapolação (iNEXT)",
      subtitle = "Diversidade de ordem q = 0 (riqueza)",
      x        = "Número de indivíduos amostrados",
      y        = "Riqueza de espécies (S)"
    ) +
    theme_classic(base_size = 13)
  print(p_inext)
} else {
  cat("\nSkipping iNEXT analysis (package not installed).\n")
}
