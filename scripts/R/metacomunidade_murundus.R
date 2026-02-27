# =============================================================================
# metacomunidade_murundus.R
# Projeto Murundus — Ecologia de Metacomunidades
# =============================================================================
# Autores: Denis S. Nogueira
# Descrição: Análise da estrutura de metacomunidade em murundus usando:
#   1. Partição de diversidade beta (turnover vs. nestedness)
#   2. Partição de variância (ambiental vs. espacial via db-RDA + MEM)
#   3. Teste de aninhamento (NODF)
#   4. Classificação segundo os paradigmas de Leibold et al. (2004)
#
# Referências principais:
#   Leibold et al. (2004). Ecology Letters 7:601-613.
#   Baselga (2010). Glob. Ecol. Biogeogr. 19:134-143.
#   Dray et al. (2006). Ecol. Model. 196:483-493.
# =============================================================================

# --- Dependências ---
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

# Pacotes opcionais
for (pkg in c("betapart", "adespatial")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    message(sprintf("Instale o pacote '%s' para análises completas.", pkg))
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
# Substitua pela sua matriz real de comunidade (murundus × espécies) e
# pelas variáveis ambientais e coordenadas espaciais.

set.seed(2024)
n_mur <- 30
n_spp <- 25

# Simular gradiente ambiental (ex.: umidade do solo, % cobertura arbórea)
env <- data.frame(
  umidade   = runif(n_mur, 20, 80),    # % umidade volumétrica
  altura    = rlnorm(n_mur, 1, 0.4),  # altura do murundu (m)
  cobertura = runif(n_mur, 5, 60)     # % cobertura arbórea
)
coords <- data.frame(
  x = runif(n_mur, 0, 1000),  # UTM leste (metros)
  y = runif(n_mur, 0, 1000)   # UTM norte (metros)
)

# Simular matriz de presença/ausência condicionada ao gradiente ambiental
niche_pos <- runif(n_spp, 20, 80)    # centro de nicho de cada espécie em umidade
niche_bw  <- runif(n_spp, 10, 30)   # largura de nicho

prob_mat <- sapply(seq_len(n_spp), function(j) {
  pmin(1, exp(-((env$umidade - niche_pos[j])^2) / (2 * niche_bw[j]^2)) * 0.9 + 0.05)
})
matriz_pa <- matrix(
  rbinom(n_mur * n_spp, 1, prob_mat),
  nrow = n_mur, ncol = n_spp,
  dimnames = list(paste0("mur", seq_len(n_mur)), paste0("sp", seq_len(n_spp)))
)
cat(sprintf("Matriz PA: %d murundus × %d espécies  |  %d ocorrências\n",
            n_mur, n_spp, sum(matriz_pa)))

# =============================================================================
# 2. DIVERSIDADE ALFA, BETA E GAMA
# =============================================================================

alpha_diversidade <- rowSums(matriz_pa)
gamma_diversidade <- sum(colSums(matriz_pa) > 0)
beta_w            <- beta_whittaker(matriz_pa)

cat(sprintf("\nα-médio = %.1f  |  γ = %d  |  β_Whittaker = %.3f\n",
            mean(alpha_diversidade), gamma_diversidade, beta_w))

p_alpha <- ggplot(data.frame(murundu = seq_len(n_mur), riqueza = alpha_diversidade),
                  aes(x = riqueza)) +
  geom_histogram(binwidth = 1, fill = "#2C7A3F", color = "white") +
  labs(title = "Distribuição de Riqueza Alpha (por murundu)",
       x = "Riqueza de espécies (S)", y = "Frequência") +
  theme_classic(base_size = 13)
print(p_alpha)

# =============================================================================
# 3. PARTIÇÃO DE DIVERSIDADE BETA (betapart)
# =============================================================================

if (requireNamespace("betapart", quietly = TRUE)) {
  library(betapart)
  beta_core  <- betapart.core(matriz_pa)
  beta_multi <- beta.multi(beta_core)

  cat("\nPartição de β-diversidade (Baselga 2010):\n")
  cat(sprintf("  β_total (Sørensen)    = %.4f\n", beta_multi$beta.SOR))
  cat(sprintf("  β_nestedness          = %.4f\n", beta_multi$beta.SNE))
  cat(sprintf("  β_turnover (Simpson)  = %.4f\n", beta_multi$beta.SIM))
  cat(sprintf("  Proporção nestedness  = %.1f%%\n",
              100 * beta_multi$beta.SNE / beta_multi$beta.SOR))
  cat(sprintf("  Proporção turnover    = %.1f%%\n",
              100 * beta_multi$beta.SIM / beta_multi$beta.SOR))

  # Mapa de dissimilaridade par-a-par
  beta_par <- beta.pair(beta_core)
  beta_df  <- as.data.frame(as.matrix(beta_par$beta.sor))
  cat("\nMatriz de dissimilaridade de Sørensen (primeiros 5 × 5):\n")
  print(round(beta_df[1:5, 1:5], 3))

} else {
  cat("\nInstale 'betapart' para partição de beta-diversidade.\n")
}

# =============================================================================
# 4. ORDENAÇÃO: NMDS
# =============================================================================

# Dissimilaridade de Bray-Curtis (ou Jaccard para P/A)
dist_jacc <- vegdist(matriz_pa, method = "jaccard", binary = TRUE)

set.seed(42)
nmds <- metaMDS(dist_jacc, k = 2, trymax = 100, trace = FALSE)
cat(sprintf("\nNMDS: stress = %.4f  (< 0.10 = bom; < 0.20 = aceitável)\n",
            nmds$stress))

nmds_df <- as.data.frame(scores(nmds, display = "sites"))
nmds_df$umidade <- env$umidade

p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = umidade)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "#D4A017", high = "#1A6B3C",
                       name = "Umidade (%)") +
  labs(
    title    = "NMDS — Composição de Espécies nos Murundus",
    subtitle = sprintf("Dissimilaridade de Jaccard  |  Stress = %.3f", nmds$stress),
    x = "NMDS1", y = "NMDS2"
  ) +
  theme_classic(base_size = 13)
print(p_nmds)

# =============================================================================
# 5. PARTIÇÃO DE VARIÂNCIA (ambiental vs. espacial)
# =============================================================================
# Abordagem db-RDA (Redundancy Analysis baseada em distâncias)
# Variáveis espaciais: coordenadas x e y (ou vetores MEM via adespatial)

# 5.1 Construir variáveis espaciais (PCNMs simples via vegan)
pcnm_coords <- pcnm(dist(coords))
pcnm_vars   <- as.data.frame(scores(pcnm_coords))
# Manter apenas PCNMs positivos (padrão)
pcnm_pos    <- pcnm_vars[, pcnm_coords$values > 0, drop = FALSE]
cat(sprintf("\nVariáveis PCNM positivas: %d\n", ncol(pcnm_pos)))

# 5.2 Selecionar variáveis ambientais relevantes por forward selection
# (simplificado — usar vegan::ordistep em análise real)
env_std <- as.data.frame(scale(env))

# Testes simples de Mantel (correlação ambiente vs. composição)
mantel_env <- mantel(dist_jacc, dist(env_std), permutations = 499)
mantel_geo <- mantel(dist_jacc, dist(coords),  permutations = 499)
cat(sprintf("\nMantel: composição × ambiente — r = %.3f, p = %.3f\n",
            mantel_env$statistic, mantel_env$signif))
cat(sprintf("Mantel: composição × distância geográfica — r = %.3f, p = %.3f\n",
            mantel_geo$statistic, mantel_geo$signif))

# 5.3 Partição de variância formal (varpart)
# Requer transformação de Hellinger da matriz de PA
hellinger_pa <- decostand(matriz_pa, method = "hellinger")

vp <- varpart(hellinger_pa, env_std, pcnm_pos[, 1:min(5, ncol(pcnm_pos))])
cat("\nPartição de variância (db-RDA):\n")
print(vp)

# Interpretação dos fracionamentos:
# [a] = variação ambiental pura → filtragem de nicho
# [b] = variação ambiental+espacial confundida → difícil de interpretar
# [c] = variação espacial pura → limitação de dispersão / neutralidade
# [d] = variação não explicada → estocasticidade residual

# =============================================================================
# 6. CLASSIFICAÇÃO POR PARADIGMA DE METACOMUNIDADE
# =============================================================================
# Regras heurísticas (Leibold & Chase 2018):
#   - [a] > [c]  &  nestedness alta → Triagem de Espécies
#   - [c] > [a]  &  turnover alto   → Dinâmica de Manchas / Neutro
#   - [a] > [c]  &  turnover alto   → Efeitos de Massa

fracao_ambiental <- vp$part$indfract$Adj.R.squared[1]  # [a]
fracao_espacial  <- vp$part$indfract$Adj.R.squared[3]  # [c]

if (requireNamespace("betapart", quietly = TRUE)) {
  nest_prop <- beta_multi$beta.SNE / beta_multi$beta.SOR
  turn_prop <- beta_multi$beta.SIM / beta_multi$beta.SOR
} else {
  nest_prop <- NA
  turn_prop <- NA
}

cat("\n--- Classificação de Paradigma ---\n")
cat(sprintf("  Variação ambiental pura [a] = %.3f\n",
            max(fracao_ambiental, 0)))
cat(sprintf("  Variação espacial pura [c]  = %.3f\n",
            max(fracao_espacial, 0)))
if (!is.na(nest_prop)) {
  cat(sprintf("  Proporção de nestedness     = %.3f\n", nest_prop))
  cat(sprintf("  Proporção de turnover       = %.3f\n", turn_prop))
  if (fracao_ambiental > fracao_espacial && nest_prop > 0.5) {
    cat("  → Paradigma dominante: TRIAGEM DE ESPÉCIES\n")
  } else if (fracao_espacial > fracao_ambiental && turn_prop > 0.5) {
    cat("  → Paradigma dominante: DINÂMICA DE MANCHAS / NEUTRO\n")
  } else if (fracao_ambiental > fracao_espacial && turn_prop > 0.5) {
    cat("  → Paradigma dominante: EFEITOS DE MASSA\n")
  } else {
    cat("  → Padrão misto — múltiplos processos atuando\n")
  }
}
