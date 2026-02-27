# =============================================================================
# biogeografia_ilhas_murundus.R
# Projeto Murundus — Biogeografia de Ilhas
# =============================================================================
# Autores: Denis S. Nogueira
# Descrição: Análise da relação espécie-área (SAR), estimativa de taxas de
#   extinção esperadas, e teste do modelo de equilíbrio dinâmico de
#   MacArthur & Wilson (1967) aplicado a murundus do Cerrado.
#
# Referências principais:
#   MacArthur & Wilson (1967). The Theory of Island Biogeography. PUP.
#   Lomolino (2000). J. Biogeogr. 27:17-26.
#   Triantis et al. (2012). J. Biogeogr. 39:215-231.
# =============================================================================

# --- Dependências ---
library(ggplot2)
library(dplyr)

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
# Substitua por seus dados reais:
#   dados <- read.csv("../../dados/murundus_field_data.csv")
# Aqui usamos dados simulados para demonstração:
dados <- simular_dados_murundus(n_murundus = 50, semente = 2024)

# =============================================================================
# 2. RELAÇÃO ESPÉCIE-ÁREA (SAR)
# =============================================================================

# 2.1 Ajustar modelo potência S = c * A^z no espaço log-log
sar <- ajustar_sar(dados$area, dados$riqueza)
cat(sprintf("SAR — c = %.3f, z = %.3f, R² = %.3f, AIC = %.1f\n",
            sar$c, sar$z, sar$R2, sar$AIC))

# 2.2 Comparar com modelos alternativos de SAR
# Modelo linear (S ~ A)
lm_linear  <- lm(riqueza ~ area, data = dados)
# Modelo logarítmico (S ~ log A) — Arrhenius alternativo
lm_log     <- lm(riqueza ~ log(area), data = dados)
# Modelo potência (já ajustado via log-log)
lm_loglog  <- lm(log(riqueza) ~ log(area), data = dados)

tabela_aic <- data.frame(
  Modelo = c("Linear", "Logarítmico", "Potência (log-log)"),
  AIC    = c(AIC(lm_linear), AIC(lm_log), AIC(lm_loglog)),
  R2     = c(summary(lm_linear)$r.squared,
             summary(lm_log)$r.squared,
             summary(lm_loglog)$r.squared)
)
tabela_aic <- tabela_aic[order(tabela_aic$AIC), ]
cat("\nComparação de modelos SAR:\n")
print(tabela_aic)

# 2.3 Visualização: SAR no espaço log-log
area_seq   <- seq(min(dados$area), max(dados$area), length.out = 200)
riqueza_pred <- sar_potencia(area_seq, sar$c, sar$z)

p_sar <- ggplot(dados, aes(x = area, y = riqueza)) +
  geom_point(size = 2, alpha = 0.7, color = "#2C7A3F") +
  geom_line(data = data.frame(area = area_seq, riqueza = riqueza_pred),
            aes(x = area, y = riqueza), color = "#E05C2A", linewidth = 1.2) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks() +
  labs(
    title    = "Relação Espécie-Área nos Murundus",
    subtitle = sprintf("S = %.2f × A^{%.3f}   (R² = %.3f)",
                       sar$c, sar$z, sar$R2),
    x        = "Área do murundu (m²) [escala log]",
    y        = "Riqueza de espécies (S) [escala log]",
    caption  = "Fonte: Projeto Murundus"
  ) +
  theme_classic(base_size = 13)

print(p_sar)
# ggsave("../../resultados/sar_murundus.png", p_sar, width = 8, height = 6, dpi = 300)

# =============================================================================
# 3. EXTINÇÃO ESPERADA POR REDUÇÃO DE ÁREA
# =============================================================================
# Cenários de desmatamento do entorno dos murundus que reduzem a área efetiva
z_est     <- sar$z
remanesc  <- c(0.90, 0.75, 0.50, 0.25, 0.10)  # frações de área remanescente

cenarios <- data.frame(
  fracao_remanescente    = remanesc,
  extincao_esperada_frac = sapply(remanesc, function(r)
    extincao_por_area(A1 = 1, A2 = r, z = z_est)),
  extincao_esperada_pct  = sapply(remanesc, function(r)
    100 * extincao_por_area(A1 = 1, A2 = r, z = z_est))
)
cat("\nExtinções esperadas pela SAR (z =", round(z_est, 3), "):\n")
print(cenarios)

# Visualização: débito de extinção
p_extinc <- ggplot(cenarios, aes(x = fracao_remanescente * 100,
                                  y = extincao_esperada_pct)) +
  geom_area(fill = "#C0392B", alpha = 0.3) +
  geom_line(color = "#C0392B", linewidth = 1.2) +
  geom_point(size = 3, color = "#C0392B") +
  scale_x_reverse() +
  labs(
    title    = "Extinção Esperada pela Redução de Área dos Murundus",
    subtitle = sprintf("Usando z = %.3f da SAR observada", z_est),
    x        = "Área remanescente (%)",
    y        = "Espécies esperadas extintas (%)"
  ) +
  theme_classic(base_size = 13)

print(p_extinc)
# ggsave("../../resultados/extincao_sar_murundus.png", p_extinc,
#        width = 8, height = 6, dpi = 300)

# =============================================================================
# 4. EFEITO DO ISOLAMENTO (MODELO DE EQUILÍBRIO DINÂMICO)
# =============================================================================
# A TIB prediz que murundus mais isolados têm menor taxa de imigração e,
# portanto, menor riqueza de equilíbrio (para mesma área).
# Testamos com regressão múltipla em escala log:
#   log(S) ~ log(A) + log(isolamento)
# O coeficiente de log(isolamento) deve ser negativo.

modelo_tib <- lm(log(riqueza) ~ log(area) + log(isolamento), data = dados)
cat("\nModelo TIB completo (área + isolamento):\n")
print(summary(modelo_tib))

# Comparar AIC: com vs. sem isolamento
delta_aic <- AIC(lm_loglog) - AIC(modelo_tib)
cat(sprintf("\nΔAIC (modelo com isolamento — sem isolamento) = %.2f\n", delta_aic))
cat("(negativo = isolamento melhora o ajuste)\n")

# =============================================================================
# 5. ANÁLISE DE INCIDÊNCIA (FUNÇÃO DE INCIDÊNCIA DE DIAMOND)
# =============================================================================
# Para cada espécie, a probabilidade de ocorrência j(S) é plotada contra
# a riqueza do murundu em que foi registrada.
# Requer: dados de presença/ausência por murundu e por espécie.
# (Exemplo com dados simulados)

set.seed(123)
n_spp <- 20
n_mur <- nrow(dados)
# Simular prob. de ocorrência dependente da riqueza do murundu
prob_ocorr <- outer(seq(0.05, 0.95, length.out = n_spp),
                    normalizar(dados$riqueza), function(p, r) pmin(p * (1 + r), 1))
matriz_pa  <- matrix(rbinom(n_spp * n_mur, 1, prob_ocorr),
                     nrow = n_spp, ncol = n_mur,
                     dimnames = list(paste0("sp", seq_len(n_spp)),
                                     paste0("mur", seq_len(n_mur))))

# Calcular função de incidência por espécie
incidencia <- data.frame(
  especie      = rownames(matriz_pa),
  n_murundus   = rowSums(matriz_pa),
  prop_ocorr   = rowMeans(matriz_pa)
)
incidencia <- incidencia[order(-incidencia$prop_ocorr), ]

cat("\nFunções de incidência (top 10 espécies):\n")
print(head(incidencia, 10))

# Classificação Diamond: espécies de super-tramp (alta incidência em ilhas
# pequenas), espécies de interior (apenas em ilhas grandes)
incidencia$categoria <- cut(incidencia$prop_ocorr,
                             breaks = c(0, 0.33, 0.67, 1),
                             labels = c("interior", "generalista", "super-tramp"))
cat("\nDistribuição por categoria Diamond:\n")
print(table(incidencia$categoria))
