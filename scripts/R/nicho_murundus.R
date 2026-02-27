# =============================================================================
# nicho_murundus.R
# Projeto Murundus — Teoria de Nicho
# =============================================================================
# Autores: Denis S. Nogueira
# Descrição: Análise do espaço de nicho para espécies associadas a murundus:
#   1. Amplitude e sobreposição de nicho (Levins/Pianka)
#   2. Estimativa de volume do hipervolume de Hutchinson
#   3. Filtragem ambiental: comparação de traços funcionais vs. nulo
#   4. Partição entre efeitos de estabilização e equalização (Chesson 2000)
#
# Referências principais:
#   Hutchinson (1957). Cold Spring Harbor Symp. 22:415-427.
#   Levins (1968). Evolution in Changing Environments. PUP.
#   Pianka (1973). Annu. Rev. Ecol. Syst. 4:53-74.
#   Chesson (2000). Annu. Rev. Ecol. Syst. 31:343-366.
# =============================================================================

# --- Dependências ---
library(ggplot2)
library(dplyr)
library(tidyr)

# Pacote hypervolume (opcional — para hipervolumes multivariados)
if (requireNamespace("hypervolume", quietly = TRUE)) {
  library(hypervolume)
  hypervolume_disp <- TRUE
} else {
  message("Instale 'hypervolume' para cálculo de hipervolumes multivariados.")
  hypervolume_disp <- FALSE
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
# Substitua por dados reais de:
#   - ocorrencias: data.frame de ocorrências com variáveis ambientais
#   - tracos:      data.frame de traços funcionais das espécies
#   - utilizacao:  matriz de espécies × recursos utilizados

set.seed(2024)
n_spp <- 15

# Simular gradiente ambiental de dois eixos:
# Eixo 1: umidade (murundu seco vs. matriz úmida)
# Eixo 2: temperatura do solo (murundu exposto vs. sombreado)
eixo1 <- seq(20, 80, length.out = 100)  # umidade (%)
eixo2 <- seq(25, 45, length.out = 100)  # temperatura (°C)

# Nicho fundamental de cada espécie: elipse no espaço 2D
niche_params <- data.frame(
  especie    = paste0("sp", seq_len(n_spp)),
  mu_umid    = runif(n_spp, 25, 75),
  mu_temp    = runif(n_spp, 27, 43),
  sd_umid    = runif(n_spp, 5, 20),
  sd_temp    = runif(n_spp, 2, 8)
)

# Simular ocorrências de cada espécie ao longo do gradiente
n_ocorr  <- 50
ocorr_list <- lapply(seq_len(n_spp), function(j) {
  data.frame(
    especie  = paste0("sp", j),
    umidade  = rnorm(n_ocorr, niche_params$mu_umid[j], niche_params$sd_umid[j]),
    temp_solo = rnorm(n_ocorr, niche_params$mu_temp[j], niche_params$sd_temp[j])
  )
})
ocorr_df <- do.call(rbind, ocorr_list)
ocorr_df <- ocorr_df[ocorr_df$umidade >= 0 & ocorr_df$umidade <= 100 &
                       ocorr_df$temp_solo >= 20 & ocorr_df$temp_solo <= 50, ]

# Simular matriz de utilização de recursos (5 recursos: tipos de solo/alimento)
n_recursos <- 5
utilizacao <- matrix(
  rpois(n_spp * n_recursos, lambda = rep(c(5, 3, 8, 2, 6), each = n_spp)),
  nrow = n_spp, ncol = n_recursos,
  dimnames = list(paste0("sp", seq_len(n_spp)),
                  paste0("recurso", seq_len(n_recursos)))
)

# =============================================================================
# 2. AMPLITUDE DE NICHO (LEVINS 1968)
# =============================================================================
amplitude_df <- data.frame(
  especie    = paste0("sp", seq_len(n_spp)),
  B_Levins   = apply(utilizacao, 1, amplitude_nicho_levins),
  n_recursos = n_recursos
)
amplitude_df$B_padronizado <- amplitude_df$B_Levins / amplitude_df$n_recursos

cat("Amplitude de Nicho de Levins (B_padronizado: 1/K = especialista; 1 = generalista):\n")
print(amplitude_df[order(-amplitude_df$B_padronizado), ])

p_amplitude <- ggplot(amplitude_df, aes(x = reorder(especie, B_padronizado),
                                         y = B_padronizado)) +
  geom_col(fill = "#2C7A3F") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#C0392B") +
  coord_flip() +
  labs(
    title    = "Amplitude de Nicho de Levins (Padronizada)",
    subtitle = "Linha tracejada = generalismo intermediário",
    x        = "Espécie",
    y        = "B padronizado (1/K → K/K)"
  ) +
  theme_classic(base_size = 12)
print(p_amplitude)

# =============================================================================
# 3. SOBREPOSIÇÃO DE NICHO (PIANKA 1973)
# =============================================================================
overlap_mat <- sobreposicao_nicho(utilizacao)

# Média de sobreposição por espécie (excluindo diagonal)
overlap_mean <- apply(overlap_mat, 1, function(x) mean(x[x < 1]))

cat("\nSobreposição média de nicho por espécie (Pianka):\n")
print(round(data.frame(especie = rownames(overlap_mat),
                        overlap_medio = overlap_mean), 3))

# Heatmap de sobreposição
overlap_long <- as.data.frame(overlap_mat) |>
  tibble::rownames_to_column("sp1") |>
  pivot_longer(-sp1, names_to = "sp2", values_to = "overlap")

p_overlap <- ggplot(overlap_long, aes(x = sp1, y = sp2, fill = overlap)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#FFFFFF", mid = "#A8D5A2", high = "#1A6B3C",
                        midpoint = 0.5, limits = c(0, 1),
                        name = "Sobreposição\n(Pianka)") +
  labs(title = "Sobreposição de Nicho par-a-par (Pianka 1973)",
       x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_overlap)

# =============================================================================
# 4. HIPERVOLUME DE HUTCHINSON
# =============================================================================
# 4.1 Estimativa simplificada via elipsoide
cat("\nVolume do hipervolume de Hutchinson (elipsoide) por espécie:\n")
vol_df <- ocorr_df |>
  group_by(especie) |>
  summarise(
    volume_elips = volume_nicho_elipsoide(cbind(umidade, temp_solo)),
    n_ocorr      = n(),
    .groups = "drop"
  )
print(vol_df)

# 4.2 Hipervolume multivariado via kernel (se pacote disponível)
if (hypervolume_disp) {
  cat("\nCalculando hipervolumes via kernel (pacote 'hypervolume')...\n")
  # Exemplo com as 3 primeiras espécies
  spp_exemplo <- c("sp1", "sp2", "sp3")
  hv_list <- lapply(spp_exemplo, function(sp) {
    dados_sp <- ocorr_df[ocorr_df$especie == sp, c("umidade", "temp_solo")]
    dados_sp <- scale(dados_sp)
    tryCatch(
      hypervolume_gaussian(as.data.frame(dados_sp), verbose = FALSE),
      error = function(e) NULL
    )
  })

  volumes_kernel <- sapply(hv_list, function(hv) {
    if (is.null(hv)) NA_real_ else get_volume(hv)
  })
  cat(sprintf("Volume do hipervolume (sp1, sp2, sp3): %.2f, %.2f, %.2f\n",
              volumes_kernel[1], volumes_kernel[2], volumes_kernel[3]))
}

# =============================================================================
# 5. FILTRAGEM AMBIENTAL vs. LIMITAÇÃO DE DISPERSÃO
# =============================================================================
# Teste: a dispersão de traços funcionais nas comunidades de murundus é
# menor que o esperado ao acaso? → filtragem ambiental
# Maior? → repulsão de nicho (over-dispersion) / limitação de dispersão

# Usar comprimento de traço (leaf length, body size, etc.) simulado
tracos <- data.frame(
  especie         = paste0("sp", seq_len(n_spp)),
  comprimento_cm  = rlnorm(n_spp, meanlog = 1.5, sdlog = 0.5),
  massa_g         = rlnorm(n_spp, meanlog = 2.0, sdlog = 0.7)
)

# Simular presença/ausência filtrada por ambiente (murundus secos vs. úmidos)
set.seed(42)
n_comunidades <- 20
pa_tracos <- matrix(
  rbinom(n_comunidades * n_spp, 1, 0.4),
  nrow = n_comunidades, ncol = n_spp,
  dimnames = list(paste0("com", seq_len(n_comunidades)),
                  paste0("sp", seq_len(n_spp)))
)

# Calcular variância de traço por comunidade
tracos_vec <- setNames(tracos$comprimento_cm, tracos$especie)
var_traco_obs <- apply(pa_tracos, 1, function(linha) {
  spp_presentes <- names(linha[linha == 1])
  if (length(spp_presentes) < 2) return(NA_real_)
  var(tracos_vec[spp_presentes])
})

# Distribuição nula (1000 permutações)
n_perm  <- 999
var_nulo <- replicate(n_perm, {
  pa_nulo <- apply(pa_tracos, 2, sample)  # aleatorizar colunas
  apply(pa_nulo, 1, function(linha) {
    spp_presentes <- names(linha[linha == 1])
    if (length(spp_presentes) < 2) return(NA_real_)
    var(tracos_vec[spp_presentes])
  })
})

# Standardized Effect Size (SES): < -2 = filtragem; > +2 = repulsão
ses_df <- data.frame(
  comunidade = names(var_traco_obs),
  var_obs    = var_traco_obs,
  var_nulo_media = rowMeans(var_nulo, na.rm = TRUE),
  var_nulo_sd    = apply(var_nulo, 1, sd, na.rm = TRUE)
)
ses_df$SES <- (ses_df$var_obs - ses_df$var_nulo_media) / ses_df$var_nulo_sd

cat("\nSES (Standardized Effect Size) de variância de traço:\n")
print(summary(ses_df$SES))
cat(sprintf("Proporção de comunidades com filtragem (SES < -2): %.1f%%\n",
            100 * mean(ses_df$SES < -2, na.rm = TRUE)))
cat(sprintf("Proporção de comunidades com repulsão (SES > +2): %.1f%%\n",
            100 * mean(ses_df$SES > 2, na.rm = TRUE)))

p_ses <- ggplot(ses_df, aes(x = SES)) +
  geom_histogram(binwidth = 0.5, fill = "#2C7A3F", color = "white") +
  geom_vline(xintercept = c(-2, 2), color = "#C0392B", linetype = "dashed") +
  labs(
    title    = "SES de Variância de Traço (Filtragem Ambiental)",
    subtitle = "SES < -2: filtragem ambiental  |  SES > +2: repulsão de nicho",
    x        = "SES",
    y        = "Frequência"
  ) +
  theme_classic(base_size = 13)
print(p_ses)
