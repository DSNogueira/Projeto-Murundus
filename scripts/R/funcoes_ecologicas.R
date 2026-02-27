# =============================================================================
# funcoes_ecologicas.R
# Projeto Murundus — Funções Utilitárias de Ecologia Teórica
# =============================================================================
# Autores: Denis S. Nogueira
# Descrição: Funções reutilizáveis para análises de biogeografia de ilhas,
#   teoria neutra, teoria de nicho, metacomunidades e macroecologia.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. RELAÇÃO ESPÉCIE-ÁREA (SAR)
# -----------------------------------------------------------------------------

#' Calcular riqueza predita pela relação potência S = c * A^z
#'
#' @param area vetor numérico de áreas (mesma unidade em todos os sítios)
#' @param c intercepto da SAR (escalar positivo)
#' @param z expoente da SAR (tipicamente 0.10–0.40)
#' @return vetor de riqueza predita
sar_potencia <- function(area, c, z) {
  c * area^z
}

#' Ajustar SAR pelo método dos mínimos quadrados (escala log-log)
#'
#' @param area vetor numérico de áreas
#' @param riqueza vetor numérico de riqueza observada
#' @return lista com coeficientes (c, z), R², AIC e objeto lm
ajustar_sar <- function(area, riqueza) {
  stopifnot(length(area) == length(riqueza), all(area > 0), all(riqueza > 0))
  modelo <- lm(log(riqueza) ~ log(area))
  coef_log <- coef(modelo)
  list(
    c   = exp(coef_log[1]),
    z   = coef_log[2],
    R2  = summary(modelo)$r.squared,
    AIC = AIC(modelo),
    lm  = modelo
  )
}

#' Calcular taxa de extinção esperada pela perda de área (SAR)
#'
#' Estima a fração de espécies perdidas quando a área é reduzida de A1 para A2,
#' usando a relação potência:  delta_S / S1 = 1 - (A2/A1)^z
#'
#' @param A1 área original
#' @param A2 área remanescente
#' @param z expoente da SAR
#' @return fração de espécies esperadas extintas (0–1)
extincao_por_area <- function(A1, A2, z) {
  stopifnot(A2 <= A1, z > 0)
  1 - (A2 / A1)^z
}

# -----------------------------------------------------------------------------
# 2. DISTRIBUIÇÃO DE ABUNDÂNCIAS DE ESPÉCIES (SAD)
# -----------------------------------------------------------------------------

#' Calcular abundâncias esperadas pelo modelo Broken Stick (MacArthur 1957)
#'
#' @param S número de espécies
#' @param N número total de indivíduos
#' @return vetor ordenado decrescente de abundâncias esperadas
sad_broken_stick <- function(S, N) {
  rank <- seq_len(S)
  n_i  <- (N / S) * sapply(rank, function(i) sum(1 / seq(i, S)))
  sort(n_i, decreasing = TRUE)
}

#' Calcular riqueza esperada pela Série Logarítmica de Fisher
#'
#' S = alpha * ln(1 + N/alpha),  onde alpha é o parâmetro de diversidade
#'
#' @param N número total de indivíduos
#' @param alpha parâmetro de diversidade de Fisher
#' @return riqueza esperada
riqueza_serie_log <- function(N, alpha) {
  alpha * log(1 + N / alpha)
}

#' Estimar alpha de Fisher iterativamente
#'
#' Resolve S = alpha * ln(1 + N/alpha) para alpha dado S e N
#'
#' @param S riqueza observada
#' @param N número total de indivíduos
#' @param tol tolerância de convergência
#' @param max_iter número máximo de iterações
#' @return estimativa de alpha
estimar_alpha_fisher <- function(S, N, tol = 1e-6, max_iter = 1000) {
  alpha <- S / log(N)  # valor inicial
  for (i in seq_len(max_iter)) {
    alpha_novo <- S / log(1 + N / alpha)
    if (abs(alpha_novo - alpha) < tol) return(alpha_novo)
    alpha <- alpha_novo
  }
  warning("Convergência não atingida em ", max_iter, " iterações")
  alpha
}

# -----------------------------------------------------------------------------
# 3. DIVERSIDADE
# -----------------------------------------------------------------------------

#' Calcular índices de diversidade clássicos
#'
#' @param abundancias vetor de abundâncias (inteiros ou proporções)
#' @return lista com Shannon (H'), Simpson (D), equitabilidade de Pielou (J')
diversidade <- function(abundancias) {
  abundancias <- abundancias[abundancias > 0]
  p  <- abundancias / sum(abundancias)
  H  <- -sum(p * log(p))
  D  <- 1 - sum(p^2)
  S  <- length(p)
  J  <- if (S > 1) H / log(S) else NA_real_
  list(Shannon = H, Simpson = D, Pielou_J = J, riqueza = S)
}

#' Calcular diversidade beta de Whittaker
#'
#' beta_W = gamma / mean(alpha) - 1
#'
#' @param matriz_comunidade matriz de sítios (linhas) × espécies (colunas)
#'        com contagens ou presença/ausência
#' @return diversidade beta de Whittaker (escalar)
beta_whittaker <- function(matriz_comunidade) {
  presenca <- matriz_comunidade > 0
  gamma    <- sum(colSums(presenca) > 0)
  alpha    <- mean(rowSums(presenca))
  gamma / alpha - 1
}

# -----------------------------------------------------------------------------
# 4. TEORIA NEUTRA (UNTB)
# -----------------------------------------------------------------------------

#' Calcular riqueza esperada pela UNTB na metacomunidade
#'
#' E[S] ≈ theta * ln(1 + J / theta)
#'
#' @param theta número fundamental de biodiversidade (θ = 2 * J_M * nu)
#' @param J tamanho da comunidade local (número de indivíduos)
#' @return riqueza esperada
riqueza_neutra <- function(theta, J) {
  theta * log(1 + J / theta)
}

#' Calcular probabilidade de Ewens (Ewens Sampling Formula) para uma amostra
#'
#' P(n_1, ..., n_S | J, theta) ∝ theta^S / (theta^{(J)}↑) * prod(1 / n_i)
#' Retorna o log-verossimilhança (para maximização numérica)
#'
#' @param abundancias vetor de abundâncias (soma = J)
#' @param theta parâmetro de diversidade da metacomunidade
#' @return log-verossimilhança (escalar)
loglik_ewens <- function(abundancias, theta) {
  abundancias <- abundancias[abundancias > 0]
  J <- sum(abundancias)
  S <- length(abundancias)
  log_pochhammer <- sum(log(theta + seq(0, J - 1)))
  S * log(theta) - log_pochhammer - sum(log(abundancias))
}

#' Estimar theta por máxima verossimilhança (Ewens)
#'
#' @param abundancias vetor de abundâncias observadas
#' @return lista com theta estimado e log-verossimilhança máxima
estimar_theta <- function(abundancias) {
  resultado <- optimize(
    f        = function(th) -loglik_ewens(abundancias, th),
    interval = c(1e-4, 1e4),
    tol      = 1e-6
  )
  list(theta = resultado$minimum, loglik = -resultado$objective)
}

#' Simular comunidade neutra (processo de morte-nascimento com imigração)
#'
#' @param J tamanho da comunidade local
#' @param theta número fundamental de biodiversidade
#' @param m taxa de imigração (0 < m < 1)
#' @param t_max número de passos de tempo
#' @param semente semente aleatória (opcional)
#' @return data.frame com composição final (espécie, abundância)
simular_comunidade_neutra <- function(J, theta, m, t_max = J * 10,
                                      semente = NULL) {
  if (!is.null(semente)) set.seed(semente)
  # Inicializar: cada indivíduo é uma espécie distinta
  comunidade  <- seq_len(J)
  proxima_spp <- J + 1L

  for (t in seq_len(t_max)) {
    # Escolhe um indivíduo para morrer
    morto <- sample.int(J, 1L)
    # Imigrante ou descendente local?
    if (runif(1) < m) {
      # Imigrante da metacomunidade (processo de Poisson com theta)
      # Com prob theta/(theta+J-1) é nova espécie; senão, cópia de existente
      prob_nova <- theta / (theta + J - 1)
      if (runif(1) < prob_nova) {
        comunidade[morto] <- proxima_spp
        proxima_spp <- proxima_spp + 1L
      } else {
        # Proporcional às abundâncias na metacomunidade (aprox. igual)
        comunidade[morto] <- comunidade[sample.int(J, 1L)]
      }
    } else {
      # Descendente local: escolha proporcional às abundâncias
      doador <- sample(setdiff(seq_len(J), morto), 1L)
      comunidade[morto] <- comunidade[doador]
    }
  }

  tab <- table(comunidade)
  data.frame(
    especie    = as.integer(names(tab)),
    abundancia = as.integer(tab)
  )
}

# -----------------------------------------------------------------------------
# 5. NICHO
# -----------------------------------------------------------------------------

#' Calcular amplitude de nicho (Levins 1968)
#'
#' B = 1 / sum(p_k^2), onde p_k = proporção do recurso k utilizado
#'
#' @param utilizacao vetor de utilização de recursos (não precisa somar 1)
#' @return amplitude de nicho de Levins (1 = especialista; K = generalista)
amplitude_nicho_levins <- function(utilizacao) {
  p <- utilizacao / sum(utilizacao)
  1 / sum(p^2)
}

#' Calcular sobreposição de nicho par-a-par (Pianka 1973)
#'
#' @param matriz_utilizacao matriz de espécies (linhas) × recursos (colunas)
#' @return matriz de sobreposição simétrica (0–1)
sobreposicao_nicho <- function(matriz_utilizacao) {
  n_spp <- nrow(matriz_utilizacao)
  props  <- sweep(matriz_utilizacao, 1, rowSums(matriz_utilizacao), "/")
  overlap <- matrix(NA_real_, n_spp, n_spp,
                    dimnames = list(rownames(props), rownames(props)))
  for (i in seq_len(n_spp)) {
    for (j in seq_len(n_spp)) {
      num <- sum(props[i, ] * props[j, ])
      den <- sqrt(sum(props[i, ]^2) * sum(props[j, ]^2))
      overlap[i, j] <- if (den > 0) num / den else 0
    }
  }
  overlap
}

#' Calcular volume do hipervolume de Hutchinson (estimativa simplificada)
#'
#' Aproximação pelo volume do hiperelipsoide definido pelos intervalos
#' [min, max] de cada variável ambiental nas ocorrências da espécie.
#'
#' @param ocorrencias matriz de ocorrências (linhas) × variáveis ambientais (colunas)
#' @return volume do hiperelipsoide (em unidades das variáveis)
volume_nicho_elipsoide <- function(ocorrencias) {
  n <- ncol(ocorrencias)
  semi_eixos <- apply(ocorrencias, 2, function(x) diff(range(x)) / 2)
  # Volume de hiperelipsoide: V_n = (pi^(n/2) / Gamma(n/2 + 1)) * prod(semi_eixos)
  vol_esfera <- pi^(n / 2) / gamma(n / 2 + 1)
  vol_esfera * prod(semi_eixos)
}

# -----------------------------------------------------------------------------
# 6. MÉTRICAS DE PAISAGEM SIMPLES
# -----------------------------------------------------------------------------

#' Calcular distância ao vizinho mais próximo (ENN) entre murundus
#'
#' @param coords matriz ou data.frame com colunas x, y (coordenadas geográficas
#'        ou métricas; se UTM, resultado em metros)
#' @return vetor de distâncias ENN para cada murundu
enn_murundus <- function(coords) {
  coords  <- as.matrix(coords[, c("x", "y")])
  n       <- nrow(coords)
  distancias <- as.matrix(dist(coords))
  diag(distancias) <- Inf
  apply(distancias, 1, min)
}

#' Calcular índice de conectividade de Hanski (1994) para manchas
#'
#' S_i = sum_{j != i} A_j * exp(-alpha * d_ij)
#'
#' @param areas vetor de áreas das manchas
#' @param distancias matriz de distâncias entre manchas (mesma unidade que 1/alpha)
#' @param alpha taxa de decaimento da dispersão (1/distância média)
#' @return vetor de índices de conectividade para cada mancha
conectividade_hanski <- function(areas, distancias, alpha) {
  n <- length(areas)
  Si <- numeric(n)
  for (i in seq_len(n)) {
    j_idx <- setdiff(seq_len(n), i)
    Si[i]  <- sum(areas[j_idx] * exp(-alpha * distancias[i, j_idx]))
  }
  Si
}

# -----------------------------------------------------------------------------
# 7. UTILITÁRIOS
# -----------------------------------------------------------------------------

#' Normalizar vetor para [0, 1]
normalizar <- function(x) (x - min(x)) / (max(x) - min(x))

#' Gerar data.frame com dados simulados de murundus para testes
#'
#' @param n_murundus número de murundus a simular
#' @param semente semente aleatória
#' @return data.frame com colunas: id, area, isolamento, riqueza, N_ind
simular_dados_murundus <- function(n_murundus = 30, semente = 42) {
  set.seed(semente)
  area      <- rlnorm(n_murundus, meanlog = 1.5, sdlog = 0.8)   # m²
  isolamento <- rexp(n_murundus, rate = 0.05)                   # metros
  z_true     <- 0.28
  c_true     <- 3.5
  riqueza    <- round(c_true * area^z_true + rnorm(n_murundus, 0, 1))
  riqueza    <- pmax(riqueza, 1L)
  N_ind      <- round(riqueza * rlnorm(n_murundus, meanlog = 2, sdlog = 0.5))
  data.frame(
    id         = seq_len(n_murundus),
    area       = area,
    isolamento = isolamento,
    riqueza    = riqueza,
    N_ind      = N_ind
  )
}
