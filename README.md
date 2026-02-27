# Projeto Murundus

Este repositório congrega scripts, dados, resultados e publicações do **Projeto Murundus** — uma investigação ecológica sobre a biodiversidade associada a murundus (termiteiros ou campos de murundus) no Cerrado brasileiro, integrando fundamentos de biogeografia de ilhas, teoria neutra, teoria de nicho, ecologia de metacomunidades, ecologia de paisagens e macroecologia numa síntese teórica e prática para a conservação.

---

## O que são murundus?

Murundus (também chamados *campos de murundus*, *termiteiros*, *covoais* ou *tussocks*) são estruturas físicas proeminentes — montículos de solo argiloso e/ou ninhos de cupins — que se elevam sobre a superfície plana dos campos e veredas do Cerrado. Estas ilhas de micro-habitat diferenciado concentram espécies de plantas, invertebrados, répteis, anfíbios e aves em densidades contrastantes com a matriz circundante, tornando-as objetos naturais de estudo para as teorias de biogeografia de ilhas, metacomunidades e ecologia de paisagens.

---

## Estrutura do Repositório

```
Projeto-Murundus/
├── README.md
├── docs/
│   └── sintese_teorica.md          # Síntese filosófica e matemática das teorias ecológicas
├── scripts/
│   └── R/
│       ├── funcoes_ecologicas.R    # Funções utilitárias: SAR, SAD, modelos neutros, nicho
│       ├── biogeografia_ilhas_murundus.R   # Relações espécie-área e dinâmica de colonização/extinção
│       ├── metacomunidade_murundus.R       # Análises de metacomunidade (RDA, NMDS, variação de Leibold)
│       ├── teoria_neutra_murundus.R        # Estimativa de θ e m; ajuste de RSA neutro
│       ├── nicho_murundus.R                # Hipervolume de Hutchinson, sobreposição e breadth
│       └── macroecologia_murundus.R        # SADs, relações de escalonamento, lei de Damuth
├── dados/
│   └── README.md                   # Descrição dos conjuntos de dados (a serem adicionados)
└── resultados/
    └── README.md                   # Descrição dos resultados gerados pelos scripts
```

---

## Síntese Teórica

O documento [`docs/sintese_teorica.md`](docs/sintese_teorica.md) apresenta uma síntese integrativa das principais teorias ecológicas relevantes para o estudo dos murundus, articulando:

- **Biogeografia de Ilhas** (MacArthur & Wilson 1967)
- **Teoria de Nicho** (Hutchinson 1957; Grinnell 1917)
- **Teoria Neutra da Biodiversidade e Biogeografia** (Hubbell 2001)
- **Ecologia de Metacomunidades** (Leibold et al. 2004)
- **Ecologia de Paisagens** (Forman & Godron 1986; Turner 1989)
- **Macroecologia** (Brown & Maurer 1989; Gaston & Blackburn 2000)

Cada teoria é apresentada com seus fundamentos filosóficos (ontológicos e epistemológicos), fundamentos matemáticos centrais, e suas implicações práticas para o monitoramento e conservação de murundus no Cerrado.

---

## Scripts R

Os scripts estão organizados por tema teórico e são executáveis de forma independente, embora compartilhem funções utilitárias definidas em `funcoes_ecologicas.R`.

### Dependências principais

```r
install.packages(c(
  "vegan",        # Ecologia de comunidades: diversidade, ordination, PERMANOVA
  "iNEXT",        # Rarefação e extrapolação de curvas de acumulação
  "untb",         # Teoria Neutra Unificada da Biodiversidade e Biogeografia
  "hypervolume",  # Hipervolume de Hutchinson multi-dimensional
  "betapart",     # Partição da diversidade beta
  "MuMIn",        # Seleção de modelos (AIC, BIC)
  "ggplot2",      # Visualização
  "dplyr",        # Manipulação de dados
  "tidyr"         # Transformação de dados
))
```

---

## Contexto Científico

Os murundus do Cerrado são análogos ecológicos funcionais às ilhas oceânicas no sentido de MacArthur & Wilson: cada murundu representa um fragmento de habitat rodeado por uma matriz de micro-habitat distinto. A variação em área, isolamento, altura e composição de solo entre murundus cria um gradiente natural que permite testar predições simultâneas de:

1. A relação espécie-área (S = cA^z)
2. O equilíbrio entre colonização e extinção
3. A partição entre processos determinísticos (nicho) e estocásticos (deriva neutra)
4. A estrutura de metacomunidades ao longo de gradientes de conectividade
5. Padrões macroecológicos de distribuição de abundâncias

Esta multiplicidade de ângulos analíticos faz dos murundus um sistema modelo privilegiado para avançar a síntese teórica em ecologia.

---

## Referências Chave

- MacArthur, R.H. & Wilson, E.O. (1967). *The Theory of Island Biogeography*. Princeton University Press.
- Hubbell, S.P. (2001). *The Unified Neutral Theory of Biodiversity and Biogeography*. Princeton University Press.
- Hutchinson, G.E. (1957). Concluding remarks. *Cold Spring Harbor Symposia on Quantitative Biology*, 22, 415–427.
- Leibold, M.A. et al. (2004). The metacommunity concept: a framework for multi-scale community ecology. *Ecology Letters*, 7, 601–613.
- Tilman, D. (1982). *Resource Competition and Community Structure*. Princeton University Press.
- Brown, J.H. & Maurer, B.A. (1989). Macroecology: the division of food and space among species on continents. *Science*, 243, 1145–1150.
- Forman, R.T.T. & Godron, M. (1986). *Landscape Ecology*. Wiley.
- Oliveira-Filho, A.T. (1992). The vegetation of Brazilian murundus — the island-effect on the plant community. *Journal of Tropical Ecology*, 8, 465–486.

---

## Licença

MIT © Denis S. Nogueira

