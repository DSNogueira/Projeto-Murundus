# Resultados do Projeto Murundus

Esta pasta armazena os outputs gerados pelos scripts de análise.

## Estrutura Prevista

```
resultados/
├── figuras/
│   ├── sar_murundus.png               # Relação espécie-área
│   ├── extincao_sar_murundus.png      # Extinção esperada por perda de área
│   ├── nmds_comunidades.png           # Ordenação NMDS das comunidades
│   ├── beta_diversidade.png           # Partição de beta-diversidade
│   ├── sad_rank_abundance.png         # SAD (rank-abundance)
│   ├── sad_preston.png                # SAD (octaves de Preston)
│   ├── damuth_densidade_massa.png     # Lei de Damuth
│   ├── mte_temperatura.png            # MTE: riqueza × temperatura
│   ├── nicho_overlap.png              # Sobreposição de nicho (heatmap)
│   └── ses_filtragem.png              # SES de variância de traço
└── tabelas/
    ├── sar_parametros.csv             # Parâmetros da SAR (c, z, R², AIC)
    ├── beta_particao.csv              # Partição de β (turnover, nestedness)
    ├── varpart_resultados.csv         # Partição de variância (ambiental vs. espacial)
    └── neutro_parametros.csv          # Parâmetros neutros (θ, m) por murundu
```

## Reproducibilidade

Todos os resultados são totalmente reproduzíveis a partir dos scripts em `scripts/R/`
com os dados em `dados/`. Para reproduzir, execute os scripts na seguinte ordem:

1. `funcoes_ecologicas.R` (funções base — carregado automaticamente)
2. `biogeografia_ilhas_murundus.R`
3. `teoria_neutra_murundus.R`
4. `metacomunidade_murundus.R`
5. `nicho_murundus.R`
6. `macroecologia_murundus.R`

Para salvar as figuras, descomente as chamadas a `ggsave()` ao final de cada script.
