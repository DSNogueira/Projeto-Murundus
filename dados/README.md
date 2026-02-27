# Dados do Projeto Murundus

Esta pasta contém (ou conterá) os conjuntos de dados brutos e processados do Projeto Murundus.

## Estrutura Prevista

```
dados/
├── campo/
│   ├── presenca_ausencia_murundus.csv   # Matriz P/A espécies × murundus
│   ├── abundancia_murundus.csv          # Abundância de cada espécie por murundu
│   ├── atributos_murundus.csv           # Área, altura, isolamento, coord. UTM
│   └── variaveis_ambientais.csv         # Umidade, temperatura, pH, %C, granulometria
├── teledeteccao/
│   ├── ndvi_murundus.csv                # NDVI por murundu (Sentinel-2)
│   └── cobertura_paisagem.csv           # Cobertura da terra na vizinhança (buffer 500m)
└── tracos_funcionais/
    └── tracos_especies.csv              # Traços funcionais das espécies registradas
```

## Formato Padrão

Todos os arquivos CSV seguem o padrão:
- Codificação: UTF-8
- Separador: vírgula (`,`)
- Cabeçalho: primeira linha
- Dados ausentes: `NA`
- Coordenadas: SIRGAS 2000 / UTM zona 22S (EPSG:31982)

## Metadados

Cada arquivo de dados terá um arquivo `.yml` com metadados correspondente:
- Responsável pela coleta
- Data(s) de coleta
- Localidade e coordenadas da área de estudo
- Protocolo de amostragem
- Unidades das variáveis
