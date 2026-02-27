# Síntese Teórica: Ecologia, Conservação e Biogeografia dos Murundus

> **Denis S. Nogueira — Projeto Murundus**
>
> *"A ecologia é a biologia dos números e dos processos que os geram."*

---

## 1. Introdução: O Problema da Diversidade

A biodiversidade — a variedade de formas de vida e suas interações — é simultaneamente o objeto central da ecologia e seu maior enigma. Por que existem tantas espécies? Por que estão distribuídas da forma como estão? Como coexistem? E como persistirão diante das perturbações antrópicas crescentes?

Essas questões não admitem resposta única. Ao longo do século XX e início do XXI, diferentes tradições teóricas ofereceram respostas parciais e, muitas vezes, aparentemente contraditórias. Este documento propõe uma síntese dessas tradições, articulando seus fundamentos filosóficos e matemáticos com aplicações práticas ao estudo e conservação dos murundus do Cerrado brasileiro.

A síntese não pretende dissolver as contradições, mas torná-las produtivas: cada teoria captura dimensões reais do fenômeno ecológico, e sua integração produz uma compreensão mais rica do que qualquer uma isoladamente.

---

## 2. Fundamentos Filosóficos

### 2.1 Ontologia Ecológica: O que existe no mundo ecológico?

Antes de modelar, é preciso decidir o que existe. As teorias ecológicas diferem em suas ontologias:

- **Teoria de Nicho**: Postula a existência de entidades funcionais distintas — espécies com atributos determinísticos (tolerâncias fisiológicas, taxas vitais) que determinam onde e como vivem. O mundo ecológico é fundamentalmente estruturado por diferenças entre espécies. Raízes em Grinnell (1917) e Elton (1927), formalização por Hutchinson (1957).

- **Teoria Neutra**: Postula que as diferenças funcionais entre espécies são ecologicamente irrelevantes — os indivíduos são equivalentes per capita, e a diversidade emerge de processos estocásticos de nascimento, morte, migração e especiação. O mundo ecológico é fundamentalmente aleatório ao nível da interação individual. Raízes em Caswell (1976), formalização por Hubbell (2001).

- **Biogeografia de Ilhas**: Postula que comunidades locais são determinadas pelo equilíbrio dinâmico entre dois processos opostos — colonização (imigração) e extinção — mediados pela geometria (área, isolamento). A história e a geometria da paisagem importam tanto quanto os atributos das espécies.

- **Metacomunidades**: Síntese explícita entre processos locais (nicho, competição, predação) e regionais (dispersão, deriva estocástica). Quatro paradigmas: dinâmica de manchas, triagem de espécies, efeitos de massa, neutro.

- **Macroecologia**: Propõe que padrões emergentes em grandes escalas — distribuições de abundâncias, relações tamanho-abundância, gradientes de diversidade latitudinal — são objetos legítimos e irredutíveis aos processos locais. Propriedades emergentes de sistemas complexos.

### 2.2 Epistemologia: Como conhecemos os padrões e processos?

- **Indutivismo popperiano**: Hipóteses derivadas de teorias são testadas contra dados. A teoria neutra, por exemplo, gera predições quantitativas para a Série Logarítmica e a Série de Preston que podem ser falsificadas.

- **Inferência por múltiplos modelos** (Anderson & Burnham 2002): Dado que múltiplas teorias geram predições sobrepostas, a abordagem AIC/BIC é epistemologicamente mais honesta do que o teste de hipótese nula único.

- **Realismo estrutural**: As equações de um modelo podem ser "verdadeiras" em sua estrutura mesmo que os parâmetros sejam estimados com incerteza. A relação espécie-área S = cA^z é estruturalmente robusta mesmo que *z* varie com o contexto.

---

## 3. Biogeografia de Ilhas

### 3.1 A Teoria do Equilíbrio (MacArthur & Wilson 1967)

A teoria do equilíbrio da biogeografia de ilhas (TIB) propõe que o número de espécies em uma ilha converge para um equilíbrio dinâmico determinado pelo balanço entre taxas de imigração (I) e extinção (E):

**Pressupostos fundamentais:**
1. A taxa de imigração de novas espécies decresce com o número de espécies já presentes na ilha (espécies disponíveis no pool regional se esgotam).
2. A taxa de extinção aumenta com o número de espécies (mais espécies = mais competição, menor N por espécie).
3. Ambas as taxas dependem do tamanho da ilha (área A) e da distância ao continente/fonte (d).

**Equação de equilíbrio:**

```
dS/dt = I(S) - E(S)
```

No equilíbrio:
```
S* = (λ_max × E_max) / (λ_max + E_max)
```

onde λ_max é a taxa máxima de imigração (S = 0) e E_max é a taxa máxima de extinção (S = P, pool regional).

**Relação Espécie-Área (SAR):**
O padrão mais robusto em ecologia insular é a relação potência:

```
S = c · A^z
```

onde *c* é uma constante que depende do taxa do organismo e da região, e *z* é o expoente de escalonamento (tipicamente 0,20–0,35 para ilhas verdadeiras; 0,12–0,20 para fragmentos de habitat).

**Linearização:** log S = log c + z · log A

**Aplicação a Murundus:** Cada murundu pode ser tratado como uma ilha de habitat. A área do murundu prediz o número de espécies (plantas, artrópodes, répteis) que o habitam. O isolamento (distância ao murundu mais próximo ou à borda da formação) modula as taxas de colonização. Um gradiente de murundus de diferentes tamanhos e distâncias é um experimento natural para testar a TIB.

### 3.2 Extensões da TIB

- **Teoria da Biogeografia de Ilhas com Triagem de Nicho** (Lomolino 2000): incorpora atributos das espécies (capacidade de dispersão, tamanho corporal) para predizer quais espécies estão presentes além de quantas.
- **Modelo de Incidência** (Diamond 1975): j(S) = probabilidade de ocorrência da espécie em uma ilha de tamanho S. Funções de incidência permitem identificar espécies de interior vs. espécies generalistas.
- **Efeitos de Área vs. Heterogeneidade de Habitat** (Triantis et al. 2012): O aumento de S com A pode ser mediado por maior variedade de micro-habitats. Murundus maiores têm mais heterogeneidade topográfica e de solo.

---

## 4. Teoria de Nicho

### 4.1 O Nicho de Hutchinson (1957)

Hutchinson formalizou o nicho como um **hipervolume n-dimensional** no espaço de variáveis ambientais (temperatura, umidade, pH, recursos, etc.) dentro do qual uma espécie pode manter populações positivas:

**Nicho Fundamental (NF):** região do espaço ambiental onde a espécie pode sobreviver e reproduzir na ausência de competidores.

**Nicho Realizado (NR):** subconjunto do NF onde a espécie realmente ocorre, após a exclusão competitiva e limitação de dispersão.

**Matematicamente:**
```
NF_i = {x ∈ ℝ^n : r_i(x) > 0}
```

onde r_i(x) é a taxa de crescimento per capita da espécie i no ponto ambiental x.

**Sobreposição de Nicho (Pianka 1973):**
```
O_ij = Σ_k (p_ik · p_jk) / sqrt(Σ_k p_ik² · Σ_k p_jk²)
```

onde p_ik é a proporção do recurso k utilizada pela espécie i.

**Amplitude de Nicho (Levins 1968):**
```
B_i = 1 / Σ_k p_ik²
```

B_i varia de 1 (especialista absoluto) a K (generalista perfeito entre K recursos).

**Aplicação a Murundus:** As condições de solo, umidade, luz e temperatura variam sistematicamente entre a superfície dos murundus (mais seco, mais exposto, solo mais compactado) e a matriz circundante (mais úmido, anegado sazonalmente). Espécies com nichos fundamentais distintos são filtradas diferencialmente entre os dois micro-habitats.

### 4.2 Modelos de Distribuição de Espécies (SDMs)

SDMs estimam a distribuição geográfica do nicho realizado:

- **MaxEnt** (Phillips et al. 2006): Maximiza entropia sujeita a restrições impostas pelas ocorrências observadas. Equivalente a modelo de regressão logística regularizado com dados de presença-background.
- **BioClim/GARP**: Abordagens de envelope climático.
- **Modelos de nicho mecanísticos** (Kearney & Porter 2009): Integradores de fisiologia e clima, mais robustos para extrapolação.

**Em murundus:** SDMs podem identificar quais atributos do murundu (área, altura, distância a cursos d'água, composição de solo) são preditores da ocorrência de espécies-chave, e projetar como o desmatamento do entorno altera a disponibilidade de nichos.

### 4.3 Filtragem de Habitat e Montagem de Comunidades

A teoria de nicho prevê que as comunidades são estruturadas por filtragem ambiental (espécies com nichos incompatíveis são excluídas) e por exclusão competitiva (espécies com nichos muito similares não coexistem — princípio de exclusão competitiva de Gause):

```
Coexistência estável ↔ diferenciação de nicho suficiente
```

A teoria moderna de coexistência (Chesson 2000) formaliza isso em:
- **Efeitos de estabilização** (nicho differentiation): promovem coexistência aumentando a intra-específica relativa à interespecífica competição.
- **Efeitos de equalização** (fitness differentiation): diferenças de aptidão (fitness) entre espécies que, sem estabilização, levariam à exclusão.

---

## 5. Teoria Neutra da Biodiversidade e Biogeografia (UNTB)

### 5.1 O Postulado de Equivalência per Capita (Hubbell 2001)

A teoria neutra de Hubbell parte de um pressuposto radicalmente simplificador: todos os indivíduos de todas as espécies são ecologicamente equivalentes per capita — têm a mesma taxa de nascimento, morte, imigração e especiação. A diversidade emerge de **deriva ecológica estocástica** (análoga à deriva genética).

**Comunidade Local de Tamanho J:**
A dinâmica é um processo de morte-nascimento com imigração. A cada passo de tempo:
1. Um indivíduo morre aleatoriamente (prob. 1/J por indivíduo).
2. Com probabilidade (1 − m), o slot é preenchido por descendente local (escolha proporcional às abundâncias locais).
3. Com probabilidade m (taxa de imigração), o slot é preenchido por imigrante da metacomunidade.

**Parâmetros:**
- **θ (número fundamental de biodiversidade):** θ = 2·J_M·ν, onde J_M é o tamanho da metacomunidade e ν é a taxa de especiação per capita. Controla a diversidade da metacomunidade.
- **m (taxa de imigração):** Fração de mortes preenchidas por imigrantes. Controla o acoplamento comunidade local ↔ metacomunidade.
- **I (número de imigrantes):** I = m·(J − 1)/(1 − m).

**Distribuição de Abundâncias de Espécies (SAD) Neutra:**

A SAD da metacomunidade é a **série logarítmica de Fisher**:
```
φ(n) = θ/n · x^n / Γ(n+1)   [approx para n grande]
```

A SAD da comunidade local é a **distribuição de Ewens** (Ewens sampling formula):
```
P(n_1, ..., n_S | J, θ) = (J! · θ^S) / (θ^(J)↑ · Π n_i · Π m_i!)
```

onde θ^(J)↑ = θ(θ+1)...(θ+J−1) é o fatorial ascendente de Pochhammer.

**Riqueza Esperada:**
```
E[S] = θ · Σ_{k=1}^{J} 1/(θ + k - 1) ≈ θ · ln(1 + J/θ)
```

### 5.2 UNTB e Biogeografia de Ilhas

A UNTB unifica a TIB dentro de um framework estocástico: a taxa de imigração m da UNTB é o análogo da taxa de colonização I da TIB. As predições diferem quantitativamente (a UNTB prediz SADs; a TIB prediz apenas S*), mas convergem qualitativamente.

**Aplicação a Murundus:** Murundus com menor conectividade (maior isolamento) devem ter menor m. A estimativa de m por máxima verossimilhança (Etienne 2005) a partir de dados de composição de espécies de múltiplos murundus permite testar se o isolamento geográfico afeta o acoplamento regional de forma quantitativamente consistente com a UNTB.

### 5.3 Limites da Neutralidade

- A UNTB gera SADs que se ajustam bem a muitos conjuntos de dados, mas os mecanismos subjacentes diferem dos mecanismos reais (a equivalência per capita é sabidamente falsa para a maioria das comunidades).
- Mcgill et al. (2006): múltiplos modelos não-neutros geram SADs indistinguíveis do modelo neutro — a SAD não é diagnóstico suficiente.
- O valor da UNTB está em sua função como **modelo nulo**: desvios da neutralidade indicam onde os processos de nicho são importantes.

---

## 6. Ecologia de Metacomunidades

### 6.1 O Conceito de Metacomunidade (Leibold et al. 2004)

Uma metacomunidade é um **conjunto de comunidades locais conectadas por dispersão**. O conceito integra dinâmica local (competição, predação, filtragem ambiental) com dinâmica regional (dispersão, colonização-extinção) em quatro paradigmas:

| Paradigma | Dispersão | Nicho/Equivalência | Dinâmica dominante |
|---|---|---|---|
| **Dinâmica de Manchas** | Alta | Espécies equivalentes | Colonização-extinção estocástica |
| **Triagem de Espécies** | Alta | Diferenciação de nicho | Filtragem ambiental determinística |
| **Efeitos de Massa** | Alta | Diferenciação de nicho | Imigração sustenta populações sink |
| **Neutro** | Variável | Equivalência per capita | Deriva estocástica + imigração |

**A dispersão medeia a transição entre paradigmas:** Com dispersão muito baixa, cada comunidade evolui independentemente (isolamento); com dispersão muito alta, a metacomunidade converge para uma comunidade única bem misturada; com dispersão intermediária, a variação espacial no ambiente e nas condições locais molda a estrutura da metacomunidade.

### 6.2 Variação de Leibold: Partição de Diversidade Beta

A diversidade beta (β-diversidade) é a variação na composição de espécies entre comunidades. Leibold & Mikkelson (2002) propõem a partição:

```
β_total = β_nestedness + β_turnover
```

- **Turnover (reposição):** Espécies de uma comunidade substituem espécies de outra.
- **Nestedness (aninhamento):** A composição da comunidade mais pobre é um subconjunto da mais rica.

**Aninhamento e Murundus:** Se a composição de espécies em murundus menores é sistematicamente um subconjunto da composição em murundus maiores, a estrutura é aninhada — consistente com extinção aleatória e colonização preferencial dos maiores. A partição de β usando o pacote `betapart` (Baselga 2010) permite testar esta hipótese.

### 6.3 Análise de Variação de Comunidade (RDA / db-RDA)

A estrutura da metacomunidade pode ser decomposta em:
1. **Variação ambiental pura** (condições locais do murundu)
2. **Variação espacial pura** (estrutura geográfica da paisagem)
3. **Variação ambiental e espacial confundida**
4. **Variação não explicada** (estocasticidade)

Esta partição é feita via **db-RDA** (Redundancy Analysis baseada em distâncias) com variáveis de Moran Eigenvector Maps (MEM) para capturar estrutura espacial (Dray et al. 2006):

```
β = [Env | Spa | Env∩Spa | Resid]
```

Quando a variação espacial pura é alta (após controle ambiental), inferimos que limitação de dispersão é importante (consistente com neutralidade ou dinâmica de manchas). Quando a variação ambiental pura domina, inferimos triagem de nicho.

---

## 7. Ecologia de Paisagens

### 7.1 Estrutura da Paisagem e Fluxos Ecológicos

A ecologia de paisagens (Forman & Godron 1986; Turner 1989) analisa a **estrutura espacial da paisagem** — o arranjo de manchas (patches), corredores e matriz — e seus efeitos sobre processos ecológicos (fluxo de energia, movimento de animais, propagação de perturbações).

**Métricas de Paisagem (FRAGSTATS):**
- **PLAND**: Percentagem da paisagem coberta por cada classe de uso
- **LPI**: Maior manchas individual (Large Patch Index)
- **PARA**: Razão perímetro:área
- **ENN**: Distância ao vizinho mais próximo (Euclidean Nearest Neighbor)
- **CONNECT**: Conectividade de manchas dentro de uma distância limiar
- **COHESION**: Coesão espacial de manchas

Para murundus: a densidade de murundus, o espaçamento médio entre eles, e o contraste entre o murundu e a matriz (campos alagados vs. campos secos vs. pastagem cultivada) são métricas fundamentais.

### 7.2 Teoria do Gradiente de Paisagem

Mais do que dividir a paisagem em classes discretas, a abordagem de gradiente (McGarigal & Cushman 2005) trata a cobertura da terra como variável contínua. Superfícies de resistência ao movimento dos organismos substituem a classificação binária habitat/não-habitat.

**Circuitscape:** Modela a conectividade funcional como fluxo elétrico em circuitos de resistência. Prediz a intensidade de migração entre manchas com base na teoria de circuitos (McRae et al. 2008):

```
I_ij = V_ij / R_ij
```

onde V_ij é a diferença de "voltagem" (densidade de imigrantes) e R_ij é a resistência efetiva entre as manchas i e j.

**Aplicação:** A resistência entre murundus depende da matriz circundante — pastagens manejadas têm resistência diferente de campos nativos, e esta diferença afeta a taxa efetiva de imigração (m na UNTB).

### 7.3 Limiares de Percolação e Extinções em Cascata

A teoria de percolação prevê que, quando a cobertura do habitat cai abaixo de um limiar crítico (~20–30% em modelos de grade), a conectividade da paisagem colapsa abruptamente — uma transição de fase. Abaixo do limiar, ilhas de habitat tornam-se funcionalmente isoladas e as taxas de extinção aumentam super-linearly:

```
p_c ≈ 0,5927 (grade quadrada, percolação de sítios)
```

Para murundus em paisagens desmatadas do Cerrado, este limiar tem implicações diretas para a viabilidade a longo prazo das populações que dependem de conectividade entre murundus.

---

## 8. Macroecologia

### 8.1 Distribuições de Abundância de Espécies (SAD)

A SAD descreve a frequência com que espécies de diferentes abundâncias ocorrem numa comunidade. Modelos concorrentes:

| Modelo | Fórmula | Mecanismo proposto |
|---|---|---|
| **Série Logarítmica** (Fisher 1943) | φ(n) = α·x^n/n | Chegada aleatória de indivíduos de um pool infinito |
| **Log-Normal** (Preston 1948) | S(R) = S_0·e^{-a²R²} | Produto de muitas variáveis independentes (TLC) |
| **Neutro** (Hubbell 2001) | Zero-Sum Multinomial | Deriva estocástica + especiação |
| **Broken Stick** (MacArthur 1957) | E[n_i] = (N/S)·Σ_{k=i}^{S} 1/k | Divisão aleatória de nicho |

**Comparação de modelos:** A log-normal é consistente com comunidades onde muitas espécies têm abundâncias intermediárias (padrão de Preston). A série logarítmica caracteriza amostras pequenas de comunidades diversas. O modelo neutro produz SADs intermediárias entre os dois extremos.

### 8.2 Relações de Escalonamento Macroecológico

**Lei de Damuth (1987):** Densidade populacional (D) decresce com massa corporal (M):
```
D ~ M^{-3/4}
```

Deriva da teoria metabólica da ecologia (Brown et al. 2004): o metabolismo basal B ~ M^{3/4}, e a energia total disponível divide-se proporcionalmente ao metabolismo.

**Teoria Metabólica da Ecologia (MTE):**
```
B = b_0 · M^{3/4} · e^{-E/(k·T)}
```

onde E é a energia de ativação (~0,65 eV), k é a constante de Boltzmann e T é a temperatura (Kelvin). Prediz taxas de crescimento populacional, tempo de geração, diversidade em gradientes de temperatura e produtividade.

**Relação Riqueza-Energia:** S ~ Produção Primária Líquida (NPP). Em murundus: NPP local do murundu (maior biomassa vegetal, solo mais fértil) vs. NPP da matriz deve predizer o diferencial de diversidade.

### 8.3 Gradientes de Diversidade

- **Gradiente latitudinal:** Aumento de riqueza em direção aos trópicos. Hipóteses: maior energia solar, maior área efetiva, maior diversificação evolutiva, maior estabilidade climática histórica.
- **Gradiente altitudinal:** Padrão humped (mid-domain effect) ou monotônico decrescente conforme o grupo taxonômico.
- **Gradiente de produtividade:** Relação hump-shaped entre produtividade e diversidade (mais divergente; ver Mittelbach et al. 2001).

---

## 9. Síntese Integrativa: Murundus como Sistema Modelo

### 9.1 Hierarquia de Processos e Escalas

A comunidade de um murundu individual é moldada por processos que operam em múltiplas escalas:

```
Escala Geográfica (continental)
  └── Pool Regional de Espécies (biogeografia histórica)
        └── Metacomunidade (dispersão + filtros regionais)
              └── Paisagem Local (conectividade, matrix quality)
                    └── Murundu Individual (filtragem ambiental, competição, deriva)
```

Cada escala impõe restrições ao conjunto de espécies que pode estar presente na escala imediatamente inferior. A teoria de metacomunidades opera na interface entre as três escalas intermediárias.

### 9.2 Diagrama de Integração Teórica

```
                    MACROECOLOGIA
                   (SADs, scaling)
                        |
          +-------------+-------------+
          |                           |
    BIOGEOGRAFIA               TEORIA NEUTRA
    DE ILHAS                   (θ, m, deriva)
    (S*, SAR)                       |
          |                         |
          +-------+    +------------+
                  |    |
            METACOMUNIDADE
            (patch dynamics,
             species sorting,
             mass effects)
                  |
          +-------+-------+
          |               |
    ECOLOGIA DE       TEORIA DE
    PAISAGENS         NICHO
    (conectividade,   (hipervolume,
     fragmentação)    coexistência)
                  |
           CONSERVAÇÃO
           (murundus do
            Cerrado)
```

### 9.3 Perguntas Integrativas para o Projeto Murundus

1. **A relação espécie-área nos murundus tem o expoente z típico de ilhas verdadeiras (z ≈ 0,30) ou de fragmentos de habitat (z ≈ 0,15)?** A resposta indica se os murundus funcionam como ilhas verdadeiras ou como manchas dentro de uma matriz permeável.

2. **A estrutura de β-diversidade é dominada por turnover ou por nestedness?** Se nestedness domina, os murundus maiores são a fonte de colonização para os menores — implicação para priorização de conservação.

3. **O ajuste do modelo neutro (UNTB) aos dados de abundância melhora ou piora com o isolamento do murundu?** Se m estimado decresce com o isolamento, os murundus mais isolados estão mais sob dominância de deriva — mais vulneráveis a extinções locais estocásticas.

4. **A variação de composição entre murundus é melhor explicada por variáveis ambientais locais (suporte para triagem de nicho) ou por estrutura espacial (suporte para limitação de dispersão/neutralidade)?** A partição de variância via db-RDA responde diretamente.

5. **A densidade de murundus na paisagem está próxima do limiar de percolação?** Se sim, pequenas reduções adicionais da cobertura podem causar colapso de conectividade e extinções em cascata.

6. **As SADs dos murundus seguem log-normal (nicho estruturado) ou série logarítmica (comunidades saturadas)?** A forma da SAD muda com o tamanho do murundu, consistente com a transição entre paradigmas metacomunidade?

---

## 10. Implicações para Conservação

### 10.1 Priorização de Áreas

- **Murundus maiores:** maior riqueza de espécies (lei de potência), maior estabilidade populacional, menor risco de extinção estocástica. Prioridade máxima.
- **Murundus conectores:** murundus que estão na rota de fluxo efetivo entre grupos maiores (identificados por Circuitscape). Críticos para manutenção da metacomunidade.
- **Murundus com alta diversidade beta local:** murundus com composição muito diferente de seus vizinhos — contribuem de forma desproporcional para a diversidade regional (γ-diversidade).

### 10.2 Estratégias de Restauração

- **Stepping-stones:** Restaurar ou proteger murundus intermediários para reduzir a distância efetiva de dispersão entre núcleos maiores.
- **Corredores:** Manter ou recuperar vegetação nativa na matriz entre murundus reduz a resistência ao movimento e aumenta m.
- **Gestão da qualidade da matriz:** A resistência de pastagens manejadas é menor que a de lavouras — o manejo extensivo da pecuária pode ser compatível com a manutenção de conectividade funcional entre murundus.

### 10.3 Monitoramento Baseado em Teoria

- **Indicadores de neutralidade vs. nicho:** Se a proporção de espécies dominantes aumenta ao longo do tempo em murundus isolados (deriva), isso sinaliza degradação funcional antes de extinções observáveis.
- **Curvas de rarefação e z:** Monitorar z (expoente SAR) ao longo do tempo — um aumento de z indica que os murundus estão se tornando mais insulares (matrix mais hostil).
- **Diversidade funcional (FD):** Redução de FD antes de redução de riqueza taxonômica — sinal precoce de filtragem ambiental intensificada ou homogeneização biótica.

---

## 11. Considerações Filosóficas Finais

A tensão entre determinismo (nicho) e estocasticidade (neutra) na ecologia reflete uma tensão filosófica mais profunda entre o **mecanicismo** (todo padrão tem causa determinística discernível) e o **probabilismo** (padrões emergem de processos fundamentalmente estocásticos). A síntese contemporânea não resolve esta tensão — ela a reformula como questão empírica: **em que escala, e para quais sistemas, cada processo domina?**

Para os murundus do Cerrado, a resposta provavelmente é: processos de nicho dominam na escala local do murundu individual (filtragem ambiental entre a superfície e a matriz), processos de deriva e dispersão dominam na escala da metacomunidade de murundus próximos, e padrões macroecológicos emergem da interação entre ambos ao longo de gradientes climáticos e históricos do Cerrado.

Esta síntese multi-escalar, multi-processo, é o que torna os murundus — e a ecologia — fascinantes.

---

## Referências

- Anderson, D.R. & Burnham, K.P. (2002). *Model Selection and Multimodel Inference*. Springer.
- Baselga, A. (2010). Partitioning the turnover and nestedness components of beta diversity. *Global Ecology and Biogeography*, 19, 134–143.
- Brown, J.H. et al. (2004). Toward a metabolic theory of ecology. *Ecology*, 85, 1771–1789.
- Caswell, H. (1976). Community structure: a neutral model analysis. *Ecological Monographs*, 46, 327–354.
- Chesson, P. (2000). Mechanisms of maintenance of species diversity. *Annual Review of Ecology and Systematics*, 31, 343–366.
- Diamond, J.M. (1975). Assembly of species communities. In: Cody & Diamond (eds.), *Ecology and Evolution of Communities*, pp. 342–444.
- Dray, S. et al. (2006). Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). *Ecological Modelling*, 196, 483–493.
- Elton, C.S. (1927). *Animal Ecology*. Sidgwick & Jackson.
- Etienne, R.S. (2005). A new sampling formula for neutral biodiversity. *Ecology Letters*, 8, 253–260.
- Fisher, R.A., Corbet, A.S. & Williams, C.B. (1943). The relation between the number of species and the number of individuals. *Journal of Animal Ecology*, 12, 42–58.
- Forman, R.T.T. & Godron, M. (1986). *Landscape Ecology*. Wiley.
- Grinnell, J. (1917). The niche-relationships of the California thrasher. *The Auk*, 34, 427–433.
- Hubbell, S.P. (2001). *The Unified Neutral Theory of Biodiversity and Biogeography*. Princeton University Press.
- Hutchinson, G.E. (1957). Concluding remarks. *Cold Spring Harbor Symposia on Quantitative Biology*, 22, 415–427.
- Kearney, M. & Porter, W. (2009). Mechanistic niche modelling: combining physiological and spatial data to predict species ranges. *Ecology Letters*, 12, 334–350.
- Leibold, M.A. et al. (2004). The metacommunity concept: a framework for multi-scale community ecology. *Ecology Letters*, 7, 601–613.
- Leibold, M.A. & Mikkelson, G.M. (2002). Coherence, species turnover, and boundary clumping: elements of meta-community structure. *Oikos*, 97, 237–250.
- Levins, R. (1968). *Evolution in Changing Environments*. Princeton University Press.
- Lomolino, M.V. (2000). Ecology's most general, yet protean pattern: the species-area relationship. *Journal of Biogeography*, 27, 17–26.
- MacArthur, R.H. (1957). On the relative abundance of bird species. *Proceedings of the National Academy of Sciences*, 43, 293–295.
- MacArthur, R.H. & Wilson, E.O. (1967). *The Theory of Island Biogeography*. Princeton University Press.
- McGarigal, K. & Cushman, S.A. (2005). The gradient concept of landscape structure. In: Wiens & Moss (eds.), *Issues and Perspectives in Landscape Ecology*, pp. 112–119.
- McGill, B.J. et al. (2006). Empirical evaluation of neutral theory. *Ecology*, 87, 1411–1423.
- McRae, B.H. et al. (2008). Using circuit theory to model connectivity in ecology, evolution, and conservation. *Ecology*, 89, 2712–2724.
- Mittelbach, G.G. et al. (2001). What is the observed relationship between species richness and productivity? *Ecology*, 82, 2381–2396.
- Oliveira-Filho, A.T. (1992). The vegetation of Brazilian murundus. *Journal of Tropical Ecology*, 8, 465–486.
- Phillips, S.J., Anderson, R.P. & Schapire, R.E. (2006). Maximum entropy modeling of species geographic distributions. *Ecological Modelling*, 190, 231–259.
- Pianka, E.R. (1973). The structure of lizard communities. *Annual Review of Ecology and Systematics*, 4, 53–74.
- Preston, F.W. (1948). The commonness and rarity of species. *Ecology*, 29, 254–283.
- Tilman, D. (1982). *Resource Competition and Community Structure*. Princeton University Press.
- Triantis, K.A. et al. (2012). The island species-area relationship: biology and statistics. *Journal of Biogeography*, 39, 215–231.
- Turner, M.G. (1989). Landscape ecology: the effect of pattern on process. *Annual Review of Ecology and Systematics*, 20, 171–197.
