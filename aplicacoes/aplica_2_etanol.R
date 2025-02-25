#remotes::install_github("PedroHCAlmeida/mixRegEM")
#install.packages("SemiPar")

library(SemiPar) # Pacote com os dados
library(mixRegEM) # Pacote com os algoritmos
library(dplyr) # Pacote para auxiliar manipulação dos resultados

# Carrega o conjunto de dados
data(ethanol)

# Define as variáveis X, y e R
df = ethanol
X = df[,c('E')] |>
  as.matrix()

y = df$NOx

R = df[,c('E', 'C')] |>
  as.matrix()

# Define o número de grupos a serem testados
g = 1:4

# Define distribuicoes
distribuicoes = c(
  "MixNormal",
  "MixT",
  "MixSN",
  "MixST",
  "MoENormal",
  "MoET",
  "MoESN",
  "MoEST"
)

# Roda algoritmos do cenário 1
modelos_cenario_1 = regEM(
  y,
  X,
  r = R,
  g = g,
  family = distribuicoes,
  tol = 1E-4,    # Critério de parada
  nuIgual = F,   # Define nu diferentes entre grupos
  varEqual = F,  # Define gamma diferentes entre grupos
  verbose = TRUE,# Exibe mensagens na tela
  max_iter = 100 # Define máximo de iterações
)

# Roda algoritmos do cenário 2
modelos_cenario_2 = regEM(
  y,
  X,
  r = R,
  g = g,
  family = distribuicoes,
  tol = 1E-4,     # Critério de parada
  nuIgual = F,    # Define nu diferentes entre grupos
  varEqual = T,   # Define gamma iguais entre grupos
  verbose = TRUE, # Exibe mensagens na tela
  max_iter = 100  # Define máximo de iterações
)

# Roda algoritmos do cenário 3
modelos_cenario_3 = regEM(
  y,
  X,
  r = R,
  g = g,
  family = distribuicoes[c(2,4,6,8)], # Define apenas distribuições que possuem nu
  tol = 1E-4,     # Critério de parada
  nuIgual = T,    # Define nu iguais entre grupos
  varEqual = F,   # Define gamma diferentes entre grupos
  verbose = TRUE, # Exibe mensagens na tela
  max_iter = 100  # Define máximo de iterações
)

# Roda algoritmos do cenário 4
modelos_cenario_4 = regEM(
  y,
  X,
  r = R,
  g = g,
  family = distribuicoes[c(2,4,6,8)],
  tol = 1E-4,     # Critério de parada
  nuIgual = T,    # Define nu iguais entre grupos
  varEqual = T,   # Define gamma iguais entre grupos
  verbose = TRUE, # Exibe mensagens na tela
  max_iter = 100  # Define máximo de iterações
)

# Combina modelos em uma lista
modelos_ethanol_list = list(
  modelos_cenario_1,
  modelos_cenario_2,
  modelos_cenario_3,
  modelos_cenario_4
)

# Combina resultados em uma tabela
modelos_ethanol =
  modelos_cenario_1$results |>
  mutate(cenario = 1) |>
  rbind(
    modelos_cenario_2$results |>
      mutate(cenario = 2)
  ) |>
  rbind(
    modelos_cenario_3$results |>
      mutate(cenario = 3)
  ) |>
  rbind(
    modelos_cenario_4$results |>
      mutate(cenario = 4)
  )


# Verifica qual o melhor modelo
best_config = modelos_ethanol[modelos_ethanol$BIC |> which.min(),]
# seleciona o melhor modelo
best_model = modelos_ethanol_list[[best_config$cenario]]$best[[best_config$family]][[best_config$g]]

# Exibe os resultados
best_model

