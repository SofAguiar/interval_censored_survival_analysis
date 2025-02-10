# Exemplo de Implementação de 
# Análise de Sobrevivência com Censura Intervalar
# para Grandes Amostras

# Sofia Moreira de Aguiar - Estatística | UFMG
# Prof. Orientador Enrico Colosimo - DEST | UFMG

# Carregando os pacotes necessários

library(tidyverse)   # Pacote para manipulação e visualização de dados
library(survival)    # Pacote para análise de sobrevivência
library(icenReg)     # Pacote para análise de sobrevivência com censura intervalar

# Criando dados sintéticos com as mesmas características dos dados originais

# Estes dados não têm valor analítico, apenas possuem o objetivo de 
# exemplificar a aplicação da análise de sobrevivência

set.seed(123) # Para garantir a reprodutibilidade do exemplo

n <- 653 # Número de indivíduos na base de dados

# Número de observações por indivíduo, variando entre 2 e 5 observações
obs_por_individuo <- sample(2:5, n, replace = TRUE)

# Expandindo os IDs para refletir múltiplas observações por indivíduo
dados <- data.frame(
  ID = rep(1:n, times = obs_por_individuo)
)

# Gerando outras variáveis com dependência entre sexo e peso ao nascer
dados <- dados %>%
  mutate(
    # Variável categórica exemplo (var_cat), representando sexo ou outro fator binário
    var_cat = rep(sample(c(0, 1), n, replace = TRUE), times = obs_por_individuo),
    
    # Variável contínua exemplo (var_cont), peso ao nascer com dependência de var_cat
    var_cont = unlist(lapply(var_cat, function(s) {
      if (s == 1) {
        sample(5500:8000, 1)  # Peso para 'var_cat' igual a 1 (por exemplo, masculino)
      } else {
        sample(1800:3500, 1)  # Peso para 'var_cat' igual a 0 (por exemplo, feminino)
      }
    })),
    
    # Idade atual dos indivíduos, gerada aleatoriamente entre 0 e 10 anos
    idade_atual = unlist(lapply(obs_por_individuo, function(x) runif(x, min = 0, max = 10))),
    
    # Risco de dupla carga (para exemplificar) com base em 'var_cat' e 'var_cont'
    risco_dupla_carga = 0.8 + 7 * var_cat + 0.015 * (var_cont - 3000)
  )

# Gerando a variável 'dupla_carga' com uma distribuição binomial, representando a presença ou não da condição
dados$dupla_carga <- rbinom(nrow(dados), 1, prob = 0.07)

# Verificando a prevalência da condição de dupla carga
cat("Prevalência de dupla carga:", mean(dados$dupla_carga) * 100, "%\n")


# Transformação dos dados para o formato necessário para análise de sobrevivência

# O objetivo aqui é transformar os dados para a estrutura que a análise de sobrevivência exige:
# A coluna 'left' representa o início do intervalo de censura e 'right' o final do intervalo.
# A variável 'cens' indica se o evento de interesse (a presença da dupla carga) foi observado ou censurado.

df_surv = dados %>% 
  group_by(ID) %>% 
  summarise(
    left = ifelse(
      any(dupla_carga == 1),
      max(idade_atual[idade_atual < min(idade_atual[dupla_carga == 1], na.rm = T)
                      & dupla_carga == 0], na.rm = T),
      max(idade_atual) # Caso não tenha ocorrência do evento, a maior idade é registrada
    ),
    right = ifelse(any(dupla_carga == 1),
                   min(idade_atual[dupla_carga == 1], na.rm = T),
                   NA), # Caso o evento não tenha ocorrido, substitui por NA
    cens = ifelse(any(dupla_carga == 1), 1, 0),  # Indica se o evento foi observado (1) ou censurado (0)
    var_cat = first(var_cat),   # Variáveis de covariáveis para análise
    var_cont = first(var_cont)
  ) %>% 
  ungroup() %>% 
  mutate(left = ifelse(left < 0, 0, left))  # Garante que 'left' não seja negativo

# Estimativa de Turnbull para a variável categórica
# O método de Turnbull estima a função de sobrevivência ajustada a partir da censura intervalar.

# Funções escritas por Colosimo e Giolo em "Análise de Sobrevivência Aplicada" (2021)

#################################################################

# Função para criar o vetor de tempos (tau) a partir dos intervalos censurados
cria.tau <- function(data, digits = 4) {
  l <- data$left
  r <- data$right
  
  # Arredonda os valores para evitar pequenas diferenças numéricas
  l <- round(l, digits = digits)
  r <- round(r, digits = digits)
  
  tau <- sort(unique(c(l, r[is.finite(r)])))  # Combina os tempos de censura
  return(tau)
}

# Função para inicializar a função de sobrevivência com base no vetor tau
S.ini <- function(tau){
  m <- length(tau)
  ekm <- survfit(Surv(tau[1:m-1], rep(1, m-1)) ~ 1)  # Estimativa de Kaplan-Meier inicial
  So <- c(1, ekm$surv)
  p <- -diff(So)  # Probabilidade de sobrevivência
  return(p)
}

# Função para construir a matriz A de intervalos de censura
cria.A <- function(data, tau){
  tau12 <- cbind(tau[-length(tau)], tau[-1])  # Cria os intervalos [tau[i], tau[i+1]]
  interv <- function(x, inf, sup) ifelse(x[1] >= inf & x[2] <= sup, 1, 0)
  A <- apply(tau12, 1, interv, inf = data$left, sup = data$right)  # Matriz de censura
  id.lin.zero <- which(apply(A == 0, 1, all))  # Filtra intervalos não observados
  if (length(id.lin.zero) > 0) A <- A[-id.lin.zero, ]
  return(A)
}

# Função de Turnbull, que realiza a estimação iterativa da função de sobrevivência
Turnbull <- function(p, A, data, eps = 1e-3, iter.max = 200, verbose = FALSE){
  n <- nrow(A)
  m <- ncol(A)
  Q <- matrix(1, m)
  iter <- 0
  repeat {
    iter <- iter + 1
    diff <- (Q - p)
    maxdiff <- max(abs(as.vector(diff)))  # Diferença máxima entre iterações
    if (verbose)
      print(maxdiff)
    if (maxdiff < eps | iter >= iter.max)
      break
    Q <- p
    C <- A %*% p
    p <- p * ((t(A) %*% (1 / C)) / n)  # Atualiza as estimativas de probabilidade
  }
  cat("Iterações = ", iter, "\n")
  cat("Máxima diferença = ", maxdiff, "\n")
  return(list(time = tau, surv = 1 - cumsum(p)))  # Retorna a função de sobrevivência
}

#############################################################

# Estimação da função de sobrevivência para 'var_cat' igual a 0 e 1 (sexos diferentes)
dat1 = df_surv[df_surv$var_cat == 0,]
dat1$right[is.na(dat1$right)] = Inf
tau = cria.tau(dat1)
p = S.ini(tau = tau)
A = cria.A(data = dat1, tau = tau)
tb1 = Turnbull(p, A, dat1)

dat1 = df_surv[df_surv$var_cat == 1,]
dat1$right[is.na(dat1$right)] = Inf
tau = cria.tau(dat1)
p = S.ini(tau = tau)
A = cria.A(data = dat1, tau = tau)
tb2 = Turnbull(p, A, dat1)

# Visualização das funções de sobrevivência para os dois grupos
par(mfrow = c(1, 1))
plot(tb1$time, tb1$surv, col = "red", type = "s", ylim = c(0, 1), xlim = c(0, 10),
     xlab = "Tempo em anos", ylab = "S(t)")
lines(tb2$time, tb2$surv, col = "blue", type = "s")
legend(0, 0.2, lty = 1, col = c("red", "blue"), c("Feminino", "Masculino"),
       bty = "n", cex = 0.9)


# Modelo de Cox para Censura Intervalar
# Análise Univariada

m = max(summary(df_surv$right))  # Máximo dos valores de censura
li = df_surv$left 
ui = ifelse(is.na(df_surv$right), m + 1000, df_surv$right)  # Ajuste para censura direita

# Ajuste do modelo de Cox para a variável categórica
fit1 = ic_sp(cbind(li, ui) ~ var_cat, model = 'ph', bs_samples = 100, data = df_surv)
summary(fit1)

# Ajuste do modelo de Cox para a variável contínua
fit2 = ic_sp(cbind(li, ui) ~ var_cont, model = 'ph', bs_samples = 100, data = df_surv)
summary(fit2)

# Seleção das variáveis com p-valor < 0,25 para a análise multivariada
fit3 = ic_sp(cbind(li, ui) ~ var_cat + var_cont, model = 'ph', bs_samples = 100, data = df_surv)
summary(fit3)

# Implementando a função 'step' para seleção de variáveis com base na log-verossimilhança

##########################

step_icph <- function(model, data) {
  vars <- attr(terms(model$formula), "term.labels")  # Variáveis do modelo
  current_model <- model
  best_llk <- current_model$llk  # Log-verossimilhança inicial
  
  # Processo iterativo de remoção de variáveis com base na log-verossimilhança
  for (var in vars) {
    formula_new <- as.formula(paste("cbind(li, ui) ~", paste(setdiff(vars, var), collapse = " + ")))
    model_new <- ic_sp(formula_new, model = "ph", bs_samples = 100, data = data)
    new_llk <- model_new$llk  # Log-verossimilhança do novo modelo
    
    # Diagnóstico e atualização do modelo
    if (new_llk > best_llk) {
      best_llk <- new_llk
      current_model <- model_new
      vars <- setdiff(vars, var)
    }
  }
  
  return(current_model)  # Retorna o modelo final
}

##########################

# Ajustando o modelo final após a seleção de variáveis
fit_final = step_icph(fit3, data = df_surv)
summary(fit_final)
