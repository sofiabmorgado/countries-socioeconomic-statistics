#IMPORTAR BIBLIIOTECAS
library(readxl)
library(xlsx)
library(cluster)
library(factoextra)
library(gridExtra)

#IMPORTAR OS DADOS
dados0 <- read_excel("C:/Users/sofia/OneDrive/Ambiente de Trabalho/FCT/estatistica_multivariada/projeto/Country-data.xlsx", col_names = TRUE); dados 
dados0[,-1]


##################################################
####ACP####
##################################################

#Calcular a matriz de covariancias
S <- cov(dados0[,-1]); S

#Avaliar variancia das variaveis para decidir se standardizamos os dados
diag(S) #apenas queremos olhar para a variancia - como as variancias são muito diferentes, queremos standardizar os dados

#Standartizar as variaveis
dados <- scale(dados0[,-1], center=TRUE, scale= TRUE); dados


#Calcular matriz de covariancias dos dados standardizados
S <- cov(dados); S

#Numero de variaveis
p <- ncol(dados); p
sum(eigen(S)$value) #é igual a sum(diag(S)) #A soma dos vetores proprios corresponde ao numero de variaveis pq a variancia entre cada variavel é 1

#Numero de exemplos
n <- nrow(dados); n

#Calcular os valores e vetores próprios
val_prop <- eigen(S)$value; val_prop
vet_prop <- eigen(S)$vector; vet_prop #em cada coluna temos 1 vetor proprio, com os valores correspondentes a cada variavel 

#Exemplo de uma combinacao linear que define o componente principal
#y1 <- eigen(S)$vectors[1,1] * x1 + eigen(S)$vectors[2,1] * x2 + eigen(S)$vectors[3,1] * x3 + eigen(S)$vectors[4,1] * x4 + eigen(S)$vectors[5,1] * x5 + eigen(S)$vectors[6,1] * x6 + eigen(S)$vectors[7,1] * x7 + eigen(S)$vectors[8,1] * x8 + eigen(S)$vectors[9,1] * x9



###########
#SELECAO DAS COMPONENTES PRINCIPAIS
###########


#ACP INFO
#Desvio padrao
dp <- sqrt(val_prop); dp #é a sqrt da diagonal de S ;)

#Proporcao de variancia retida 
proporcao <- eigen(S)$value/sum(eigen(S)$value); proporcao #que cada um retem (multiplicar por 100 para percentagem)

#Percentagem de variancia retida 
percentagem <- eigen(S)$value/sum(eigen(S)$value) * 100; percentagem

#Soma das proporcoes - proporcao retida comulativa
sum_proporcao <- cumsum(proporcao); sum_proporcao

#O numero da CP
k <- 1:9; k#k é o numero de variveis

#Matriz com os dados
PCA_info <- matrix(c(k, val_prop, dp, proporcao, sum_proporcao), 9 , 5); PCA_info


#Critério de Kaiser
#Em vari´aveis padronizadas, uma componente com valor pr´oprio inferior a 1 indica que n~ao retem
#a informa¸c~ao (variabilidade) equivalente a uma das vari´aveis originais (com vari^ancia 1, por
#estarem padronizadas). Assim, de acordo com este crit´erio, em vari´aveis padronizadas, s~ao
#retidas as CP com valor proprio superior a 1.
#Assim, apenas retemos as primeiras 3 CP

#Scree plot: como podemos ver nao ha um ponto entre as 3 primeiras CP onde a inclinação mude acentuadamente
#Gráficos
plot(k, val_prop, type="b", main = "Scree Plot", xlab = "Componentes principais", ylab = "Valor próprio") 
plot(k, sum_proporcao, type="b", main = "Variabilidade retida cumulativa", xlab = "Componentes principais", ylab = "Proporção de variabilidade retida cumulativa")


#Portanto, selecionamos as primeiras 3 componentes principais
val_prop_final <- eigen(S)$value[1:3]; val_prop_final
vet_prop_final <- eigen(S)$vector[,1:3]; vet_prop_final
percentagem_final <- sum(val_prop_final)/sum(eigen(S)$value) * 100; percentagem_final



###########
#ESTUDO DA ADQUABILIDADE DA ACP
###########


#PROPORCAO DE VARIAVEIS CORRELACIONADAS
num <- (sum(abs(S)>0.5)-p)/2
#numero de entradas na matriz correspondentes a correlacoes (tira se a diagonal principal -p, e divide se por 2 pois e simetrica)
den <- p * (p-1)/2 
num/den #percentagem de correlacoes bivariadas= 27.78%


#CORRELACOES PARCIAIS
pcor <- corpcor::cor2pcor(S); round(pcor, 2)
  #quantas sao correlacionadas pela correlacoes parciais
num_parcial <- (sum(abs(pcor)<0.5) )/2 #aqui n se tira as diag principal pq ela e sempre sup a 1, mas divide se por 2
  #numero de entradas na matriz correspondentes a correlacoes (tira se a diagonal principal)
den_parcial <- p * (p-1)/2 
num_parcial/den_parcial #percentagem de correlacoes bivariadas parciais= 88.89%




#TESTE DE ESFERICIDADE
#Valores elevados do U_asterisco e baixos valores de p value, pois aqui as variaveis tem correlacoes = 0
X <- dados
p <- ncol(X); p
n <- nrow(X); n

U <- p^p * det(S) / sum(diag(S))^p ; U
U_asterisco <- -(n - 1 - (2 * p^2 + p +2) / (6 * p) )* log(U) ; U_asterisco

#p_value
p_value <- 1 - pchisq(U_asterisco, p*(p+1)/2 -1); p_value #portanto, aqui utilizar o ACP faz sentido


#Teste de esfericidade BUILT IN - Mauchly test - Podemos confirmar os valores obtidos no teste de esfericidade
#mauchly.test(lm(dados~1))

#O resutado sugere que a ACP é adequada, uma vez que não existem evidências a favor de uma estrutura de covariância esférica


###########
#CRIAR NOVA MATRIZ COM CP
###########

  #Aqui usamos os dados normalizados

cp1 <- list()
cp2 <- list()
cp3 <- list()

for (i in 1:n) {
  v <- dados[i,]
  vet_1 <- vet_prop_final[,1]
  cp1[i] <- vet_1[1] * v[1] + vet_1[2] * v[2] + vet_1[3] * v[3] + vet_1[4] * v[4] + vet_1[5] * v[5] + vet_1[6] * v[6] + vet_1[7] * v[7] + vet_1[8] * v[8] + vet_1[9] * v[9]
  
  vet_2 <- vet_prop_final[,2]
  cp2[i] <- vet_2[1] * v[1] + vet_2[2] * v[2] + vet_2[3] * v[3] + vet_2[4] * v[4] + vet_2[5] * v[5] + vet_2[6] * v[6] + vet_2[7] * v[7] + vet_2[8] * v[8] + vet_2[9] * v[9]
  
  vet_3 <- vet_prop_final[,3]
  cp3[i] <- vet_3[1] * v[1] + vet_3[2] * v[2] + vet_3[3] * v[3] + vet_3[4] * v[4] + vet_3[5] * v[5] + vet_3[6] * v[6] + vet_3[7] * v[7] + vet_3[8] * v[8] + vet_3[9] * v[9]
  
}

cp <- matrix(c(cp1, cp2, cp3), n, 3); cp




##################################################
####CLUSTERS####
##################################################


###########
#HIERARCHICHAL CLUSTERING
###########

#usar cp como dataframe

#Criar matriz de distancias 
d <- dist(cp, method = "euclidean") #tiramos coluna com nome das cidades
d


#Metodo de ligacao simples -cria 2 clusteres 1 com muito poucos paises, pouca diferenca tambem da distancia
fit2 <- hclust(d, method = "single")
plot(fit2, hang = 0.1, labels = dados0$country)
rect.hclust(fit2, k = 2, border = "red")

#grafico melhor
fviz_dend(fit2, rect = TRUE, cex = 0.1, k_colors = c("#00AFBB","#2E9FDF"))


#Metodo de ligacao completa -cria 2 clusteres 1 com muito poucos paises, pouca diferenca tambem da distancia
fit3 <- hclust(d, method = "complete")
plot(fit3, hang = 0.1, labels = dados0$country, cex=0.5)
rect.hclust(fit3, k = 2, border = "red")

#grafico melhor
fviz_dend(fit3, rect = TRUE, cex = 0.1, k_colors = c("#00AFBB","#2E9FDF"))


#Metodo de ligacao Ward #duvidas usar ?hclust
fit1 <- hclust(d, method = "ward.D")
plot(fit1, hang = 0.1, labels = dados0$country)
rect.hclust(fit1, k = 2, border = "red")

#grafico melhor
fviz_dend(fit1, rect = TRUE, cex = 0.1, k_colors = c("#00AFBB","#2E9FDF"))



###########
#K-MEANS
###########

Z <- as.data.frame(cp); Z #dados[,-1]
#Z <- cp; Z

#Algoritmo K-MEANS
set.seed(1)
c <- 3 #numero de centroides
fitkm <- kmeans(Z, centers = c, nstart = 3); fitkm #nstart é quantas vezes ele corre o algoritmo!


#Dados do k-means
fitkm #cada centroide é um vetor de 9 variaveis
fitkm$centers
#Coordenadas dos centroides -medias 
centroides <- aggregate(Z, list(fitkm$cluster), mean); centroides
#avaliar a dispersao dentro de cada cluster
fitkm$withinss #o 3 e o ultimo sao os mais heterogeneos 
#avaliar dispersao entre clusters 
fitkm$betweenss 




#DETERMINACAO DO MELHOR NUMERO DE CLUSTERS
#Elbow method
fviz_nbclust(Z, kmeans, nstart = 3, method = "wss") #dificil de vizualizar o cotovelo, talvez se veja no k = 2 ou 5

#Silhouette method
fviz_nbclust(Z, kmeans, nstart = 3, method = "silhouette") #avaliar o silhouete medio, queremos o maior possivel, e tao proximo de 1 quanto possivel pois assim os nossos objetos estao mais bem agregados

#Silhouette score para cada objeto 
k3 <- kmeans(Z, centers = 3, nstart = 3); k3
tab <- silhouette(k3$cluster, dist(Z)) #assim temos o valor de silhouete para cada objeto; quanto maior, mais bem colocado esta o objeto no cluster
max(tab[,3]) #da nos o valor maximo de tab, ou seja, o maior silhouette score
plot(tab) #Representar graficamente o indice silhouette para cada uma das observacoes #informa sobre cada clister, qual o numero de objetos e o seu silhouette medio (1: 4 | 0.21), podemos observar quais os exemplos mais bem colocados 
#para maximizar o silhouette score, selecionamos numero de clusters = 2


#Criar dataframe com clusters ja 
dd <- cbind(dados0[,1], cp, cluster = fitkm$cluster); dd
dd$cluster
dd[, 5]




#GRAFICOS
library(scatterplot3d)
scatterplot3d(dd[,2:4], pch = 1, lwd = 0.5, color=dd[,5], main = "Clusters com K-means", 
              xlab = "CP1", ylab = "CP2", zlab = "CP3", grid = FALSE, box = FALSE)

text(dd, labels = dd[,1],cex= 0.7, col = "steelblue")

#3D PLOT 1
library(plotly)
library(dplyr)

p <- plot_ly(dd, x=dd[,2], y=dd[,3], z=dd[,4], color = dd[,5], colors = c('#BF382A', '#0C4B8E', "goldenrod")) 
p <- p %>% layout(scene = list(xaxis = list(title = 'CP1'),
                                   yaxis = list(title = 'CP2'),
                                   zaxis = list(title = 'CP3')))
print(p)


"""
dd$cluster[which(dd$cluster == 1)] <- 'Cluster 1'
dd$cluster[which(dd$cluster == 2)] <- 'Cluster 2'
dd$cluster[which(dd$cluster == 3)] <- 'Cluster 3'
"""

##gg3D
library(ggplot2)
library(cli)
library(devtools)
library(usethis)
library(gg3D)

#devtools::install_github("AckerDWM/gg3D", force = TRUE)

ggplot(dd, aes(x=dd[,2], y=dd[,3], z=dd[,5], color=dd[,5])) + 
  theme_void() +
  axes_3D() +
  stat_3D()


#PLOT THE CLUSTERS in 2D
fviz_cluster(fitkm, axes = c(1, 2), data = Z)
+ ggtitle("k=2")#para representar graficamente ele usa os 2 primeiros PCA



#2D plot
fviz_cluster(k3, data = Z)
#+ ggtitle("k=5")

#2D ggplot
dd <- as.data.frame(dd); dd
dd2 <- as.numeric(dd[,2]); dd2
dd3 <- as.numeric(dd[,3]); dd3
dd5 <- as.numeric(dd[,5]); dd5

library(ggplot2)
ggplot(dd, aes(x=dd[,2], y=dd[,3], color=dd[,5])) +
  geom_point()


#Get sample from the clusters
set.seed(123)
g = sample(1:167, 25, replace=FALSE); g

print("Sample from data")
for (i in g){
  print(dd[i,])
}
