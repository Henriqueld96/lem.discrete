#b,c pares discordantesd
lem.exact.mcnemar <- function(df){
  b = df[1,2]
  c = df[2,1]
  m1 = (1/2)^(b+c)
  i = 1
  while (i<=b) {
    m1 = m1 + choose(b+c,i)*(1/2)^(b+c)
    i = i+1
  }
  m1



  m2 = 0
  i = b
  while (i<=b+c) {
    m2 = m2 + choose(b+c,i)*(1/2)^(b+c)
    i = i+1
  }
  m2

  pvalor = 2*min(c(m1,m2))

  re <- cbind(m1,m2,pvalor)
  colnames(re) <- c("T1", "T2", "Pvalor")
  return(re)
}

#n1,n2 sao marginais, x1 valor superior esquerdo da tabela, y marginal primeira coluna
lem.exact.fisher <- function(df){
  n1 = sum(df[1,])
  n2 = sum(df[2,])
  x1 = df[1,1]
  y = sum(df[,1])
  teste <- c(0:y)
  s=c()
  obs = choose(n1, x1)*choose(n2,  y-x1)/choose(n1+n2, y)
  for (i in teste) {
    temp = choose(n1, i)*choose(n2,  y-i)/choose(n1+n2, y)
    s=append(s, temp, after = length(s))
  }
  pvalor = sum(s[which(s<=obs)])

  re <- cbind(pvalor)
  colnames(re) <- c("Pvalor")
  return(pvalor)
}


#Retorna a tabela com os valores esperados de uma dada tabela
#df = tabela
lem.expected <- function(df){
  m = data.frame()
  n = sum(df)

  for (i in 1:dim(df)[1]){
    for (j in 1:dim(df)[2]){
      m[i,j] = sum(df[i,])*sum(df[,j])/n
    }
  }
  return(m)
}


#criar os vetores "n" e "e" usando as caselas (por linha) das tabelas observada e esperada, respectivamente
lem.qpearson <- function(df){
  esperado = lem.expected(df)
  vetordf = as.vector(t(df))
  vetoresp = as.vector(t(esperado))

  s = c()
  for (i in 1:length(vetordf)){
    s = append(s, values = ((vetordf[i]-vetoresp[i])^2)/vetoresp[i], after = length(s))
  }
  estat_teste = sum(s)
  gl= (dim(df)[1]-1)*(dim(df)[2]-1)
  pvalor <- pchisq(estat_teste, df=gl, lower.tail = F)
  re <- cbind(estat_teste, pvalor)
  colnames(re) <- c("Estat Teste", "Pvalor")
  return(re)
}

#Calcula a matriz de variancia covariancia para o teste MH (dep: lem.expected)
#df = tabela
lem.cov.mh <- function(df){
  esperado <- lem.expected(df)
  m = data.frame()
  n = sum(df)
  colunas = c()
  for (i in 1:dim(df)[1]){
    for (j in 1:dim(df)[2]){
      colunas = append(colunas, paste("X",i,j, sep=""))
      t <- c()
      for (i2 in 1:dim(df)[1]){
        for (j2 in 1:dim(df)[2]){

          if (i==i2){
            deltai=1
          }else{
            deltai=0
          }
          if (j==j2){
            deltaj=1
          }else{
            deltaj=0
          }

          conta = (esperado[i,j]*(n*deltai-sum(df[i2,]))*(n*deltaj-sum(df[,j2])))/(n*(n-1))

          t = append(t, conta, after=length(t))
        }
      }
      m = rbind(m,t)
    }
  }
  colnames(m) = colunas
  rownames(m) = colunas
  return(t(m))
}

#Calcula a estat do teste de MH e o seu p-valor (dep: lem.expected, lem.cov.mh)
#df = tabela, a = matriz contraste
lem.qmh <- function(df, a){
  esperado <- lem.expected(df)

  a <- as.matrix(a)
  v <- as.matrix(lem.cov.mh(df))
  nm <- as.matrix(as.vector(t(df))-as.vector(t(esperado)))

  estat_teste <- t(nm)%*%t(a)%*%solve(a%*%v%*%t(a))%*%a%*%nm

  gl= (dim(df)[1]-1)*(dim(df)[2]-1)

  pvalor <- pchisq(estat_teste[1], df=gl, lower.tail = F)


  re <- cbind(estat_teste[1], pvalor)
  colnames(re) <- c("Estat Teste", "Pvalor")
  return(re)
}

#Calcula a estat do teste de MH com variavel ordinal e o seu p-valor (dep: lem.expected, lem.cov.mh)
#df = tabela, a = matriz contraste
lem.qmh.ord <- function(df, a){
  esperado <- lem.expected(df)

  a <- as.matrix(a)
  v <- as.matrix(lem.cov.mh(df))
  nm <- as.matrix(as.vector(t(df))-as.vector(t(esperado)))

  estat_teste <- t(nm)%*%t(a)%*%solve(a%*%v%*%t(a))%*%a%*%nm

  gl= (dim(df)[1]-1)

  pvalor <- pchisq(estat_teste[1], df=gl, lower.tail = F)


  re <- cbind(estat_teste[1], pvalor)
  colnames(re) <- c("Estat Teste", "Pvalor")
  return(re)
}

#calcula teste expandido de mh
#df é as tabelas separadas por extratos, uma embaixo da outra
#ntabs = numero de extratos, a = matriz contraste
#exemplo tabela ex3a-lista2
#df <- rbind(c(6,6,3),
#             c(2,7,6),
#             c(4,7,6),
#             c(2,4,11),
#             c(11,19,6),
#             c(6,12,17))
lem.qemh <- function(df, ntabs, a){

  colunas = dim(df)[1]/ntabs
  inicio = seq(from=1, to=dim(df)[1], by=colunas)
  fim = inicio+colunas-1

  pt <- 0
  st <- 0
  tt <- 0
  for (t in 1:length(inicio)){
    #calculando sum (n-m)'a' e a(n-m)
    esperado <- lem.expected(df[inicio[t]:fim[t],])
    nm <- as.matrix(as.vector(t(df[inicio[t]:fim[t],]))-as.vector(t(esperado)))
    pt =  pt + t(nm)%*%t(a)
    tt =  tt + a%*%nm

    #calculando sum (aVa')
    v <- as.matrix(lem.cov.mh(df[inicio[t]:fim[t],]))
    st = st + a%*%v%*%t(a)

  }

  estat_teste <- pt%*%solve(st)%*%tt

  gl= (dim(df[inicio[1]:fim[1],])[1]-1)*(dim(df[inicio[1]:fim[1],])[2]-1)

  pvalor <- pchisq(estat_teste[1], df=gl, lower.tail = F)


  re <- cbind(estat_teste[1], pvalor)
  colnames(re) <- c("Estat Teste", "Pvalor")
  return(re)
}

