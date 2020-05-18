\name{lem.f.retas}
\alias{lem.f.retas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{lem.f.retas
%%  ~~function to do ... ~~
}
\description{
Calcula a estatística F na forma matricial
}
\usage{
lem.f.retas(C, B, X, Y, df1, df2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{C}{
Matriz de contraste
}
  \item{B}{
Vetor coluna estimadores B
}
  \item{X}{
Matriz X do modelo
}
  \item{Y}{
Vetor coluna com variavel resposta
}
  \item{df1}{
Numero de linhas em C
}
  \item{df2}{
Numero de amostras - número de estimadores B
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (C, B, X, Y, df1, df2)
{
    EstatF <- ((t(C \%*\% B) \%*\% solve(C \%*\% (solve(t(X) \%*\% X)) \%*\%
        t(C)) \%*\% (C \%*\% B))/df1)/((t(Y) \%*\% Y - t(B) \%*\% t(X) \%*\%
        Y)/df2)
    pvalor = 1 - pf(EstatF, df1, df2)
    re <- cbind(EstatF, pvalor)
    colnames(re) <- c("Estat Teste", "Pvalor")
    return(re)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
