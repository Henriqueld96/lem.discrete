\name{lem.exact.mcnemar}
\alias{lem.exact.mcnemar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{lem.exact.mcnemar
%%  ~~function to do ... ~~
}
\description{
Teste de McNemar exato (b+c < 20) bilateral (H1:p1 diff p2)}
\usage{
lem.exact.mcnemar(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{Tabela 2x2 com pares concordantes e discordantes
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Retorna as duas contas de combinatoria, T1 e T2, e o p-valor do teste (2*min{T1,T2})
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
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
