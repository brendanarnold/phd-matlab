function P=PofChisqd(chisqd,n)

P=2^-(n/2)/gamma(n/2)*sqrt(chisqd)^(n-2)*exp(-chisqd/2);