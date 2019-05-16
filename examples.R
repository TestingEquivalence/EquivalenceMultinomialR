source("asymptotic_test.R")
source("bootstrap_test.R")

# examples of apllication of the equivalence test to the real data sets
# we test equivalence of the observed counting frequencies to
# the fully specified multinomial distribution (categorial distribution) 
        

# 1. example: Mendel's inheritance law
# -------------------------------------------------

# counting observed by Mendel's (of pea plant seed)
p = c(315, 101, 108, 32)
# theoretic distribution, stated by Mendel's inheritance law
q = c(9 / 16, 3 / 16, 3 / 16, 1 / 16)
# number of observations
n=sum(p)
# smoothing parameter for the smooth total variation distance
b= 0.01 / sqrt(n)
# level
alpha=0.05

# test returns the minimum tolerance parameter epsilon,
# for which the equivalence between the theoretical and observed distribution can be shown
# at the significance level alpha
example1_asympt=asymptotic_test(p/n,q,n,b,alpha)
example1_bst=bootstrap_test(p/n,q,n,b,alpha)


# 2. example: Benford's law for Apple daily returns
# -------------------------------------------------

# counting of the first digit of the Apple daily returns
p = c(2001, 1359, 872, 625, 468, 417, 306, 251, 208)
# theoretic distribution, stated by Benford's law 
q=c(1:9)
q =log10(1 + 1 / q)
q=q/sum(q)

# number of observations
n=sum(p)
# smoothing parameter for the smooth total variation distance
b= 0.01 / sqrt(n)
# level
alpha=0.05

# test returns the minimum tolerance parameter epsilon,
# for which the equivalence between the theoretical and observed distribution can be shown
# at the significance level alpha
example2_asympt=asymptotic_test(p/n,q,n,b,alpha)
example2_bst=bootstrap_test(p/n,q,n,b,alpha)
