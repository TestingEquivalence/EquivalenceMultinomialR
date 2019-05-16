This R project provides two tests to show equivalence of observed counting frequencies to a fully specified multinomial distribution.  
Thus, it is possible to show equivalence of the observed empirical data to the theoretical distribution. 

The package is based on the article:
Vladimir Ostrovski, “Testing equivalence of multinomial distributions”, Statistics and Probability Letters 124 (2017), 77–82.

Two examples are available in the script "examples.R“.
The program is rewritten in R and optimized for better performance and understanding. Particularly the optimization in the bootstrap
test is reworked considerable for a stable performance. Hence the results for real data sets deviate slightly from these, 
published in the article. However, the deviations are very small and do not matter for any application.

Let d denote the total variation distance and let q be the fully specified multinomial distribution.
Let p be the true underlying distribution of the observed data. 
The distribution p is equivalent to q if d(p,q)<e, where e>0 is a tolerance parameter.
The equivalence test problem is formally stated by 
H0={d(p,q) >=e } against H1={d(p,q)<e}.
The true distribution p is unknown. Instead we observe the counting frequencies p_n, where n is the number of observations. 

The goal is to reject the hypothesis of the non-equivalency H0 at a significance level alpha based on the counting frequencies p_n.

The project provides two tests for that purpose: 
the asymptotic test in "asymptotic_test.R";
the bootstrap test in "bootstrap_test.R". 

Both tests return the smallest tolerance parameter e, for which H0 can be rejected and hence equivalence can be shown.

The asymptotic test is based on the asymptotic distribution of the test statistic.
Therefore, the asymptotic test needs some sufficiently large number of observations. 
It should be used carefully because the test is approximate and may be anti-conservative at some points. 
To obtain a conservative test reducing of alpha (usually halving) or 
slight shrinkage of the tolerance parameter e may be appropriate. 

The bootstrap test is based on the re-sampling method called bootstrap. 
For the small sample sizes, the bootstrap test is more precise and reliable than the asymptotic test. 
However, it should be used carefully because the test is approximate and may be anti-conservative at some points. 
To obtain a conservative test reducing of alpha (usually halving) 
or slight shrinkage of the tolerance parameter e may be appropriate. 
The slight shrinkage of the tolerance parameter is preferable 
because even small shrinkage of the tolerance parameter is often sufficient.
