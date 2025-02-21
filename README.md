# Multivariate Normal Random Deviates

Multiplies uncorrelated normal deviates by a Cholesky factor of the covariance matrix to generate 
multivariate normal deviates with a specified covariance matrix. Sample output:
```
 #obs:     1000000

true means and standard deviations
Variable  1: Mean =  -5.00000000, StdDev =  10.00000000
Variable  2: Mean =   0.00000000, StdDev =  20.00000000
Variable  3: Mean =   5.00000000, StdDev =  30.00000000

empirical means and standard deviations
Variable  1: Mean =  -4.99696615, StdDev =   9.99908359
Variable  2: Mean =   0.00859363, StdDev =  19.99599058
Variable  3: Mean =   5.00536297, StdDev =  30.01204086

true correlation matrix:
    1.000000    0.500000    0.300000
    0.500000    1.000000    0.400000
    0.300000    0.400000    1.000000

empirical correlation matrix:
    1.000000    0.499548    0.299604
    0.499548    1.000000    0.400492
    0.299604    0.400492    1.000000

maxval(abs(xcorr - emp_corr)):     0.000492

empirical means and standard deviations
Variable  1: Mean =  -0.01287349, StdDev =   9.99949712
Variable  2: Mean =   0.00482072, StdDev =  20.00958772
Variable  3: Mean =  -0.02900624, StdDev =  30.02427445

empirical correlation matrix:
    1.000000    0.500036    0.299948
    0.500036    1.000000    0.400134
    0.299948    0.400134    1.000000

maxval(abs(xcorr - emp_corr)):     0.000134
```
