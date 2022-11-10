# ACT: Analysis of Contingency Tables
The purpose of this code is to offer a method to **analyze contingency tables and their residuals**.

This is an implementation of the method "Analysis of Contingency Tables" presented in García-Pérez, M.A., Núñez-Antón, V. & Alcalá-Quintana, R. *Analysis of residuals in contingency tables: Another nail in the coffin of conditional approaches to significance testing*. Behav Res 47, 147–161 (2015). https://doi.org/10.3758/s13428-014-0472-0. (NB.: I have no connection to the authors.)

The method avoids losing control of type I error rates, which typically happens with the very common method "performing an omnibus test first (e.g. a chi-square test of independence), and if the test is significant, analyzing the residuals". García-Pérez et al. explain in their paper why it should be avoided, and mention possible remedies.

There are more complete implementations of their method available in R and Matlab, developed by García-Pérez et al., and available for download here: https://static-content.springer.com/esm/art%3A10.3758%2Fs13428-014-0472-0/MediaObjects/13428_2014_472_MOESM1_ESM.zip

For the moment, the Python code available here only implements their method for testing independence. Moreover, it does not implement various "safety checks" the original authors use in their own code. So, until implemented here, you should take care of:
 - using a two-way contingency table;
 - using an alpha level $a$ where $0 < a \leqslant 0.05 $
 - using more than 1,000 replicates for the bootstrap (the authors recommend at least 30,000).
 - not using a table containing NaNs;
 - not using a table with empty rows or columns (i.e. summing to zero);

Numpy, scipy, and statsmodels are required to use the code, though the versions of these packages mentioned in the requirements.txt file are not set in stone. I have not tested the code against various versions of these packages -yet I guess it is quite likely to work with more recent versions, and probably with some older versions as well.

## Usage example
Using the ```ACT_I``` function available in the ```ACT.py``` file:

```
import numpy as np
observed = np.array([[ 1 , 7, 15, 12, 12, 14],
                  [ 1, 16, 22, 31, 32, 27],
                  [ 7, 14, 25, 28, 46, 44],
                  [13, 19, 34, 45, 63, 72]]
               )
alpha_level = 0.05 #should always be > 0 and <=0.05
residual_type = "ADJ" #‘ADJ’ for adjusted residuals, otherwise will use moment-corrected residuals
n_bootstrap = 50000 #number of replicates to generate during bootstrapping, García-Pérez et al. recommend a number of 30,000 at least
result = ACT_I(observed, alpha_level, residual_type, n_bootstrap)
print(result)
```

It returns a dictionary, formatted according to the convention used in the R implementation of this method. For a standard usage, the values of direct interest will generally be **Famwise_Significant** and **OmnibusHypothesis**. The documentation attached to the R and Matlab code gives additional useful information about the output (see [this link](https://static-content.springer.com/esm/art%3A10.3758%2Fs13428-014-0472-0/MediaObjects/13428_2014_472_MOESM1_ESM.zip)).

Output:
```
{'Problem': 'Omnibus test of independence and ADJ residual analysis',
 'InputTable': array([[ 1,  7, 15, 12, 12, 14],
        [ 1, 16, 22, 31, 32, 27],
        [ 7, 14, 25, 28, 46, 44],
        [13, 19, 34, 45, 63, 72]]),
 'NominalTestSize': 0.05,
 'NumReplicates': 50000,
 'ValidReplicates': 50000,
 'ExpectedFrequencies': array([[ 2.23666667,  5.69333333,  9.76      , 11.79333333, 15.555     ,
         15.96166667],
        [ 4.73      , 12.04      , 20.64      , 24.94      , 32.895     ,
         33.755     ],
        [ 6.01333333, 15.30666667, 26.24      , 31.70666667, 41.82      ,
         42.91333333],
        [ 9.02      , 22.96      , 39.36      , 47.56      , 62.73      ,
         64.37      ]]),
 'Residuals': array([[-0.88888376,  0.60679071,  1.93084752,  0.07069462, -1.10181252,
         -0.60289419],
        [-1.972219  ,  1.3527679 ,  0.36864605,  1.52490404, -0.20405401,
         -1.52719759],
        [ 0.48090235, -0.41146515, -0.30983645, -0.85979232,  0.87849342,
          0.22646775],
        [ 1.75778399, -1.12995049, -1.2135889 , -0.53807892,  0.0514188 ,
          1.44088905]]),
 'Cellwise_CriticalValue': 1.959963984540054,
 'Cellwise_Significant': array([[False, False, False, False, False, False],
        [ True, False, False, False, False, False],
        [False, False, False, False, False, False],
        [False, False, False, False, False, False]]),
 'Cellwise_ExactTestSize': 0.0485925,
 'Famwise_AlphaStar': 0.002155596168895355,
 'Famwise_CriticalValue': 3.067912650623893,
 'Famwise_Significant': array([[False, False, False, False, False, False],
        [False, False, False, False, False, False],
        [False, False, False, False, False, False],
        [False, False, False, False, False, False]]),
 'Famwise_ExactTestSize': 0.05002,
 'OmnibusHypothesis': 'Not rejected'}
 ```

## Development roadmap

- automatically checking the validity of inputs
- adding support for tests of homogeneity and tests of fit
- improving documentation
- making it available as a package on pipy
- generating better formatted reports
