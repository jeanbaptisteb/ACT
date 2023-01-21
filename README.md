# ACT: Analysis of Contingency Tables
The purpose of this code is to offer a statistical method to **analyze contingency tables and their residuals**, with a bootstrap correction for multiple testing. 

**This version is still under development**, so any use in production should be tested against its original implementation in R (available for download [here](https://static-content.springer.com/esm/art%3A10.3758%2Fs13428-014-0472-0/MediaObjects/13428_2014_472_MOESM1_ESM.zip)).


##### Table of Contents  
- [Motivation](#introduction-and-motivation)  
- [Usage example](#usage-example)
- [Interpreting the output](#interpreting-the-output)
- [Requirements](#requirements)
- [Limitations](#limitations)
- [Development roadmap](#development-roadmap)

## Motivation
This is a Python implementation of the method "Analysis of Contingency Tables" (ACT) presented in García-Pérez, M.A., Núñez-Antón, V. & Alcalá-Quintana, R. *Analysis of residuals in contingency tables: Another nail in the coffin of conditional approaches to significance testing*. Behav Res 47, 147–161 (2015). https://doi.org/10.3758/s13428-014-0472-0.

The ACT method avoids losing control of type I error rates, which typically happens with the very common method "performing an omnibus test on the contingency table first, and if the test is significant, analyzing the residuals". García-Pérez et al. explain in their paper why it should be avoided, and suggest this alternative method as a possible remedy.

Note that there are other implementations of the ACT method, developed by García-Pérez et al. for R and Matlab, and available for download here: https://static-content.springer.com/esm/art%3A10.3758%2Fs13428-014-0472-0/MediaObjects/13428_2014_472_MOESM1_ESM.zip

## Usage example
The example below uses the ```ACT_I``` function available in the ```ACT.py``` file.
The function takes 5 parameters:
- `observed`: the contingency table to analyse (type: numpy array)
- `alpha`: significance level. Should always be > 0 and < 0.5. Default 0.05
- `Rtype`: type of residuals to use in the analysis. "ADJ" for adjusted residuals, "MC" for moment-corrected residuals. Default "ADJ".
- `nrep`: number of replicates to generate during bootstrapping. Default 30000, as recommended by García-Pérez et al.
- `seed`: seed state used to generate simulated tables. The purpose of this parameter is to make reproduction of the results easier. Default 0. *NB: this value is not present as a parameter in the R implementation of ACT.*


```
import numpy as np
#here is our contingency table
observed = np.array(
                 [[ 1 , 7, 15, 12, 12, 14],
                  [ 1, 16, 22, 31, 32, 27],
                  [ 7, 14, 25, 28, 46, 44],
                  [13, 19, 34, 45, 63, 72]]
               ) 
result = ACT_I(observed, alpha=0.05 , Rtype="ADJ" , nrep=50000, seed=0)
print(result)
```

**Output:**
```
{'Problem': 'Omnibus test of independence and ADJ residual analysis',
 'Seed': 0,
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
 'Cellwise_ExactTestSize': 0.048640833333333335,
 'Famwise_AlphaStar': 0.002218134534779443,
 'Famwise_CriticalValue': 3.059355919579737,
 'Famwise_Significant': array([[False, False, False, False, False, False],
        [False, False, False, False, False, False],
        [False, False, False, False, False, False],
        [False, False, False, False, False, False]]),
 'Famwise_ExactTestSize': 0.05,
 'OmnibusHypothesis': 'Not rejected'}
 ```
## Interpreting the output

The function returns a dictionary, formatted according to the convention used in the R implementation of this method. (Except for the 'Seed' value, which is not available in the R version of ACT).

**TL;DR**

For a standard use, the values of direct interest will generally be **Cellwise_Significant**, **OmnibusHypothesis**, and **Famwise_Significant**. The **Cellwise_Significant** value is an array, showing which cells have significant residuals (`True`)  or non-significant residuals (`False`). The  **OmnibusHypothesis** value shows if we can reject (or not) the null hypothesis (based on if any **Famwise_Significant** is `True`). The example above shows that none of the residuals is significant, so the omnibus test is non-significant (so we cannot reject the null hypothesis).

**Detailed explanation of the output**

- Problem:  Description of the analysis, mentioning the type of residuals defined in the Rtype parameter.
- Seed: seed state used to generate simulated tables. If exact reproduction of the results is needed, pass this value in the 'seed' parameter of the function when you call it again. *NB: this value is not present in the output of the R implementation of ACT.*
- InputTable: The array passed in the 'observed'  parameter.
- NominalTestSize: The alpha level defined in the 'alpha' parameter.
- NumReplicates: Number of replicates generated, defined by the 'nrep' parameter.
- ValidReplicates: Number of valid replicates. Invalid replicates may appear when the input table is too sparse, ending up with empty rows or columns. This kind of table is excluded from the analysis. A remedy may be to merge some rows or columns in the original table.
- ExpectedFrequencies: Expected cell frequencies under independence.
- Residuals: Cell residuals. Use statsmodels adjusted residuals when Rtype is 'ADJ' (statsmodels.stats.contingency_tables.Table.standardized_resids)
- Cellwise_CriticalValue: Critical value used in cellwise significance tests of residuals.
- Cellwise_Significant: a two-dimensional array of booleans, indicating which residuals are significant in cellwise tests (True == significant, False == non-significant)
- Cellwise_ExactTestSize: [TODO]
- Famwise_AlphaStar: Value of α* for familywise tests, defined through bootstraping
- Famwise_CriticalValue: Critical value to be used in the familywise significance test
- **Famwise_Significant**: Two-dimensional array of booleans, indicating which residuals are significant (or not) for the familywise test (True == significant, False == non-significant)
- Famwise_ExactTestSize: [TODO]
- **OmnibusHypothesis**: 'Rejected' or 'Not rejected'. Answers the question 'Is the omnibus hypothesis rejected?'. It is rejected if at least one item of Famwise_Significant is True.
    
The documentation attached to the R and Matlab code developed by García-Pérez et al. gives additional useful information about the output (see [this link](https://static-content.springer.com/esm/art%3A10.3758%2Fs13428-014-0472-0/MediaObjects/13428_2014_472_MOESM1_ESM.zip)).

 
## Requirements
Developed with Python 3.10.4. No tests have been performed against other versions of Python, though support for ulterior versions is on the roadmap.

Numpy, scipy, and statsmodels are required, though the versions of these packages mentioned in the [requirements.txt](https://github.com/jeanbaptisteb/ACT/blob/main/requirements.txt) file are not set in stone. The code hasn't been tested against various versions of these packages -yet it is quite likely to work with more recent versions, and probably with some older versions as well.

## Limitations
- This is a beta version, which still requires more extensive testing againt García-Pérez et al. implementation in R.
- For the moment, the Python code available here only implements García-Pérez et al.'s method for testing independence, not the other tests they mention (homogeneity...).

## Development roadmap

- testing more extensively against García-Pérez et al. implementation in R
- adding support for tests of homogeneity and tests of fit
- improving documentation
- testing support for Python > 3.10.4
- testing support for numpy > 1.22.0, scipy > 1.8.0, statsmodels > 0.13.2
- making it available as a package on pypi
- improving execution time
- generating better formatted reports? (heatmaps, pdf, etc.)
- make reproduction of results easier between R and Python
