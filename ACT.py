# -*- coding: utf-8 -*-
"""Python implementation of the method presented in García-Pérez, M.A., Núñez-Antón, V. & Alcalá-Quintana, R. 
*Analysis of residuals in contingency tables: Another nail in the coffin of conditional approaches to significance testing*. 
Behav Res 47, 147–161 (2015). 
https://doi.org/10.3758/s13428-014-0472-0. """

import scipy
import statsmodels.api as sm
import numpy as np
from warnings import warn

def ACT_I(observed, alpha=0.05, Rtype="ADJ", nrep=30000):    
    """
    Testing a contingency table and its residuals for independence.
    
    Parameters
    ----------
    observed: numpy.array
        Observed contingency table. Must be two dimensional. Cannot contain NaN, and cannot contain empty rows or columns (i.e. summing to 0).
    alpha: float, default 0.05
        Significance level, between 0 and 0.5.
    Rtype: str, default 'ADJ'
        One of 'ADJ' (for adjusted residuals) or 'MC' (for moment-corrected residuals)
    nrep: int, default 30000
        Number of replicates to generate in the bootstrap procedure. Must be > 0. nrep >= 30000 is recommended.
    
    Returns
    ----------
    dict
        returns a dictionary with the following keys:
            - Problem: str
                Description of the analysis, mentionning the type of residuals defined in the Rtype parameter
            - InputTable: numpy.array
                The array passed in the 'observed'  parameter.
            - NominalTestSize: float
                The alpha level defined in the 'alpha' parameter.
            - NumReplicates: int
                Number of replicates generated, defined by the 'nrep' parameter.
            - ValidReplicates: int
                Number of valid replicates. Invalid replicates may appear when the input table is too sparse, ending up with empty rows or columns. This kind of table is excluded from the analysis. A remedy may be to merge some rows or columns in the original table.
            - ExpectedFrequencies: numpy.array
                Expected cell frequencies under independence.
            - Residuals: numpy.array
                Cell residuals. Use statsmodels adjusted residuals when Rtype is 'ADJ' (statsmodels.stats.contingency_tables.Table.standardized_resids)
            - Cellwise_CriticalValue: float
                Critical value used in cellwise significance tests of residuals.
            - Cellwise_Significant: numpy array
                Two-dimensional array of booleans, indicating which residuals are significant in cellwise tests (True == significant, False == non-significant)
            - Cellwise_ExactTestSize: float
            
            - Famwise_AlphaStar: float
                Value of α* for familywise tests, defined through bootstraping
            - Famwise_CriticalValue: float
                Critical value to be used in the familywise significance test
            - Famwise_Significant: numpy.array
                Two-dimensional array of booleans, indicating which residuals are significant for the familywise test (True == significant, False == non-significant)
            - Famwise_ExactTestSize: float
            
            - OmnibusHypothesis: str
                'Rejected' or 'Not rejected'. Answer the question 'Is the omnibus hypothesis rejected?'. It is rejected if at least one item of Famwise_Significant is True
    
    References
    ----------
    [1] García-Pérez, M.A., Núñez-Antón, V. & Alcalá-Quintana, R. 
    *Analysis of residuals in contingency tables: Another nail in the coffin of conditional approaches to significance testing*. 
    Behav Res 47, 147–161 (2015). 
    https://doi.org/10.3758/s13428-014-0472-0.
    """
    #TODO: implement various checks for the inputs    
    #check the input table for possible errors
    if isinstance(observed, np.ndarray) == False:
        raise  TypeError("'observed' must be a numpy array")
    if len(observed.shape) > 2 or np.size(observed) < 4:
        raise  ValueError("The contingency table is not a two-way table")
    if np.isnan(observed).any():
        raise ValueError("The contingency table cannot contain missing values")
    if 0 in observed.sum(axis=0) or 0 in observed.sum(axis=1):
        raise ValueError("The table cannot have empty rows or columns")
    #check the other parameters
    if Rtype.lower() not in ["adj", "mc"]:
        raise  ValueError("""Invalid type of residuals: 'Rtype' must be 'MC' for 'moment corrected residuals'
               or 'ADJ' for 'adjusted residuals'""")    
    if (isinstance(alpha, float) == False) or alpha <= 0 or alpha > 0.5:
        raise  ValueError('Invalid alpha (0 < alpha <= 0.5)')
    if isinstance(nrep, float):
        nrep = int(nrep)
    elif isinstance(nrep, int) == False or nrep < 0:
        raise  ValueError("nrep must be a positive integer")
    
    #warnings
    if nrep < 30000:
        warn("Consider using at least 30,000 replicates", stacklevel=2)
        
    if (np.product(observed.shape) * nrep) <1000:
        warn("consider increasing the number of replicates to at least "+str(round((1000/np.product(observed.shape))+1, 0)) + " to get enough valid replicates",
             stacklevel=2)
    nfil, ncolumn = observed.shape
    margin_col = observed.sum(axis=0)
    margin_row = observed.sum(axis=1).reshape(-1,1)
    n = observed.sum()
    table = sm.stats.Table(observed)
    expected = table.fittedvalues
    probabilities = table.independence_probabilities # expected/n
    if Rtype=='ADJ':
        variances = (1-margin_row/n)*(1-margin_col/n)
        residuals = table.standardized_resids
    else:
        variances = np.zeros(shape=observed.shape)
        variances.fill( (ncolumn-1)*(nfil-1)/(ncolumn*nfil) )
        residuals = (observed-expected)/np.sqrt(expected*variances) 
    rng = np.random.default_rng()
    simulations = rng.multinomial(n, 
                            probabilities.flatten(), 
                            size=nrep)
    
    simulations = [ sim.reshape(observed.shape) for sim in simulations]
    sim_expected = [np.zeros(shape=observed.shape) for i in range(nrep)]
    sim_variances = [np.zeros(shape=observed.shape) for i in range(nrep)]
    sim_residuals = [np.zeros(shape=observed.shape) for i in range(nrep)]
    for index, sim in enumerate(simulations):
        sim_table = sm.stats.Table(sim)
        sim_expected[index]  = sim_table.fittedvalues
        if Rtype=='ADJ':
            sim_margin_col = sim.sum(axis=0)
            sim_margin_row = sim.sum(axis=1).reshape(-1,1)
            sim_variances[index] = (1-sim_margin_row/n)*(1-sim_margin_col/n)
            sim_residuals[index] = sim_table.standardized_resids
        else:
            sim_variances[index].fill( (ncolumn-1)*(nfil-1)/(ncolumn*nfil) )
            sim_residuals[index] = (sim-sim_expected[index])/np.sqrt(sim_expected[index]*sim_variances[index]) 
    toKeep = np.isfinite([r.sum() for r in sim_residuals])
    valid = len(toKeep)
    if valid==0 :
        raise ValueError('Table is too sparse to produce valid replicates; consider merging rows or columns')
    elif valid <= nrep/2:
        warn('Table seems to be too sparse; consider merging rows or columns')    
    sim_residuals = [s for i, s in enumerate(sim_residuals) if toKeep[i] == True]
    total = np.size(sim_residuals)    
    zmin = scipy.stats.norm.ppf(1-alpha) 
    zmax = 10
    z_residuals = scipy.stats.norm.ppf(1-alpha/2)
    signif_residual =(abs(residuals)>z_residuals)
    #check if there are enough valid replicates
    if total > 1000:
        for i in range(25):
            z_omnibus = (zmax+zmin)/2
            type1 = np.mean([True in (abs(sim) > z_omnibus) for sim in sim_residuals])
            if type1>alpha:
                zmin = z_omnibus
            else:
                zmax = z_omnibus
        
        type1_cell = np.mean(np.absolute(sim_residuals)>z_residuals)        
        alpha_star = 2*scipy.stats.norm.cdf(-z_omnibus)    
        signif_omnibus = (abs(residuals)>z_omnibus)
        
        if signif_omnibus.any():
            omnibus_test = 'Rejected'
        else:
            omnibus_test = "Not rejected"
        report = {"Problem": " ".join(['Omnibus test of independence and'
                                           ,Rtype,'residual analysis']),
                        "InputTable":  observed,
                        "NominalTestSize":  alpha,
                        "NumReplicates":  nrep,
                        "ValidReplicates":  valid,
                        "ExpectedFrequencies":  expected,
                        "Residuals":  residuals,
                        "Cellwise_CriticalValue":  z_residuals,
                        "Cellwise_Significant":  signif_residual,
                        "Cellwise_ExactTestSize":  type1_cell,                
                        "Famwise_AlphaStar":  alpha_star,
                        "Famwise_CriticalValue":  z_omnibus,
                        "Famwise_Significant":  signif_omnibus,
                        "Famwise_ExactTestSize":  type1,
                        "OmnibusHypothesis":  omnibus_test,                        
                        }
    else:
        report = {"Problem": " ".join(['Omnibus test of independence and'
                                           ,Rtype,'residual analysis']),
                        "InputTable":  observed,
                        "NominalTestSize":  alpha,
                        "NumReplicates":  nrep,
                        "ValidReplicates":  valid,
                        "ExpectedFrequencies":  expected,
                        "Residuals":  residuals,
                        "Cellwise_CriticalValue":  z_residuals,
                        "Cellwise_Significant":  signif_residual,
                        "Cellwise_ExactTestSize":  "Not computed; insufficent valid replicates",                
                        "Famwise_AlphaStar":  "Not computed; insufficent valid replicates",
                        "Famwise_CriticalValue":  "Not computed; insufficent valid replicates",
                        "Famwise_Significant":  "Not computed; insufficent valid replicates",
                        "Famwise_ExactTestSize":  "Not computed; insufficent valid replicates",
                        "OmnibusHypothesis":  "Not computed; insufficent valid replicates",                        
                        }
    return report
