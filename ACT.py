# -*- coding: utf-8 -*-
"""Python implementation of the method presented in García-Pérez, M.A., Núñez-Antón, V. & Alcalá-Quintana, R. 
*Analysis of residuals in contingency tables: Another nail in the coffin of conditional approaches to significance testing*. 
Behav Res 47, 147–161 (2015). 
https://doi.org/10.3758/s13428-014-0472-0. """

import scipy
import statsmodels.api as sm
import numpy as np

def ACT_I(observed, alpha, Rtype, nrep):
    #TODO: implement various checks for the inputs    
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
        return 'Table is too sparse to produce valid replicates; consider merging rows or columns'
    elif valid <= nrep/2:
        print('Table seems to be too sparse; consider merging rows or columns')
    
    sim_residuals = [s for i, s in enumerate(sim_residuals) if toKeep[i] == True]
    total = len(sim_residuals) 
    if total < 1000: #TODO: to change later, to allow returning certain value even if total < 1000
        return "You should specify a number of replicates > 1000"
    zmin = scipy.stats.norm.ppf(1-alpha) 
    zmax = 10
    for i in range(25):
        z_omnibus = (zmax+zmin)/2
        type1 = np.mean([True in (abs(sim) > z_omnibus) for sim in sim_residuals])
        if type1>alpha:
            zmin = z_omnibus
        else:
            zmax = z_omnibus
    z_residuals = scipy.stats.norm.ppf(1-alpha/2)
    type1_cell = np.mean(np.absolute(sim_residuals)>z_residuals)        
    alpha_star = 2*scipy.stats.norm.cdf(-z_omnibus)

    signif_omnibus = (abs(residuals)>z_omnibus)
    signif_residual =(abs(residuals)>z_residuals)
    if signif_omnibus.any():
        omnibus_test = 'Rejected'
    else:
        omnibus_test = "Not rejected"
    #TODO: If < 1000, some value may be returned, but not all of them should not be returned
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
                    "OmnibusHypothesis":  omnibus_test
                    }
    return report
