import numpy as np
import pandas as pd
##from numpy import log, sum, abs, min, sign, arange
import random
##from random
import datetime
##from datetime.datetime import now
from scipy.stats import pearsonr

# GET FITNESS

def _get_multifitness(genes, nwt, nko, dwt, dko, obj):
    
    """\
    Parameters
    ----------
    genes : Genes of a chromosome
    nwt : Number of observations of Wild Type
    nko : Number of observations of Knock-out
    dwt : Data of observations of Wild Type
    dko : Data of observations of Knock-out
    obj : Array of arrays of fold-changes of the metabolites in all murine models

    Notes
    -----
    
     The geneSet will now be 13 characters (the same as mouse models) plus
     the 0. Each character represents a mouse model and let's take the cases in the
     Knock-Out dataset. Now for the wild type genes, I'm going to
     consider all those other than '0'. In the case of Knock-Outs, the '0'
     indicates that it is a sample that is not classified 8not optimized) by any other group.
     It could be outliers or simply not from the groups considered.
    
     I need to order what happens if in some solution they are not represented
     all models, what values does it return? How to operate with them? 

    """

    den = dwt.loc[[True if x != '0' else False for x in genes[:nwt]], :].mean()
    list_of_correlations = []
    
    for i, j in zip(list('ABCDEFGHIJKLM'), np.arange(13)):
        num = dko.loc[[True if x == i else False for x in genes[nwt:]], :].mean()
        myfc = np.log2(num / den)
        p, _ = pearsonr(myfc, obj[j])
        list_of_correlations.append(p)
        
    """\
    Resultados de traducción
    The function returns a single value that has to be a summary of the
    behavior of all samples.
    
    The most direct output would be the average correlation, waiting to reach a
    maximum which in turn would maximize all submodels.
    
    Another option is to return the minimum of all correlations, but this
    option can force the algorithm to search for a single model as a solution
    and the rest classify it as '0'.
    
    We could also think of some kind of penalty 
    """
    
    
    return None
    

def _get_fitness(genes, nwt, nko, dwt, dko, obj, l = 0.01):
    
    """\
    Parameters
    ----------
    genes : Genes of a chromosome
    nwt : Number of observations of Wild Type
    nko : Number of observations of Knock-out
    dwt : Data of observations of Wild Type
    dko : Data of observations of Knock-out
    obj : Array of fold-changes of metabolites in murine model
    l : float Usually 0.5 / length(fc1). If length(fc1) = 50 then l = 0.01
    
    Returns
    -------
    p : float
        The aptitude of a genes of a chromosome (solution) is the correlation between fc1 and fc2. 
        Initially the aptitude (fitness) of a chromosome is the Pearson Correlation between the 
        fold-change of Knock-out select samples and Wild Type selected samples with the fold-change of
        the murine model analyzed (a vector of ¿50? fold-changes).
    
    Notes
    -----
    This function evaluates the Pearson Correlation 
    between fc1 and fc2 with a penalty for different sign 
    in the values. This function uses pearsonr from scipy and
    sign from numpy to evaluate the aptitude.
    
    The penalty has been named l.
    
    Examples
    --------

    """
    
    num = dko.loc[[True if x == '1' else False for x in genes[nwt:]], :].mean()
    den = dwt.loc[[True if x == '1' else False for x in genes[:nwt]], :].mean()

    myfc = np.log2(num / den)
    
    p, _ = pearsonr(myfc, obj)
    penalty = sum(np.sign(myfc) * np.sign(obj) < 1) * l
    p = p - penalty
    
    return p  


def _get_fitness_robust(genes, nwt, nko, dwt, dko, obj, l = 0.01):
    
    """\
    Parameters
    ----------
    genes : Genes of a chromosome
    nwt : Number of observations of Wild Type
    nko : Number of observations of Knock-out
    dwt : Data of observations of Wild Type
    dko : Data of observations of Knock-out
    obj : Array of fold-changes of metabolites in murine model
    l : float Usually 0.5 / length(fc1). If length(fc1) = 50 then l = 0.01
    
    Returns
    -------
    p : float
        The aptitude of a genes of a chromosome (solution) is the correlation between fc1 and fc2. 
        Initially the aptitude (fitness) of a chromosome is the Pearson Correlation between the 
        fold-change of Knock-out select samples and Wild Type selected samples with the fold-change of
        the murine model analyzed (a vector of ¿50? fold-changes).
    
    Notes
    -----
    This function evaluates the Pearson Correlation 
    between fc1 and fc2 with a penalty for different sign 
    in the values. This function uses pearsonr from scipy and
    sign from numpy to evaluate the aptitude.
    
    The penalty has been named l.
    
    Examples
    --------

    """
    
    num = dko.loc[[True if x == '1' else False for x in genes[nwt:]], :].median()
    den = dwt.loc[[True if x == '1' else False for x in genes[:nwt]], :].median()

    myfc = np.log2(num / den)
    
    p, _ = pearsonr(myfc, obj)
    penalty = sum(np.sign(myfc) * np.sign(obj) < 1) * l
    p = p - penalty
    
    return p  



def _get_fitness_l(genes, nwt, nko, dwt, dko, obj):
    
    """\
    Parameters
    ----------
    genes : Genes of a chromosome
    nwt : Number of observations of Wild Type
    nko : Number of observations of Knock-out
    dwt : Data of observations of Wild Type
    dko : Data of observations of Knock-out
    obj : Array of fold-changes of metabolites in murine model
    
    Returns
    -------
    p : float
        The aptitude of a genes of a chromosome (solution) is the correlation between fc1 and fc2. 
        Initially the aptitude (fitness) of a chromosome is the Pearson Correlation between the 
        fold-change of Knock-out select samples and Wild Type selected samples with the fold-change of
        the murine model analyzed (a vector of ¿50? fold-changes).
    
    Notes
    -----
    This function evaluates the Pearson Correlation 
    between fc1 and fc2 with a penalty for different sign 
    in the values. This function uses pearsonr from scipy and
    sign from numpy to evaluate the aptitude.
    
    The penalty has been named l.
    
    Examples
    --------

    """
    
    l = 0.5/len(obj)
    
    num = dko.loc[[True if x == '1' else False for x in genes[nwt:]], :].mean()
    den = dwt.loc[[True if x == '1' else False for x in genes[:nwt]], :].mean()

    myfc = np.log2(num / den)
    
    p, _ = pearsonr(myfc, obj)
    penalty = sum(np.sign(myfc) * np.sign(obj) < 1) * l
    p = p - penalty
    
    return p  
 
