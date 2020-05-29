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
    
    El geneSet será ahora de 13 caracteres (los mismos que modelos de ratón) más
    el 0. Cada caracter representa un modelo de ratón y cogerenmos los casos en el 
    dataset de los Knock-Out. Ahora bien, para los genes de los wild type, voy a 
    considerar a todos los que no sean '0'. En el caso de los Knock-Out, el '0' 
    indica que es una muestra que no se clasifica 8no optimiza) ningún otro grupo. 
    Podría ser outliers o simplemente no ser de los grupos considerados.
    
    Me falta ordenar que ocurre si en alguna solución no están representados
    todos los modelos ¿qué valores devuelve? ¿cómo operar con ellos?

    """

    den = dwt.loc[[True if x != '0' else False for x in genes[:nwt]], :].mean()
    list_of_correlations = []
    
    for i, j in zip(list('ABCDEFGHIJKLM'), np.arange(13)):
        num = dko.loc[[True if x == i else False for x in genes[nwt:]], :].mean()
        myfc = np.log2(num / den)
        p, _ = pearsonr(myfc, obj[j])
        list_of_correlations.append(p)
        
    """\
    La función devuelve un único valor que tiene que ser un resumen del 
    comportamiento de todas las muestras.
    
    La salida más directa sería la correlación promedio, esperando en alcanzar un 
    máximo que a su vez maximizaría todas los submodelos.
    
    Otra opción es devolver el mínimo de todas las correlaciones, pero esta
    opción puede obligar al algoritmo a buscar como sólución un único modelo
    y el resto clasificarlo como '0'.
    
    Podríamos pensar también en algún tipo de penalización
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
 
