import numpy as np
import pandas as pd

from ..tools.tools import chromosome

def _crossover(chromosome_a, chromosome_b, get_fitness, nwt, nko, dwt, dko, obj, k = 1):
    
    """\
    Parameters
    ----------
    chromosome_a: Chromosome
    chromosome_b: Chromosome
    k: k-point. If k = 1 the method will be single-point crossover, 
                If k = 2 the method will be two-point crossover. 
                If k = N (with N the length of chromosome) the method will be the uniform crossover. Its value is 1 by default.
    
    Returns
    -------
    Two new chromosomes
    
    Notes
    -----
    [1]
    
    References
    ----------
    [1] https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)
    
    Examples
    --------
    >>> A = _get_father(nwt = dwt.shape[0], nko = dko.shape[0], 
                        dwt = dwt, dko = dko,
                        geneSet = '01', 
                        get_fitness = _get_fitness, obj = np.array(FC.loc[:, murine_models_in[4]]))

    >>> B = _get_father(nwt = dwt.shape[0], nko = dko.shape[0], 
                        dwt = dwt, dko = dko,
                        geneSet = '01', 
                        get_fitness = _get_fitness, obj = np.array(FC.loc[:, murine_models_in[4]]))

    >>> A.print()
    011000111111010 ... 001000010001110 || 001010011000110 ... 011001001111111 	  43	 336	0.03363686825576087

    >>> B.print()
    100101111101111 ... 010001000100101 || 010101001001010 ... 111010000101010 	  42	 345	0.1410663212679354

    >>> A1, B1 = _crossover(A, B, get_fitness = _get_fitness,  
                        nwt = dwt.shape[0], nko = dko.shape[0], dwt = dwt, dko = dko, 
                        obj = np.array(FC.loc[:, murine_models_in[4]]), k = 1)

    >>> A1.print()
    011000111111010 ... 001000010001110 || 001010011000110 ... 111010000101010 	  43	 327	0.007288743389698732

    >>> B1.print()
    100101111101111 ... 010001000100101 || 010101001001010 ... 011001001111111 	  42	 354	0.15370276754212897

    >>> A1, B1 = _crossover(A, B, get_fitness = _get_fitness,  
                        nwt = dwt.shape[0], nko = dko.shape[0], dwt = dwt, dko = dko, 
                        obj = np.array(FC.loc[:, murine_models_in[4]]), k = 772)

    >>> A1.print()
    110000111101111 ... 000000010001111 || 011111001001110 ... 011000001111111 	  43	 355	0.017131591198246282

    >>> B1.print()
    001101111111010 ... 011001000100100 || 000000011000010 ... 111011000101010 	  42	 326	0.13633106854522628 
    """
    
    a = chromosome_a.Genes
    b = chromosome_b.Genes
    
    if k < 1:
        k = 1
        
    if k > len(a):
        k = len(a)
    
    cross = np.sort(np.random.choice(a = np.arange(len(a)), size = k, replace = False))

    a1 = list(a).copy()
    b1 = list(b).copy()
    for idx in cross:
        a1[idx:], b1[idx:] = b1[idx:].copy(), a1[idx:].copy()

    fitness_a1 = get_fitness(''.join(a1), len(chromosome_a.WT_Genes), len(chromosome_a.KO_Genes), dwt, dko, obj = obj)
    fitness_b1 = get_fitness(''.join(b1), len(chromosome_b.WT_Genes), len(chromosome_b.KO_Genes), dwt, dko, obj = obj)
    
    return chromosome(''.join(a1), fitness_a1, nwt, nko), chromosome(''.join(b1), fitness_b1, nwt, nko)
   
    
    
def _arithmetic_crossover(chromosome_a, chromosome_b, get_fitness, nwt, nko, dwt, dko, obj, arithmetic = 'and'):
    
    """\
    Notes
    -----
    For this specific genotype (GeneSet = '01') the arithmetic crossover techniques are available.
    These techniques are:
      [1] AND
      [2] OR
      [3] XOR
    Arithmetic crossover - some arithmetic operation is performed to make a new offspring [1]
    
    References
    ----------
    [1] https://www.obitko.com/tutorials/genetic-algorithms/crossover-mutation.php
    """

    a = chromosome_a.Genes
    b = chromosome_b.Genes
    
    c_and = []
    c_or = []
    c_xor = []
    
    for i in np.arange(len(a)):
        
        ai = int(list(a)[i])
        bi = int(list(b)[i])
        
        c_and.append(str(ai and bi))
        c_or.append(str(ai or bi))
        c_xor.append(str(xor(ai, bi)))
        
    c_and = ''.join(c_and)
    c_or = ''.join(c_or)
    c_xor = ''.join(c_xor)
    
    fitness_c_and = get_fitness(c_and, len(chromosome_a.WT_Genes), len(chromosome_a.KO_Genes), dwt, dko, obj = obj)
    fitness_c_or = get_fitness(c_or, len(chromosome_a.WT_Genes), len(chromosome_a.KO_Genes), dwt, dko, obj = obj)
    fitness_c_xor = get_fitness(c_xor, len(chromosome_a.WT_Genes), len(chromosome_a.KO_Genes), dwt, dko, obj = obj)
     
    return chromosome(c_and, fitnes_c_and, nwt, nko), chromosome(c_or, fitnes_c_or, nwt, nko), chromosome(c_xor, fitnes_c_xor, nwt, nko)
