import numpy as np
import pandas as pd

from ..tools.tools import chromosome

def _mutation(father, geneSet, get_fitness, nwt, nko, dwt, dko, obj, type = 'bsm', p = None):        

    """\
    Parameters
    ----------
    father : Genes of the father
    geneSet : Set of genes availables. Usually '0' and '1'
    get_fitness : Function to improve
    nwt : Number of No NAFLD Samples
    nko : Number of NAFLD Samples
    dwt : Data of No NAFLD Samples
    dko : Data of NAFLD Samples
    obj : The fold-changes of murine model to adjust
    type : Type of mutation process, 'bsm' by default. Five types are allowed.
           bf : bitflip
           bsm : bitstringmutation
           sw : swap
           in : inversion
           sc : scramble
    p : Parameter to modify the posibility of change in a genes
    
    Returns
    -------
    A new mutated chromosome by the type of mutation process
    
    Notes
    -----
    The flipbit (fb) mutation is an aggressive mutation because change all genes. 
    If in a specific gen the alelo value is '0' this process will change it to '1'.
    
    fb: '0101010101' -> '1010101010'
    
    The bitstringmutation (bsm) mutation process will change a gen with a 
    probability 'p'. If p = None then p = 1/len(genes) [1].
    
    bsm: '0101010101' -> '0101110101'
    
    References
    ----------
    
    [1] https://en.wikipedia.org/wiki/Mutation_(genetic_algorithm)
    [2] https://www.tutorialspoint.com/genetic_algorithms/genetic_algorithms_mutation.htm
    
    Examples
    --------
    
    """
    
    son = list(father)

    if type == 'fb': # flipbit
        
        for i in np.arange(0, len(father)):
            son[i] = '0' if son[i] == '1' else '1'
    
    if type == 'bsm': # 'bitstringmutation'
      
        mutation_probability = np.random.rand(len(father))
        p = (1 / len(father)) if p == None else p
        
        for i in np.arange(0, len(father)):
            if mutation_probability[i] < p:
                son[i] = '0' if son[i] == '1' else '1'
                
    if type == 'sw': # swap
        p1, p2 = np.sort(np.random.choice(a = np.arange(len(father)), size = 2, replace = False))
        son[p1], son[p2] = son[p2], son[p1]
        
    if type == 'in': # inversion
        p1, p2 = np.sort(np.random.choice(a = np.arange(len(father)), size = 2, replace = False))
        son[p1:p2] = reversed(son[p1:p2])
        
    if type == 'sc': # scramble
        p1, p2 = np.sort(np.random.choice(a = np.arange(10), size = 2, replace = False))
        son[p1:p2] = np.random.choice(a = son[p1:p2], size = len(son[p1:p2]), replace = False)
        
                
    genes = ''.join(son)
    fitness = get_fitness(genes, nwt, nko, dwt, dko, obj = obj)
    
    return chromosome(genes, fitness, nwt, nko)
