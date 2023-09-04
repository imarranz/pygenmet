import numpy as np
import pandas as pd

from ..tools.tools import chromosome

def _replacement(type = 'rws'):
    
    """\
    Parameters
    ----------
    type : Type of replacement
           rws : Replace Worst Strategy
           rts : Restricted Tournament Selection
           wams : Worst Among Most Similar Replacement

    Notes
    -----    
    The _replacement function is ised to change the old generation population 
    with the new generation. Four process have been evaluated [1].
    
    1) Replace Worst Strategy (RWS). The worst element of the population is 
    replaced if the child improves it. It offers high selective pressure, 
    even when its parents are chosen randomly.
    
    2) Restricted Tournament Selection (RTS) [2].
    
    3) Replace the Worst Among Most Similar Replacement (WAMS) [3]. 
    The worst chromosome of the set of the N (N = 3,...) parents most similar 
    to the offspring.
    
    4) Deterministic Crowding Algorithm (Deterministic Crowding, DC) [4]. 
    To facilitate the comparison we will use in our experiments. 
    A variant of the DC in which a single descendant will be generated for 
    each cross, which will replace the most similar parent if it improves it.  
    
    References
    ----------
    [1] https://sci2s.ugr.es/keel/pdf/keel/congreso/4-diversidadfinal2_daniel_molina.pdf
    [2] G.   Harik.   Finding   multimodal   solutions   using restricted  tournament  selection. Proc.  6th  Int.  Conf. Genetic Algorithms, páginas 24-31, 1995.
    [3] W. Cedeño and V. Vemuri. Multi-niche crowding in genetic algorithms and its application to the assembly of dna restriction-fragments. Evolutionary Computation, 2(4):321-345, 1995.    
    [4] S.W. Mahfoud. Crowding and preselection revised. Parallel Problem Solving from Nature 2, páginas 27-36, 1992.
    """
    
    return 5
    

def _selection(chromosomes, size = 1, N = 6, trace = False, method = 'random'):
    
    """\
    THERE ARE TWO TYPES OF SELECTION. THE ONE THAT SELECTS FUTURE PARENTS 
    AND SELECTS THEM PARENTS TO CROSS THEM
    
    MATCHING TECHNIQUES:
    
    Parents can be selected in a way that maintains the diversity of 
    the population
    
    1) Prohibition of crossing based on ancestry. An individual cannot 
    mate with himself, or with his parents, or with his children, 
    or with his siblings.
    
    2) Prohibition of incest. Two parents intersect if their Hamming distance 
    is above a certain threshold
    
    3) Varied pairing. One individual crosses paths with another who is 
    quite different. Hamming distance
    
    THE ONE THAT SELECTS A NEW GENERATION:
    
    1) Random Selection (RS)
    
    2) Tournament Selection (TS): choose the best individual
    fitness among N randomly selected individuals (N = 2, 3,...)
    The selection by tournament constitutes a very extended selection 
    procedure for parents and in which. The idea is to randomly pick a number 
    of individuals from the population, tournament size,
    (with or without replacement), select the best individual from this group, 
    and repeat the process until the number of individuals selected matches 
    the size of the population. Usually the tournament size is 2, and in this 
    case a probabilistic version has been used in which it is allowed
    the selection of individuals without necessarily being the best.

    3) Linear Rank Selection (LRS): the population is ordered based on 
    its fitness and a selection probability is associated with each 
    individual that depends on their order.
    
    4) Roulette Selection (RS): assigns a selection probability proportional 
    to the individual's fitness value Baker (1987) introduces a method called 
    universal stochastic sampling, which uses a single turn of the roulette 
    wheel being the circular sectors proportional to the objective function. 
    The individuals are selected from markers, equally spaced and with 
    a random start. (genetic algorithms.pdf)
    
    5) Elitist: In the elitist selection model, the best individual 
    in the population is forced to at time t, be selected as the parent. 
    
    Parameters
    ----------
    chromosomes : List of chromosomes of the population
    size : Number of individuals to select from the population. 
           This parameter is not considered with method = 'elitist'.
    N : Number of Genes to print (if trace is True)
    trace : Should be the chromosomes and the selection printed on the screen
    method: The method used for selection: 
            ['random', 'tournament', 'linear', 'roulette', 'elitist']
    
    Returns
    -------
    fitness : The fitness of each chromosome (individual)
    psel : Probabilities of selecction of each chromosome. 
           pesl(x) = f(x) / sum(f(y)), for each x in y
    selecction : List of selected chromosomes
    
    Notes
    -----
    The aim of this function is to evaluate the fitness of each population 
    (by the chromosome structure) and to asign a probability to select 
    for a new population (next generation).
    
    Examples
    --------
    >>> A_ = []
    >>> A_.append(chromosome(genes = '01010001010110100100', nwt = 10, nko = 10, fitness = 2.3454))
    >>> A_.append(chromosome(genes = '01011011110110100100', nwt = 10, nko = 10, fitness = 1.4401))
    >>> A_.append(chromosome(genes = '01010001010100000100', nwt = 10, nko = 10, fitness = 0.9254))
    >>> A_.append(chromosome(genes = '01110101010110100100', nwt = 10, nko = 10, fitness = 7.1104))
    >>> A_.append(chromosome(genes = '01010101110111100100', nwt = 10, nko = 10, fitness = 2.1494))
    >>> A_.append(chromosome(genes = '01010001110000100100', nwt = 10, nko = 10, fitness = 1.9954))
    >>> A_.append(chromosome(genes = '00000000010110100100', nwt = 10, nko = 10, fitness = 1.7323))
    >>> A_.append(chromosome(genes = '01010000111110100111', nwt = 10, nko = 10, fitness = 1.5002))
    >>> A_.append(chromosome(genes = '01110011010110100111', nwt = 10, nko = 10, fitness = 2.4119))
    >>> A_.append(chromosome(genes = '01010011011110101100', nwt = 10, nko = 10, fitness = 2.1414))

    >>> selection, fit, ps, sel_mean, sel_maximun = _selection(A_, size = 4, N = 8, trace = True, method = 'roulette')

    |---------------------------------------------------------------------------------------------|
    |-------------- fenotypes ---------------||--- genotypes ---||--- fitness ---||- Probability -|
    01010001...01000101 || 01101001...10100100         4        4           2.3454       0.09912204
    01011011...01101111 || 01101001...10100100         7        4           1.4401       0.07244166
    01010001...01000101 || 01000001...00000100         4        2           0.9254       0.05727278
    01110101...11010101 || 01101001...10100100         6        4           7.1104       0.23955286
    01010101...01010111 || 01111001...11100100         6        5           2.1494       0.09334567
    01010001...01000111 || 00001001...00100100         5        2           1.9954       0.08880708
    00000000...00000001 || 01101001...10100100         1        4           1.7323       0.08105318
    01010000...01000011 || 11101001...10100111         4        7           1.5002       0.07421289
    01110011...11001101 || 01101001...10100111         6        6           2.4119       0.10108189
    01010011...01001101 || 11101011...10101100         5        6           2.1414       0.09310990
    |---------------------------------------------------------------------------------------------|
    |-------------- fenotypes ---------------||--- genotypes ---||--- fitness ---||- Probability -|
    01110101...11010101 || 01101001...10100100         6        4           7.1104       0.23955286
    01010001...01000101 || 01101001...10100100         4        4           2.3454       0.09912204
    01010001...01000111 || 00001001...00100100         5        2           1.9954       0.08880708
    01010011...01001101 || 11101011...10101100         5        6           2.1414       0.09310990

    >>> selection, fit, ps, sel_mean, sel_maximun = _selection(A_, size = 4, N = 8, trace = True, method = 'linear')

    |---------------------------------------------------------------------------------------------|
    |-------------- fenotypes ---------------||--- genotypes ---||--- fitness ---||- Probability -|
    01010001...01000101 || 01101001...10100100         4        4           2.3454       0.15555555
    01011011...01101111 || 01101001...10100100         7        4           1.4401       0.02222222
    01010001...01000101 || 01000001...00000100         4        2           0.9254              0.0
    01110101...11010101 || 01101001...10100100         6        4           7.1104              0.2
    01010101...01010111 || 01111001...11100100         6        5           2.1494       0.13333333
    01010001...01000111 || 00001001...00100100         5        2           1.9954       0.08888888
    00000000...00000001 || 01101001...10100100         1        4           1.7323       0.06666666
    01010000...01000011 || 11101001...10100111         4        7           1.5002       0.04444444
    01110011...11001101 || 01101001...10100111         6        6           2.4119       0.17777777
    01010011...01001101 || 11101011...10101100         5        6           2.1414       0.11111111
    |---------------------------------------------------------------------------------------------|
    |-------------- fenotypes ---------------||--- genotypes ---||--- fitness ---||- Probability -|
    01011011...01101111 || 01101001...10100100         7        4           1.4401       0.02222222
    01010000...01000011 || 11101001...10100111         4        7           1.5002       0.04444444
    01110101...11010101 || 01101001...10100100         6        4           7.1104              0.2
    01010001...01000101 || 01101001...10100100         4        4           2.3454       0.15555555

    >>> selection, fit, ps, sel_mean, sel_maximun = _selection(A_, size = 4, N = 8, trace = True, method = 'random')

    |---------------------------------------------------------------------------------------------|
    |-------------- fenotypes ---------------||--- genotypes ---||--- fitness ---||- Probability -|
    01010001...01000101 || 01101001...10100100         4        4           2.3454              0.1
    01011011...01101111 || 01101001...10100100         7        4           1.4401              0.1
    01010001...01000101 || 01000001...00000100         4        2           0.9254              0.1
    01110101...11010101 || 01101001...10100100         6        4           7.1104              0.1
    01010101...01010111 || 01111001...11100100         6        5           2.1494              0.1
    01010001...01000111 || 00001001...00100100         5        2           1.9954              0.1
    00000000...00000001 || 01101001...10100100         1        4           1.7323              0.1
    01010000...01000011 || 11101001...10100111         4        7           1.5002              0.1
    01110011...11001101 || 01101001...10100111         6        6           2.4119              0.1
    01010011...01001101 || 11101011...10101100         5        6           2.1414              0.1
    |---------------------------------------------------------------------------------------------|
    |-------------- fenotypes ---------------||--- genotypes ---||--- fitness ---||- Probability -|
    01010101...01010111 || 01111001...11100100         6        5           2.1494              0.1
    01010001...01000101 || 01000001...00000100         4        2           0.9254              0.1
    01010011...01001101 || 11101011...10101100         5        6           2.1414              0.1
    01010001...01000101 || 01101001...10100100         4        4           2.3454              0.1
    
    >>> selection, fit, ps, sel_mean, sel_maximun = _selection(A_, size = 4, N = 8, trace = True, method = 'elitist')

    |---------------------------------------------------------------------------------------------|
    |-------------- fenotypes ---------------||--- genotypes ---||--- fitness ---||- Probability -|
    01010001...01000101 || 01101001...10100100         4        4           2.3454       ----------
    01011011...01101111 || 01101001...10100100         7        4           1.4401       ----------
    01010001...01000101 || 01000001...00000100         4        2           0.9254       ----------
    01110101...11010101 || 01101001...10100100         6        4           7.1104       ----------
    01010101...01010111 || 01111001...11100100         6        5           2.1494       ----------
    01010001...01000111 || 00001001...00100100         5        2           1.9954       ----------
    00000000...00000001 || 01101001...10100100         1        4           1.7323       ----------
    01010000...01000011 || 11101001...10100111         4        7           1.5002       ----------
    01110011...11001101 || 01101001...10100111         6        6           2.4119       ----------
    01010011...01001101 || 11101011...10101100         5        6           2.1414       ----------
    |---------------------------------------------------------------------------------------------|
    |-------------- fenotypes ---------------||--- genotypes ---||--- fitness ---||- Probability -|
    01010101...01010111 || 01111001...11100100         6        5           2.1494       ----------
    01010001...01000101 || 01101001...10100100         4        4           2.3454       ----------
    01110011...11001101 || 01101001...10100111         6        6           2.4119       ----------
    01110101...11010101 || 01101001...10100100         6        4           7.1104       ----------
    
    >>> selection, fit, ps, sel_mean, sel_maximun = _selection(A_, size = 4, N = 8, trace = True, method = 'tournament')

    |---------------------------------------------------------------------------------------------|
    |-------------- fenotypes ---------------||--- genotypes ---||--- fitness ---||- Probability -|
    01010001...01000101 || 01101001...10100100         4        4           2.3454       ----------
    01011011...01101111 || 01101001...10100100         7        4           1.4401       ----------
    01010001...01000101 || 01000001...00000100         4        2           0.9254       ----------
    01110101...11010101 || 01101001...10100100         6        4           7.1104       ----------
    01010101...01010111 || 01111001...11100100         6        5           2.1494       ----------
    01010001...01000111 || 00001001...00100100         5        2           1.9954       ----------
    00000000...00000001 || 01101001...10100100         1        4           1.7323       ----------
    01010000...01000011 || 11101001...10100111         4        7           1.5002       ----------
    01110011...11001101 || 01101001...10100111         6        6           2.4119       ----------
    01010011...01001101 || 11101011...10101100         5        6           2.1414       ----------
    |---------------------------------------------------------------------------------------------|
    |-------------- fenotypes ---------------||--- genotypes ---||--- fitness ---||- Probability -|
    00000000...00000001 || 01101001...10100100         1        4           1.7323       ----------
    01010000...01000011 || 11101001...10100111         4        7           1.5002       ----------
    01110011...11001101 || 01101001...10100111         6        6           2.4119       ----------
    01010101...01010111 || 01111001...11100100         6        5           2.1494       ----------
    01110101...11010101 || 01101001...10100100         6        4           7.1104       ----------


    """
    
    if method == 'roulette':
    
        fitness = [x.Fitness for x in chromosomes]
        # In our analysis, the fitness values can be negative values. For this reason a new value has been performed.
        # The fitness_ value is the fitness value + the minimun value of all values of fitness. This procedure
        # will be zero for minimun value of fitness. A positive value is added and then all values of fitness_ will be
        # positive.
        fitness_ = fitness + np.abs(np.min(fitness)) + np.abs(np.min(fitness))/10 # f(x) min = 0
        psel = fitness_ / np.sum(fitness_)
        selection = np.random.choice(a = np.arange(len(chromosomes)), size = size, p = psel, replace = False)
        
    if method == 'random':
        
        selection = np.random.choice(a = np.arange(len(chromosomes)), size = size, p = None, replace = False)
        fitness = [x.Fitness for x in chromosomes]
        psel = [1/len(chromosomes) for x in np.arange(len(chromosomes))]
        
    if method == 'linear':
        
        fitness = np.array([x.Fitness for x in chromosomes])
        temp = fitness.argsort()
        ranks = np.empty_like(temp)
        ranks[temp] = np.arange(len(fitness))
        psel = ranks / np.sum(ranks)        
        selection = np.random.choice(a = np.arange(len(chromosomes)), size = size, p = psel, replace = False)     
        
    if method == 'tournament':

        fitness = np.array([x.Fitness for x in chromosomes])        
        # Only for print purposes
        psel = ['-'*10 for x in np.arange(len(chromosomes))]

        st = 2 # The size of tournaments is 2 by default        
        tournaments_ = np.random.choice(np.arange(len(chromosomes)), size = len(chromosomes), replace = False)
        tournaments_ = [tournaments_[i:i+st] for i in np.arange(0, len(chromosomes), st)]  
        selection = []
        for t in tournaments_:
            B_ = [chromosomes[i] for i in t]
            selection.append(t[np.argmax([x.Fitness for x in B_])])
        
    if method == 'elitist':

        fitness = np.array([x.Fitness for x in chromosomes])        
        # Only for print purposes
        psel = ['-'*10 for x in np.arange(len(chromosomes))]
        selection = np.argsort(fitness)[-size:]

    
    if trace:
            
        k1 = len("-"*(N*4 + 10)) - len(" fenotypes ")
        
        header_to_print = "|" + "-"*((N*4 + 10) + 51) + "|\n" +\
                          "|" + "-"*((k1 //2)-1) + " fenotypes " + "-"*((k1 - k1 // 2)-1) + "|" +\
                          "|--- genotypes ---|" +\
                          "|--- fitness ---|" +\
                          "|- Probability -|"
        
        print(header_to_print)
        
        for i in np.arange(len(chromosomes)):

            chromosome = chromosomes[i]
            N_WT = sum(1 for x in chromosome.WT_Genes if x == '1')
            N_KO = sum(1 for x in chromosome.KO_Genes if x == '1')

            print("{}...{} || {}...{} {:>9}{:>9}{:>17}{:>17}".format(chromosome.WT_Genes[:N], 
                                                                                   chromosome.WT_Genes[-N:],
                                                                                   chromosome.KO_Genes[:N],
                                                                                   chromosome.KO_Genes[-N:],
                                                                                   str(N_WT), str(N_KO), 
                                                                                   str(chromosome.Fitness)[:10],
                                                                                   str(psel[i])[:10]))
        print(header_to_print)        

        for i in selection:     
            
            chromosome = chromosomes[i]
            N_WT = sum(1 for x in chromosome.WT_Genes if x == '1')
            N_KO = sum(1 for x in chromosome.KO_Genes if x == '1')
            
            print("{}...{} || {}...{} {:>9}{:>9}{:>17}{:>17}".format(chromosome.WT_Genes[:N], 
                                                                                   chromosome.WT_Genes[-N:],
                                                                                   chromosome.KO_Genes[:N],
                                                                                   chromosome.KO_Genes[-N:],
                                                                                   str(N_WT), str(N_KO), 
                                                                                   str(chromosome.Fitness)[:10],
                                                                                   str(psel[i])[:10]))
    
    return [chromosomes[index] for index in selection], fitness, psel, np.min(fitness), np.max(fitness)
 
