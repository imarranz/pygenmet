import numpy as np
import pandas as pd
import random
##from datetime.datetime import now
import datetime
# xor function is loading to evaluate the hamming distance
from operator import xor
import json
    
from ..fitness.fitness import *

# CLASS CHROMOSOME

class chromosome:

    """\
    The class chromosome has a copy of its genes, its fitness and 
    the number of No NAFLD samples and NAFLD samples for print purposes.
    
    Attributes
    ----------
    
    Methods
    -------
    
    print()  
      Prints the chromosome
    
    to_json()   
        
    Notes
    -----
    The representation of a chromosome is the following way
        
    0110101101 ... 0110101101 || 0110110011 ... 0110110011     6     6           0.8945        
    -------------------------    -------------------------   ----- -----         ------
                A                            B                 C     D              E
                    
    Description of the representation:
      A] Wild Type Genes
      B] Knock-Out Genes
      C] Number of Wild Type cases selected 
      D] Number of Knock-out cases selected
      E] The fitness value. The value of the adaptation for this chromosome
        
    Examples
    --------
    >>> CH = chromosome('01101011010110110011', 0.8945, 10, 10)
    >>> print("Genes: {}".format(CH.Genes))
    Genes: 01101011010110110011        
    >>> print("WT Genes: {}".format(CH.WT_Genes))
    WT Genes: 0110101101        
    >>> print("KO Genes: {}".format(CH.KO_Genes))
    KO Genes: 0110110011        
    >>> print("Fitness: {}".format(CH.Fitness))
    Fitness: 0.8945        
    >>> CH.print() 
    0110101101 ... 0110101101 || 0110110011 ... 0110110011 	   6	   6	         0.8945        
    """

    # https://realpython.com/documenting-python-code/
    
    
    def __init__(self, genes, fitness, nwt, nko):
    
        self.Genes = genes
        self.Fitness = fitness
        self.WT_Genes = genes[:nwt]
        self.KO_Genes = genes[nwt:]    
        
    def print(self, N = 15):
        
        """\
        The N parameter is only for print purposes. Its value is 15 by default.
        """
    
        N_WT = sum(1 for x in self.WT_Genes if x == '1')
        N_KO = sum(1 for x in self.KO_Genes if x == '1')
    
        print("{} ... {} || {} ... {} \t{:>4}\t{:>4}\t{:>15}".format(self.WT_Genes[:N], 
                                                                     self.WT_Genes[-N:], 
                                                                     self.KO_Genes[:N], 
                                                                     self.KO_Genes[-N:],
                                                                     N_WT, N_KO, self.Fitness))      
        
    def to_json(self, path_or_buf):
        
        """\
        Export the structure of chromosome to json file
        """
        
        chromosome_to_json = {}
        chromosome_to_json['chromosome'] = []
        chromosome_to_json['chromosome'].append({
            'Genes': self.Genes, 
            'Fitness': self.Fitness, 
            'WT_Genes': self.WT_Genes,
            'KO_Genes': self.KO_Genes})
        
        with open(path_or_buf, 'w') as outfile:
            json.dump(chromosome_to_json, outfile)



def _read_chromosome(path_or_buf):   
    
    """\
    Parameters
    ----------
    path_or_buf :  string or file handle. File path or object.
    
    Returns
    -------
    This function returns a chromosome
    
    Examples
    --------
    >>> A = chromosome(genes = '00100101001001011000111010110010', nko=16, nwt=16, fitness=0.98)
    >>> A.print(5)
    00100 ... 00101 || 10001 ... 10010 	   6	   8	           0.98
    >>> A.to_json('prueba.json')
    >>> B = _read_chromosome('prueba.json')
    >>> B.print(10)
    0010010100 ... 0100100101 || 1000111010 ... 1010110010 	   6	   8	           0.98    
    """
    
    with open(path_or_buf) as json_file:
        data = json.load(json_file)
        for p in data['chromosome']:   
            ch = chromosome(genes = p['Genes'], 
                            fitness = p['Fitness'],
                            nwt = len(p['WT_Genes']),
                            nko = len(p['KO_Genes']))
    return ch    

def _get_father(nwt, nko, dwt, dko, geneSet, get_fitness, obj, random_state = None):

    """\
    Parameters
    ----------
    nwt : Number of observations of Wild Type
    nko : Number of observations of Knock-out
    dwt : Data of observations of Wild Type
    dko : Data of observations of Knock-out
    geneSet: Usually a list as '01'
    get_fitness: The fitness function to use
    obj : Array of fold-changes of metabolites in murine model
    random_state: If you don't mention the random_state in the code, then whenever you execute your code a new random value

    Returns
    -------
    This function returns a chromosome
    
    Examples
    --------    
    >>> _get_father(nwt = nwt, nko = nko, dwt = dwt, dko = dko, geneSet = '01', get_fitness = _get_fitness, obj = FC.FC).print()
    011101101010000 ... 000000001000001 || 101110011001010 ... 011110110000011 	  37	 261	-0.04568152517704738
    
    >>> _get_father(nwt = nwt, nko = nko, dwt = dwt, dko = dko, geneSet = '01', get_fitness = _get_fitness, obj = FC.FC).print()
    000101011001011 ... 101100010101101 || 010000101000100 ... 110111010001010 	  46	 252	-0.07192848665245967
    
    >>> _get_father(nwt = nwt, nko = nko, dwt = dwt, dko = dko, geneSet = '01', get_fitness = _get_fitness, obj = FC.FC).print()
    110101100111011 ... 000110110110111 || 110010000110100 ... 110010011101011 	  49	 267	-0.11718322152992591
    """
    
    genes = []
    
    if random_state is None:
    
        genes = random.choices(geneSet, k = nwt + nko)
        genes = ''.join(genes)
        fitness = get_fitness(genes, nwt, nko, dwt, dko, obj)
        
    else:
        
        random.seed(random_state)
        genes = random.choices(geneSet, k = nwt + nko)
        genes = ''.join(genes)
        fitness = get_fitness(genes, nwt, nko, dwt, dko, obj)
    
    return chromosome(genes, fitness, nwt, nko)




def _get_population(individuals, nwt, nko, dwt, dko, geneSet, get_fitness, obj, random_state = None):
    
    """\
    Parameters
    ----------
    individuals : Number of individual of first generation or initial population
    nwt : Number of observations of Wild Type
    nko : Number of observations of Knock-out
    dwt : Data of observations of Wild Type
    dko : Data of observations of Knock-out
    geneSet : Usually a list as '01'
    get_fitness : The fitness function to use
    obj : Array of fold-changes of metabolites in murine model
    random_state : If you don't mention the random_state in the code, then whenever you execute your code a new random value
    
    Returns
    -------
    population_ : A list of the initial population (chromosomes) to initialize the algorithm
    
    Notes
    -----
    
    Examples
    --------
    >>> nwt = 15
    >>> nko = 25
    >>> dwt = pd.DataFrame(np.random.rand(nwt,10))
    >>> dko = pd.DataFrame(np.random.rand(nko,10))
    >>> A = _get_population(individuals = 10, 
                            nwt = nwt, nko = nko, dwt = dwt, dko = dko, 
                            geneSet = '01', get_fitness = _get_fitness, obj = [1,0,1,-1,1,0,0,1,0.5, 0.5])

    >>> [A_.print() for A_ in A];
    110101101011111 ... 110101101011111 || 101110000011110 ... 111101111110101 	  11	  16	0.27010320047631353
    110001101000111 ... 110001101000111 || 111011001000100 ... 001001100010011 	   8	  12	0.4153593539526977
    110110110000111 ... 110110110000111 || 111001000001001 ... 010011110110101 	   9	  13	0.24120026823063395
    001001101100010 ... 001001101100010 || 101000001100111 ... 001110111010010 	   6	  12	0.5976992970173265
    110110101000111 ... 110110101000111 || 100011000000110 ... 001100011111011 	   9	  12	0.4327563879291974
    010110111011011 ... 010110111011011 || 101101110101111 ... 011110010010000 	  10	  13	0.19109063637942147
    011011111001110 ... 011011111001110 || 000001110011000 ... 110001111010111 	  10	  13	0.19291545113399683
    100010011100001 ... 100010011100001 || 011101010010100 ... 101001111100011 	   6	  14	0.3838922708272582
    110001111111010 ... 110001111111010 || 110111010001010 ... 010100000000100 	  10	   9	0.04833875255418793
    001111011101100 ... 001111011101100 || 100011000100101 ... 001010101000110 	   9	  10	0.6180057110908423    
    """
    
    population_ = []
    
    if random_state is None:
    
        for i in np.arange(individuals):
            genes = []
            genes = random.choices(geneSet, k = nwt + nko)
            genes = ''.join(genes)
            fitness = get_fitness(genes, nwt, nko, dwt, dko, obj)
            population_.append(chromosome(genes, fitness, nwt, nko))
    
    else:
        
        random.seed(random_state)
        for i in np.arange(individuals):
            genes = []
            genes = random.choices(geneSet, k = nwt + nko)
            genes = ''.join(genes)
            fitness = get_fitness(genes, nwt, nko, dwt, dko, obj)
            population_.append(chromosome(genes, fitness, nwt, nko))
    
    return population_

def _show_partial_solution(chromosome, starttime = None, N = 15):

    """\
    Parameters
    ----------
    chromosome : Chromosome to show
    starttime : This parameter show the initial time which the process 
                has been initialized. If 'None' the function will assign a 
                specific time.
    N : Number of genes to print. Only for print purposes.
    
    Returns
    -------
    Print in console or a file the partial solution (or a specific solution)
    
    Notes
    -----
    The aim of this function is to print a partial solution found, for examples, 
    after a mutation o crossover process. This solution allowscto trace the solutions found by the algorithm.
    
    Examples
    --------    
    >>> ch = chromosome(genes = '01010001010110100100', nko = 10, nwt = 10, fitness = 2.3454)
    >>>> show_partial_solution(ch, starttime = None)
    0101000101...0101000101 || 0110100100...0110100100    4    4    2.3454    time:      0.0
    >>>> ch = chromosome(genes = '01010001010110100100', nko = 10, nwt = 10, fitness = 2.3454)
    >>>> starttime = datetime.datetime.now()
    >>>> show_partial_solution(ch, starttime = starttime)
    0101000101...0101000101 || 0110100100...0110100100    4    4    2.3454    time:    1.688
    """

    if starttime is None:
        starttime = datetime.datetime.now()
    
    difference = (datetime.datetime.now() - starttime).total_seconds()
    
    N_WT = sum(1 for x in chromosome.WT_Genes if x == '1')
    N_KO = sum(1 for x in chromosome.KO_Genes if x == '1')
    
    print("{}...{} || {}...{} \t{:>4}\t{:>4}\t{:>21}\t time: {:>8}".format(chromosome.WT_Genes[:N], 
                                                                           chromosome.WT_Genes[-N:],
                                                                           chromosome.KO_Genes[:N],
                                                                           chromosome.KO_Genes[-N:],
                                                                           N_WT, N_KO, 
                                                                           chromosome.Fitness, 
                                                                           difference))

def _read_metabolic_model_data(con, murine_model, number_metabolites_model = 50, criterion = 'fc'):
    
    """\
    Parameters
    ----------
    con : SQLAlchemy connectable(engine/connection), database string URI,
        or sqlite3 DBAPI2 connection using SQLAlchemy makes it 
        possible to use any DB supported by that library
    murine_model : Murine model to select data in the Data Bases. For example: 'MAT1A', 'GNMT', ..
    number_metabolites_model : The best N (default is 50) metabolites to separate Wild Type and Knock-Out in murine model
    criterion : Criterion for metabolites selection. There are two options, by fold-change or p-value.
    
    The aim of this function is evaluate and select the best metabolites by fold-change criterion 
    (select the meabolites with the highest fold-change) or by p-value (select the meabolites with the lowest p-values).
    Moreover, this function will select those metabolites analyzed in human serum.
        
    Returns
    -------
    dmice : The data of the murine model 
    dwt : The data of Wild Type observations (human cohort)
    dko : The data of Knock-Out observations (human cohort)
    FC : A dataframe with the fold-change observed in the murine model.
    
    Notes
    -----
    This function evaluates a murine model and return the fold-change observed in the murine model and
    the data of the best 50 metabolites in No NAFLD observations and the data in NAFLD observations.
    
    Inside the definition of the function, the 535 samples of previous analysis [1], [2] and [3] are selected. 
    
    [1] Alonso, C., Fernández-Ramos, D., Varela-Rey, M., Martínez-Arranz, I., Navasa, N., Van Liempd, S. M., … Mato, J. M. (2017). 
        Metabolomic Identification of Subtypes of Nonalcoholic Steatohepatitis. Gastroenterology. 
        https://doi.org/10.1053/j.gastro.2017.01.015
    [2] Iruarrizaga-Lejarreta, M., Varela-Rey, M., Fernández-Ramos, D., Martínez-Arranz, I., Delgado, T. C., Simon, J., … Mato, J. M. (2017). 
        Role of aramchol in steatohepatitis and fibrosis in mice. Hepatology Communications. 
        https://doi.org/10.1002/hep4.1107
    [3] Morrison, M. C., Verschuren, L., Salic, K., Verheij, J., Menke, A., Wielinga, P. Y., … Kleemann, R. (2018). 
        Obeticholic Acid Modulates Serum Metabolites and Gene Signatures Characteristic of Human NASH and Attenuates Inflammation and Fibrosis Progression in Ldlr-/-.Leiden Mice. 
        Hepatology Communications, 2(12), 1513–1532. https://doi.org/10.1002/hep4.1270
    
    Examples
    --------
    """
    
    KO_human_samples = ('OL-09_0001','OL-09_0002','OL-09_0003','OL-09_0004','OL-09_0005','OL-09_0006','OL-09_0007','OL-09_0008','OL-09_0009','OL-09_0011',
                        'OL-09_0012','OL-09_0013',
          'OL-09_0014','OL-09_0015','OL-09_0016','OL-09_0017','OL-09_0018','OL-09_0019','OL-09_0020','OL-09_0021','OL-09_0022','OL-09_0023','OL-09_0024','OL-09_0027',
          'OL-09_0029','OL-09_0030','OL-09_0031','OL-09_0032','OL-09_0033','OL-09_0034','OL-09_0035','OL-09_0036','OL-09_0037','OL-09_0038','OL-09_0039','OL-09_0040',
          'OL-09_0041','OL-09_0042','OL-09_0043','OL-09_0044','OL-09_0045','OL-09_0046','OL-09_0047','OL-09_0048','OL-09_0095','OL-09_0096','OL-09_0097','OL-09_0098',
          'OL-09_0099','OL-09_0100','OL-09_0101','OL-09_0102','OL-09_0103','OL-09_0105','OL-09_0106','OL-09_0107','OL-09_0108','OL-09_0109','OL-09_0110','OL-09_0111',
          'OL-09_0112','OL-09_0113','OL-09_0115','OL-09_0116','OL-09_0118','OL-09_0119','OL-09_0121','OL-09_0122','OL-09_0123','OL-09_0124','OL-09_0127','OL-09_0128',
          'OL-09_0129','OL-09_0130','OL-09_0131','OL-09_0132','OL-09_0133','OL-09_0134','OL-09_0135','OL-09_0136','OL-09_0137','OL-09_0139','OL-09_0140','OL-09_0141',
          'OL-09_0142','OL-09_0191','OL-09_0192','OL-09_0193','OL-09_0194','OL-09_0195','OL-09_0196','OL-09_0197','OL-09_0198','OL-09_0199','OL-09_0200','OL-09_0201',
          'OL-09_0202','OL-09_0203','OL-09_0204','OL-09_0205','OL-09_0206','OL-09_0207','OL-09_0208','OL-09_0209','OL-09_0210','OL-09_0211','OL-09_0212','OL-09_0213',
          'OL-09_0215','OL-09_0216','OL-09_0217','OL-09_0218','OL-09_0220','OL-09_0221','OL-09_0222','OL-09_0223','OL-09_0224','OL-09_0225','OL-09_0226','OL-09_0229',
          'OL-09_0230','OL-09_0231','OL-09_0233','OL-09_0234','OL-09_0235','OL-09_0288','OL-09_0289','OL-09_0290','OL-09_0291','OL-09_0292','OL-09_0293','OL-09_0294',
          'OL-09_0295','OL-09_0296','OL-09_0299','OL-09_0300','OL-09_0301','OL-09_0302','OL-09_0303','OL-09_0304','OL-09_0305','OL-09_0306','OL-09_0307','OL-09_0308',
          'OL-09_0309','OL-09_0310','OL-09_0311','OL-09_0312','OL-09_0313','OL-09_0314','OL-09_0316','OL-09_0320','OL-09_0323','OL-09_0324','OL-09_0325','OL-09_0326',
          'OL-09_0327','OL-09_0332','OL-09_0334','OL-09_0335','OL-09_0336','OL-09_0337','OL-09_0338','OL-09_0339','OL-09_0340','OL-09_0383','OL-10_0001','OL-10_0002',
          'OL-10_0003','OL-10_0004','OL-10_0005','OL-10_0013','OL-10_0034','OL-10_0035','OL-10_0036','OL-10_0037','OL-10_0038','OL-10_0039','OL-10_0040','OL-10_0041',
          'OL-10_0043','OL-10_0044','OL-10_0045','OL-10_0046','OL-10_0047','OL-10_0048','OL-10_0049','OL-10_0050','OL-10_0051','OL-10_0053','OL-10_0054','OL-10_0055',
          'OL-10_0056','OL-10_0057','OL-10_0058','OL-10_0059','OL-10_0060','OL-10_0061','OL-10_0062','OL-10_0063','OL-10_0064','OL-10_0065','OL-10_0066','OL-10_0067',
          'OL-10_0068','OL-10_0069','OL-10_0071','OL-10_0072','OL-10_0074','OL-10_0075','OL-10_0076','OL-10_0077','OL-10_0078','OL-10_0079','OL-10_0080','OL-10_0081',
          'OL-10_0082','OL-10_0083','OL-10_0084','OL-10_0085','OL-10_0087','OL-10_0088','OL-10_0089','OL-10_0090','OL-15-0064','OL-15-0066','OL-15-0068','OL-15-0074',
          'OL-15-0076','OL-15-0078','OL-15-0080','OL-15-0084','OL-15-0086','OL-15-0090','OL-15-0096','OL-15-0100','OL-15-0102','OL-15-0108','OL-15-0122','OL-15-0124',
          'OL-15-0285','OL-15-0286','OL-15-0287','OL-15-0288','OL-15-0290','OL-15-0291','OL-15-0292','OL-15-0293','OL-15-0294','OL-15-0297','OL-15-0299','OL-15-0300',
          'OL-15-0301','OL-15-0303','OL-15-0304','OL-15-0305','OL-15-0306','OL-15-0307','OL-15-0308','OL-15-0309','OL-15-0310','OL-15-0314','OL-15-0315','OL-15-0316',
          'OL-15-0318','OL-15-0319','OL-15-0321','OL-15-0323','OL-15-0324','OL-15-0326','OL-15-0327','OL-15-0329','OL-15-0332','OL-15-0334','OL-15-0337','OL-15-0338',
          'OL-15-0339','OL-15-0340','OL-15-0341','OL-15-0342','OL-15-0343','OL-15-0345','OL-15-0346','OL-15-0348','OL-15-0350','OL-15-0351','OL-15-0352','OL-15-0353',
          'Owl-CBH-001','Owl-CBH-002','Owl-CBH-003','Owl-CBH-006','Owl-CBH-008','Owl-CBH-009','Owl-CBH-013','Owl-CBH-014','Owl-CBH-015','Owl-CBH-016','Owl-CBH-020',
          'Owl-CBH-023','Owl-CBH-024','Owl-CBH-027','Owl-CBH-028','Owl-CBH-032','Owl-CBH-033','Owl-CBH-041','Owl-CBH-155','Owl-CBH-156','Owl-CBH-157','Owl-CBH-158',
          'Owl-CBH-159','Owl-CBH-160','Owl-CBH-161','Owl-CBH-162','Owl-CBH-163','Owl-CBH-164','Owl-CBH-165','Owl-CBH-166','Owl-CBH-167','Owl-CBH-168','Owl-CBH-170',
          'Owl-CBH-171','Owl-CBH-172','Owl-CBH-173','Owl-CBH-174','Owl-CBH-175','Owl-CBH-176','Owl-CBH-177','Owl-CBH-178','Owl-CBH-179','Owl-CBH-180','Owl-CBH-181',
          'Owl-CBH-182','Owl-CBH-183','Owl-CBH-185','Owl-CBH-186','Owl-CBH-187','Owl-CBH-188','Owl-CBH-189','Owl-CBH-192','Owl-CBH-193','Owl-CBH-196','Owl-CBH-198',
          'Owl-CBH-200','Owl-CBH-203','Owl-CBH-205','Owl-CBH-206','Owl-CBH-207','Owl-CBH-210','Owl-CBH-212','Owl-CBH-215','Owl-CBH-216','Owl-CBH-219','Owl-CBH-220',
          'Owl-CBH-221','Owl-CBH-223','Owl-CBH-224','Owl-CBH-225','Owl-CBH-226','Owl-CBH-227','Owl-CBH-228','Owl-CBH-230','Owl-CBH-231','Owl-CBH-232','Owl-CBH-234',
          'Owl-CBH-235','Owl-CBH-236','Owl-CBH-237','Owl-CBH-238','Owl-CBH-241','Owl-CBH-242','Owl-CBH-243','Owl-CBH-244','Owl-CBH-245','Owl-CBH-246','Owl-CBH-247',
          'Owl-CBH-248','Owl-CBH-249','Owl-CBH-250','Owl-CBH-251','Owl-CBH-252','Owl-CBH-253','Owl-CBH-254','Owl-CBH-255','Owl-CBH-256','Owl-CBH-258','Owl-CBH-261',
          'Owl-CBH-262','Owl-CBH-282','Owl-CBH-283','Owl-CBH-284','Owl-CBH-285','Owl-CBH-286','Owl-CBH-287','Owl-CBH-288','Owl-CBH-289','Owl-CBH-290','Owl-CBH-291',
          'Owl-CBH-292','Owl-CBH-293','Owl-CBH-294','Owl-CBH-295','Owl-CBH-296','Owl-CBH-297','Owl-CBH-298','Owl-CBH-299','Owl-CBH-300','Owl-CBH-301','Owl-CBH-302',
          'Owl-CBH-303','Owl-CBH-304','Owl-CBH-305','Owl-CBH-306','Owl-CBH-307','Owl-CBH-308','Owl-CBH-309','Owl-CBH-310','Owl-CBH-311','Owl-CBH-312','Owl-MDC-047',
          'Owl-MDC-197','Owl-MDC-199','Owl-MDC-202','Owl-MDC-208','Owl-MDC-211','Owl-MDC-214','Owl-MDC-215','Owl-MDC-218','Owl-MDC-221','Owl-MDC-224','Owl-MDC-227',
          'Owl-MDC-230','Owl-MDC-233','Owl-MDC-236','Owl-MDC-239','Owl-MDC-242','Owl-MDC-245','Owl-MDC-248','Owl-MDC-251','Owl-MDC-254','Owl-MDC-259','Owl-MDC-265',
          'Owl-MDC-268','Owl-MDC-271','Owl-MDC-274','Owl-MDC-286','Owl-MDC-289','Owl-MDC-292','Owl-MDC-295','Owl-MDC-469','Owl-MDC-472','Owl-MDC-475','Owl-MDC-478',
          'Owl-MDC-481','Owl-MDC-485','Owl-MDC-489','Owl-MDC-495','Owl-MDC-498','Owl-MDC-501','Owl-MDC-504','Owl-MDC-508','Owl-MDC-511','Owl-MDC-513','Owl-MDC-517',
          'Owl-MDC-520','Owl-MDC-523','Owl-MDC-526','Owl-MDC-528','Owl-MDC-536','Owl-MDC-539','Owl-MDC-543','Owl-MDC-545','Owl-MDC-548','Owl-MDC-551','Owl-MDC-554',
          'Owl-MDC-557','Owl-MDC-563','Owl-MDC-566','Owl-MDC-569','Owl-MDC-572','Owl-MDC-575','Owl-MDC-578','Owl-MDC-580','Owl-MDC-583','Owl-MDC-588','Owl-MDC-592',
          'Owl-MDC-595','Owl-MDC-601','Owl-MDC-604','Owl-MDC-611','Owl-MDC-616','Owl-MDC-624','Owl-MDC-627','Owl-MDC-630','Owl-MDC-633','Owl-MDC-641','Owl-MDC-646',
          'Owl-MDC-649','Owl-MDC-655','Owl-MDC-658','Owl-MDC-663','Owl-MDC-675','Owl-MDC-681','Owl-MDC-684','Owl-MDC-690','Owl-MDC-693','Owl-MDC-698','Owl-MDC-703',
          'Owl-MDC-705','Owl-MDC-715','Owl-MDC-720','Owl-MDC-723','Owl-MDC-726','S-CR-3-27','S-E-1-398','S-E-1-412','S-E-1-433','S-E-1-47','S-E-6-330','S-E-6-331',
          'S-E-6-332','S-E-6-333','S-E-6-334','S-E-N-1-387','S-E-N-1-393','S-E-N-1-394','S-E-N-1-402','S-E-N-1-406','S-N-1-50','S-N-1-51','S-N-6-335','S-N-6-336',
          'S-N-6-337','S-N-6-338','S-N-6-339')
    
    # Randomize samples
    
    np.random.seed(7)
    KO_human_samples = tuple(np.random.choice(a = KO_human_samples, 
                                              size = len(KO_human_samples), 
                                              replace = False))

    WT_human_samples = ('OL-09_0237','Owl-CBH-045','OL-09_0028','OL-09_0050','OL-09_0051','OL-09_0117','OL-09_0238','OL-09_0239',
                        'OL-09_0214','Owl-CBH-208','Owl-CBH-213','OL-09_0317','OL-09_0318','OL-09_0322','OL-09_0328','OL-09_0329',
                        'OL-09_0330','OL-09_0331','OL-09_0333','Owl-CBH-004','Owl-CBH-005','Owl-CBH-007','Owl-CBH-011','Owl-CBH-012',
                        'Owl-CBH-017','Owl-CBH-018','Owl-CBH-019','Owl-CBH-021','Owl-CBH-025','Owl-CBH-026','Owl-CBH-029','Owl-CBH-030',
                        'Owl-CBH-034','Owl-CBH-035','Owl-CBH-036','Owl-CBH-037','Owl-CBH-038','Owl-CBH-039','Owl-CBH-040','Owl-CBH-042',
                        'Owl-CBH-043','Owl-CBH-044','OL-09_0114','OL-09_0297','OL-10_0070','OL-10_0073','OL-10_0086','Owl-CBH-010','Owl-CBH-022',
                        'Owl-CBH-031','OL-09_0010','OL-09_0025','OL-09_0026','OL-09_0104','OL-09_0125','OL-09_0138','OL-09_0219','OL-09_0227',
                        'OL-09_0228','OL-09_0232','OL-09_0236','OL-09_0298','OL-10_0042','OL-10_0052','Owl-CBH-191','Owl-CBH-194','Owl-CBH-195',
                        'Owl-CBH-197','Owl-CBH-199','Owl-CBH-201','Owl-CBH-202','Owl-CBH-204','Owl-CBH-209','Owl-CBH-211','Owl-CBH-214','Owl-CBH-217',
                        'Owl-CBH-218','Owl-CBH-222','Owl-CBH-229','Owl-CBH-233','Owl-CBH-239','Owl-CBH-240','Owl-CBH-257','Owl-CBH-259','Owl-CBH-260',
                        'S-C-O-6-325','S-C-O-6-326','S-C-O-6-327','S-C-O-6-328','S-C-O-6-329')
    
    dmice = pd.read_sql_query(sql = 
        f"""
        --==== =======================================================
        --==== ANÁLISIS CON LOS DATOS DE LOS MODELOS MURINOS
        SELECT [ms].[COD.OWL], [dt].[metabolite.id], [ms].[species], [ms].[matrix], 
               [ms].[global.murine.model], [ms].[murine.model], [ms].[class],
               [dt].[intensity], [lg].[simple.name], [lg].[shape], [lg].[size]
          FROM ms
          LEFT JOIN dt ON [ms].[COD.OWL] = [dt].[COD.OWL]
          LEFT JOIN vl ON [vl].[metabolite.id] = [dt].[metabolite.id]
          LEFT JOIN lg ON [lg].[simple.name] = [vl].[simple.name]
        --==== SELECCIONAMOS SÓLO UN MODELO
         WHERE [ms].[global.murine.model] IN ('{murine_model}')
               AND
        --==== QUE TENGA TIPO KNOCK-OUT Y WILD TYPE       
               [ms].[class] IN ("Wild Type", "Knock-out")
               AND
        --==== SÓLO SERUM
               [ms].[matrix] IN ("Serum")
        --==== FIN DE LA QUERY
        --==== =======================================================
        """, con = con)    
    
    from scipy.stats import ttest_ind

    def pvalue(a, b):
        _, p = ttest_ind(a, b)
        p = -np.log10(p)
        return p
    
    dmice = dmice\
        .groupby(['global.murine.model', 'metabolite.id', 'simple.name', 'matrix', 'class'])['intensity']\
        .apply(list).unstack('class')\
        .assign(**{'lk': lambda df: df['Knock-out'].str.len(), 
                  'lw': lambda df: df['Wild Type'].str.len()})\
        .query("lk > 1 & lw > 1")\
        .assign(medialk = lambda df: list(map(np.mean, df['Knock-out'].tolist())))\
        .assign(medialw = lambda df: list(map(np.mean, df['Wild Type'].tolist())))\
        .assign(fold_change = lambda df: np.log2(df['medialk'] / df['medialw']))\
        .assign(pvalue = lambda df: list(map(pvalue, df['Knock-out'].tolist(), df['Wild Type'].tolist())))\
        .drop(labels = ['Knock-out', 'Wild Type', 'medialk', 'medialw', 'lk', 'lw'], axis = 1).copy()

    FC = pd.DataFrame({'Metabolites': dmice.reset_index()['metabolite.id'], 
                       'FC': dmice.reset_index()['fold_change']})
    
    dwt = pd.read_sql_query(sql = 
        f"""
        --==== ============================================
        --==== START THE QUERY
        --==== ============================================
        SELECT * FROM dt
         WHERE [dt].[COD.OWL] IN {WT_human_samples}
        --==== ============================================
        """, con = con)

    dwt = dwt.pivot(index = 'COD.OWL', columns = 'metabolite.id', values = 'intensity').copy()
    dwt = dwt.fillna(dwt.mean()).copy()
    dwt.dropna(axis = 1, inplace = True)

    dko = pd.read_sql_query(sql = 
        f"""
        --==== ============================================
        --==== START THE QUERY
        --==== ============================================
        SELECT * FROM dt
         WHERE [dt].[COD.OWL] IN {KO_human_samples}
        --==== ============================================
        """, con = con)

    dko = dko.pivot(index = 'COD.OWL', columns = 'metabolite.id', values = 'intensity').copy()
    dko = dko.fillna(dko.mean()).copy()
    dko.dropna(axis = 1, inplace = True)
    
    paso = True
    N = number_metabolites_model
    NMAX = number_metabolites_model

    while paso:
        
        if criterion == 'fc':
            metabolites = dmice\
                .assign(abs_fold_change = lambda df: np.abs(df['fold_change']))\
                .sort_values(by = 'abs_fold_change', ascending=False)\
                .head(NMAX).reset_index()['metabolite.id']
            metabolites = metabolites.sort_values()
            
        if criterion == 'pv':
            metabolites = dmice\
                .sort_values(by = 'pvalue', ascending=False)\
                .head(NMAX).reset_index()['metabolite.id']
            metabolites = metabolites.sort_values()

        N1 = len(dwt.loc[:, list(set(metabolites) & set(dwt.columns))].reindex().dropna(axis = 1).columns)
        N2 = len(dko.loc[:, list(set(metabolites) & set(dko.columns))].reindex().dropna(axis = 1).columns)

        if (N1 < N) | (N2 < N):
            NMAX += 1
            if NMAX > 500:
                metabolites = list(set(metabolites) & set(dwt.columns) & set(dko.columns))
                print("There are not {} metabolites between murine data and human data. Only {} metabolites will be shown.".format(N, len(metabolites)))
                dwt = dwt.loc[:, metabolites].dropna(axis = 1).copy()
                dko = dko.loc[:, metabolites].dropna(axis = 1).copy()
                paso = False
        else:
            metabolites = list(set(metabolites) & set(dwt.columns) & set(dko.columns))
            dwt = dwt.loc[:, metabolites].dropna(axis = 1).copy()
            dko = dko.loc[:, metabolites].dropna(axis = 1).copy()
            paso = False
    
    AUX = FC.loc[[True if x in dko.columns else False for x in FC['Metabolites']], :].copy()
    AUX.set_index('Metabolites', inplace = True)
    AUX = AUX.loc[np.array(dko.columns), :].copy()
    AUX.reset_index(inplace = True)
    FC = AUX
    
    return dmice, dwt, dko, FC    


def _hamming(chromosome_a, chromosome_b):
    
    """\
    Parameters
    ----------
    chromosome_a : Chromosome A
    chromosome_b : Chromosome B
    
    Returns
    -------
    This function returns the Hamming distance and a incest parameter between two parents.
    The Hamming Distance is defined as the number of different genes. If the Hamming Distance
    is 0 then both parents are identically. 
    
    The incest parameter is defined as a one minus the proportion of different genes between the length of the chromosome. 
    A incest parameter of 1 implies the Hamming Distance is zero and therefore both parents are identically. 
    On the other hand, a incest parameter of 0 implies the Hamming Distance is N (length of chromosome) and therefore 
    both parents are totally different.
    
    This parameter can be used to select parents. For example, two parents with incest <= 0.2 can be crossed but if 
    they have an incest parameter > 0.2 can not be crossed
    
    Notes
    -----
    The goal of this function is to generate a criterion to crossover the parents.
    
    References
    ----------    
    [1] https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.distance.hamming.html#scipy.spatial.distance.hamming
    [2] https://en.wikipedia.org/wiki/Hamming_distance
    
    Examples
    --------    
    >>> ch_a = chromosome(genes = '01010001010110100100', nko = 10, nwt = 10, fitness = 2.3454)
    >>> ch_b = chromosome(genes = '01010001010110100100', nko = 10, nwt = 10, fitness = 2.3454)
    >>> ch_c = chromosome(genes = '11110100111000110101', nko = 10, nwt = 10, fitness = 1.9921)
    >>> _hamming(ch_a, ch_b)
    (0, 1.0)
    >>> _hamming(ch_a, ch_c)
    (10, 0.5)
    """
    
    distance = 0
    a = chromosome_a.Genes
    b = chromosome_b.Genes
    
    distance = 0

    for i in np.arange(len(a)):
        
        ai = list(a)[i]
        bi = list(b)[i]
        distance += xor(int(ai), int(bi))

    incest = 1 - distance / len(a)
        
    return distance, incest    

def _grouping_by(L, criterion = 'random'):

    """\
    Parameters
    ----------
    L: List of chromosomes
    criterion: Criteria for grouping parents
    
    Returns
    -------
    This function returns a list of couples of chromosomes to crossover.
    The couples are generated by hamming distance (if criterion = 'hamming') or 
    are generated randomly (if criterion = 'random')
 
    _grouping : List of list of pairs of chromosomes. These pairs will be the parents for the new generation.
    
    Examples
    --------
    >>> ch_1 = chromosome('01101011010110110011', 0.8945, 10, 10)
    >>> ch_2 = chromosome('01101011010110110011', 0.7211, 10, 10)
    >>> ch_3 = chromosome('01101011010110110011', 0.5114, 10, 10)
    >>> ch_4 = chromosome('01101011010110110011', 0.1295, 10, 10)
    >>> _grouping_by([ch_1, ch_2, ch_3, ch_4], criterion = 'random')
    [[<genetics2.tools.tools.chromosome at 0x7f3dc8ba1ef0>,
      <genetics2.tools.tools.chromosome at 0x7f3dc8ba1940>],
     [<genetics2.tools.tools.chromosome at 0x7f3dc8ba19e8>,
      <genetics2.tools.tools.chromosome at 0x7f3dc8ba1080>]]
    """
    
    if criterion == 'hamming':
    
        _grouping = []
        L_aux = L.copy()

        while len(L_aux) > 0:
            ch1 = L_aux[0]
            maximum_distance = 0
            for ch2 in L_aux[1:]:
                distance, _ = _hamming(ch1, ch2)
                if distance >= maximum_distance:
                    maximum_distance = distance
                    couple = [ch1, ch2]
            _grouping.append(couple)
            position_ch1 = L_aux.index(ch1)
            del L_aux[position_ch1] 
            position_ch2 = L_aux.index(ch2)
            del L_aux[position_ch2] 
            
    if criterion == 'random':
        
        _grouping = []

        aux = np.random.choice(a = np.arange(len(L)), 
                               size = len(L), 
                               replace = False)

        for i in np.arange(len(L))[::2]:
            ch1 = L[aux[i]]
            ch2 = L[aux[i+1]]
            couple = [ch1, ch2]
            _grouping.append(couple)
            
    return _grouping

def _run_genetic_algorithm(initial_population = 10, 
                           con = None, # Connection with DDBB
                           murine_model = 'MAT1A', # 'GNMT', 'TNO', 'MCD', ...
                           number_metabolites_model = 50,
                           size = 10,
                           selecction_criterion = 'roulette', #['random', 'tournament', 'linear', 'roulette', 'elitist']
                           trace = False,
                           N = 15,
                           geneSet = '01',
                           get_fitness = None, # Fitness or adaptation function. Usually '_get_fitness'
                           mutation_criterion = 'bsm', #['bf', 'bsm', 'sw', 'in', 'sc']
                           p = None,
                           grouping_criterion = 'random', #['random', 'hamming']
                           generations_without_change = 100):
                            
    """\
    Parameters
    ----------
    initial_population : A number of a list of chromosomes
                         If it is a number _get_population(individuals = initial_population, nwt, nko, dwt, dko, geneSet, get_fitness, obj)
                         else
    con : DataBase connection
    murine_model : 
    metabolites_model :
    size : Number of individuals to select from the population. This parameter is not considered with method = 'elitist'.
    N : Number of Genes to print (if trace is True)
    trace : Should be the chromosomes and the selection printed on the screen
    selecction_criterion: The method used for selection: ['random', 'tournament', 'linear', 'roulette', 'elitist']
    geneSet : Set of genes availables. Usually '0' and '1'
    get_fitness : Function to improve
    nwt : Number of No NAFLD Samples. From _read_metabolic_model_data(con, murine_model, metabolites_model)
    nko : Number of NAFLD Samples. From _read_metabolic_model_data(con, murine_model, metabolites_model)
    dwt : Data of No NAFLD Samples. From _read_metabolic_model_data(con, murine_model, metabolites_model)
    dko : Data of NAFLD Samples. From _read_metabolic_model_data(con, murine_model, metabolites_model)
    
    obj : The fold-changes of murine model to adjust
    mutation_criterion : Type of mutation process, 'bsm' by default. Five types are allowed.
                         bf : bitflip
                         bsm : bitstringmutation
                         sw : swap
                         in : inversion
                         sc : scramble
    p : Parameter to modify the posibility of change in a genes
    grouping_criterion = The method used for grouping parents: ['random', 'hamming']
    
    Notes
    -----
    
    Results
    -------
    >>> a, b, c, d, e, f = _run_genetic_algorithm(initial_population = 10, 
                                            con = db, 
                                            number_metabolites_model = 50, 
                                            generations_without_change = 500, 
                                            trace = False, 
                                            selecction_criterion = 'linear')
    >>> [x.print() for x in e]
    100000001000001 ... 000000100110100 || 100011111110001 ... 111001111110010 	  21	 230	0.8861236068307543
    100000001000001 ... 000000100110100 || 100011111110001 ... 111001111110010 	  21	 231	0.885523197610578
    100000001000001 ... 000000100110100 || 100011111110001 ... 111001111110010 	  21	 230	0.8861236068307543
    100000001000001 ... 000000100110100 || 100011111110001 ... 111001111110010 	  21	 229	0.8858149441931368
    100000001000001 ... 000000100110100 || 100011111110001 ... 111001111110010 	  21	 230	0.8861236068307543
    100000001000001 ... 000000100110100 || 100011111110001 ... 111001111110010 	  21	 230	0.8861236068307543
    100000001000001 ... 000000100110100 || 100011111110001 ... 111001111110010 	  21	 230	0.8861236068307543
    100000001000001 ... 000000100110100 || 100011111110001 ... 111001111110010 	  21	 230	0.8861236068307543
    100000001000001 ... 000000100110100 || 100011111110001 ... 111001111110010 	  21	 230	0.8861236068307543
    100000001000001 ... 000000100110100 || 100011111110001 ... 111001111110010 	  21	 230	0.8861236068307543
    """
    
    dmice, dwt, dko, FC = _read_metabolic_model_data(con, murine_model, number_metabolites_model = number_metabolites_model, criterion = 'pv')
    
    nwt = dwt.shape[0]
    nko = dko.shape[0]
    
    FC.set_index(['Metabolites'], inplace = True)
    FC.rename(columns = {'FC': murine_model}, inplace = True)
    
    horaInicio = datetime.datetime.now()
    results_to_dataframe = []
    generations = 0
    incremento = 0
    
    initial_population if isinstance(initial_population, int) else len(initial_population)
    
    if isinstance(initial_population, int):
    
        population = _get_population(individuals = initial_population, 
                                     nwt = nwt, 
                                     nko = nko, 
                                     dwt = dwt, 
                                     dko = dko, 
                                     geneSet = geneSet, 
                                     get_fitness = get_fitness, 
                                     obj = np.array(FC.loc[:, murine_model]))
    else:
        
        population = initial_population
        
    while generations < generations_without_change:

        L = _grouping_by(population, criterion = grouping_criterion)
        Next_Generation = []

        for i in np.arange(len(L)):

            A, B = _crossover(L[i][0], L[i][1], 
                              get_fitness = get_fitness,  
                              nwt = nwt, 
                              nko = nko, 
                              dwt = dwt, 
                              dko = dko, 
                              obj = np.array(FC.loc[:, murine_model]), 
                              k = 2)

            A = _mutation(A.Genes, geneSet = geneSet, 
                          get_fitness = get_fitness, 
                          nwt = nwt, 
                          nko = nko, 
                          dwt = dwt, 
                          dko = dko, 
                          obj = np.array(FC.loc[:, murine_model]), 
                          type = mutation_criterion, 
                          p = p)

            B = _mutation(B.Genes, geneSet = geneSet, 
                          get_fitness = get_fitness, 
                          nwt = nwt, 
                          nko = nko, 
                          dwt = dwt, 
                          dko = dko, 
                          obj = np.array(FC.loc[:, murine_model]), 
                          type = mutation_criterion, 
                          p = p)

            Next_Generation.append(A)
            Next_Generation.append(B)

        A_ = np.concatenate([population, Next_Generation])

        population, fit, ps, sel_mean, sel_maximun = _selection(A_, 
                                                                size = initial_population if isinstance(initial_population, int) else len(initial_population), 
                                                                N = 15, 
                                                                trace = False, 
                                                                method = selecction_criterion)   

        if sel_mean > incremento:
            
            results_to_dataframe.append([(datetime.datetime.now() - horaInicio).total_seconds(), generations, sel_mean, sel_maximun])
            if trace:
                print("{}".format([(datetime.datetime.now() - horaInicio).total_seconds(), generations, sel_mean, sel_maximun]))

            incremento = sel_mean
            generations = 0
            
        else:
            
            generations += 1

    results_to_dataframe = pd.DataFrame(results_to_dataframe, columns = ['time', 'generations', 'fitness (mean)', 'fitness (max)'])
    
    return dmice, dwt, dko, FC, population, results_to_dataframe


# End of tools.py
