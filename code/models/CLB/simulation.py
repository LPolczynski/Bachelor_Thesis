from joblib import Parallel, delayed
from cellclass import *

def simulation(init_number, duration):
    population = []
    old_daughters = [cell('init_species_cyclins','init_species_growth',duration,0) for i in range(init_number)]
    new_daughters = []
    while old_daughters:
        results = Parallel(n_jobs=-1)(delayed(c.simulate)() for c in old_daughters)
        #results = [old_daughters[0].simulate()]
        population.extend([res[0] for res in results])
        daughters = [res[1] for res in results]
        for d_list in daughters:
            new_daughters.extend(d_list)
        print('collected: ',len(new_daughters))
        old_daughters=new_daughters
        new_daughters=[]
    return population