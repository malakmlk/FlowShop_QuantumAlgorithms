import sys
from typing import Reversible

sys.path.insert(1,"C:\\Users\\malak\\FlowShop_QuantumAlgorithms\\VQA")

from qiskit_optimization.algorithms import CplexOptimizer,GurobiOptimizer
from instance_taillard_generator import instance_generator
from fspmodel.timeIndex import FSPTimeIndexform
from fspmodel.FSPasTSP import FSPasTSP
import random

for i in range(20):
    N = random.randint(2,20)
    M = random.randint(2,8)
    pm = instance_generator(N,M,random.randint(1,1),random.randint(2,3))    
    print(N,M,pm)
    T = sum(pm[m][i] for m in range(M) for i in range(N))
    """
    print("Time Index")
    # TimeIndex Test
    fsp = FSPTimeIndexform(T,M,pm,N,2)
    mdl1 = fsp.to_quadratic_program_time_threshold()
    mdl2 = fsp.quadratic_pn([1,1])
    print(CplexOptimizer().solve(mdl2))
    """
    print("FSP AS TSP")
    # FSP AS TSP Test
    for i in range(1,5):
        fsp = FSPasTSP(M,pm,N,i)
        mdl1 = fsp.to_QUBO()
        print(CplexOptimizer().solve(mdl1))