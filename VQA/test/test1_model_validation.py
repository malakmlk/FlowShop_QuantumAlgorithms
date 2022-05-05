import sys
from typing import Reversible

sys.path.insert(1,"C:\\Users\\malak\\FlowShop_QuantumAlgorithms\\VQA")


from instance_taillard_generator import instance_generator
from fspmodel.timeIndex import FSPTimeIndexform
from fspmodel.FSPasTSP import FSPasTSP
import random


for i in range(20):
    N = random.randint(3,6)
    M = random.randint(2,5)
    pm = instance_generator(N,M,random.randint(1,1),random.randint(2,3))    
    print(N,M,pm)
    T = sum(pm[m][i] for m in range(M) for i in range(N))
    # TimeIndex Test
    fsp = FSPTimeIndexform(T,M,pm,N,2)
    mdl1 = fsp.to_quadratic_program_time_threshold()
    mdl2 = fsp.quadratic_pn([1,1])
    print(T)
    for i in range(9):
        solution = [random.randint(0,1) for i in range(N*M*T + T)]
        print(solution)
        print("total cost")
        print(mdl2.objective.evaluate(solution))
        print("feasability")
        constraint = mdl1.get_feasibility_info(solution)
        print(len(constraint[2]))
        print("objective evaluation")
        obj_eval = 0
        for i in range(T):
            obj_eval += i*solution[N*M*T+i]
        print(obj_eval)
    """
    # TSP approximation Test
    for i in range(1,5):
        fsp = FSPasTSP(M,pm,N,i)
        mdl1 = fsp.quadratic_program()
        qubo2 = fsp.to_QUBO()
        for i in range(30):
            solution = [random.randint(0,1) for i in range(N*N)]
            print(solution)
            print("total_cost")
            print(qubo2.objective.evaluate(solution))
            print("feasability")
            constraint = mdl1.get_feasibility_info(solution)
            print(len(constraint[2]))
            print("objective evaluation")
            print(mdl1.objective.evaluate(solution))
    """
    