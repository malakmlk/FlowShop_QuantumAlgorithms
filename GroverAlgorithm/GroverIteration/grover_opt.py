# %%
import matplotlib.pyplot as plt
import math
import random
from qiskit import IBMQ
from math import sqrt
# importing Qiskit
from qiskit.utils import QuantumInstance, algorithm_globals
# importing backend 
from typing import List
from qiskit.utils import QuantumInstance
# importing Grover 
from GroverFSP import GroverFlowSop
import numpy as np

# %%
class GroverOptimizerFSP :
    def __init__(
        self,
        num_qubit_jobs : int,
        num_machines : int,
        process_time : List[List[int]],
        upper : int,
        quantum_instance : QuantumInstance,
        num_iteration  
        )-> None :
        """
        Args : 
        num_qubit_jobs :number of qubit per job,
        num_machines : nombre of machine ,
        process_time : processing times of jobs under each machine List[List[int]],
        upper : upper bound,
        quantum_instance : QuantumInstance,
        num_iteration  
        """
        self.num_qubit_jobs = num_qubit_jobs
        self.num_machines = num_machines
        self.process_time = process_time
        self.upper = upper
        self.quantum_instance = quantum_instance
        self.num_iteration = num_iteration

    def solve(self):
        """
            GAS solver 
        """
        # Optimum tracking variables
        optimum_found = False
        optimum_permu = math.inf
        optimum_value = math.inf
        trace = []
        # Grover ciruit parameters
        upperBound = self.upper
        n = self.num_qubit_jobs
        N = 2**n
        m = self.num_machines
        q = self.qubits_cmj_Estimation()
        pm = self.process_time       
        quantum_instance = self.quantum_instance
        # solotions tracking
        num_solutions = 2**(N*n)
        schedule_measurement = []
        # Grover result
        grover_output = {}
        # Stop the algorithm if we hit the rotation max
        r = 0 
        max_r = int(math.ceil(100 * math.pi / 4))
        algo_tracking = []
        it = 0
        while not optimum_found :
            r_m = 1
            impovement = False
            it+=1
            print(it)
            # Iterate until we don't get improvement
            nb_no_improvement = 0
            while not impovement :
                # The required number of rotation 
                nb_no_improvement += 1
                r_count = random.randint(1,r_m)
                r += r_count
                print("startGrover")
                # Create Grover FSP Circuit
                grover = GroverFlowSop(n,q,pm,m,upperBound + 1,r,quantum_instance)
                # Execute Grover 
                grover_output = grover.execute()
                print(grover_output)
                # Choose a random solution from grover results
                output = self.get_sol(grover_output)
                schedule = self.convert_solution_int(output,n)
                cmax = self.calculate_makespan(pm,np.array(schedule))
                algo_tracking.append(
                    {
                        "schedule" : schedule,
                        "cmax" : cmax,
                        "rotation" : r,
                        "iteration" : it,
                    }
                )
                print(algo_tracking)
                if (cmax < optimum_value 
                    and self.feasibility_check(schedule)):
                    optimum_value = cmax
                    optimum_permu = schedule
                    impovement = True
                    upperBound = cmax
                    # solution trace 
                    trace.append((optimum_permu,optimum_value))
                else:
                    m = int(np.ceil(min(r_m*8/7,2**(n*N/2))))   
                    if schedule not in schedule_measurement :
                        schedule_measurement.append(schedule)
                    if (nb_no_improvement >= self.num_iteration 
                        or r >= max_r
                        or len(schedule_measurement) == num_solutions):
                        impovement = True
                        optimum_found = True
                    print(nb_no_improvement >= self.num_iteration,r >= max_r,len(schedule_measurement) == num_solutions)

        return (trace,algo_tracking)
                           
    def qubits_cmj_Estimation(self) ->int :
        """
            Estimate the required number of qubits for the circuit
         """
        n = 2**self.num_qubit_jobs
        m = self.num_machines
        return math.ceil(math.log2(sum(self.process_time[i][j] for i in range(m) for j in range(n) )))

    def calculate_makespan(self,a, seq):
        a = np.transpose(a)
        # Order the jobs (rows) in order of the sequence
        a = a[seq]

        b = np.zeros(a.shape)
        jobs = a.shape[0]
        macs = a.shape[1]

        b[0, 0] = a[0, 0]
        # Build first row
        for i in range(1, macs):
            b[0, i] = a[0, i] + b[0, i - 1]
        # Build first column
        for i in range(1, jobs):
            b[i, 0] = a[i, 0] + b[i - 1, 0]
        # Build the rest
        for i in range(1, jobs):
            for j in range(1, macs):
                b[i, j] = a[i, j] + (b[i - 1, j] if b[i - 1, j] > b[i, j - 1] else b[i, j - 1])

        return int(b[-1, -1])
 
    
    def convert_solution_int(self,bin_s,n):
        " Return converted solution "
        s = []
        for i in range(len(bin_s)//n):
            s.append(int(bin_s[i*n : i*n + n],2))
        return s
    
    def get_sol(self,a):  
        return algorithm_globals.random.choice(list(a.keys()), 1, p=list(a.values()))[0]

    def feasibility_check(self,a : List):
        return len(set(a)) == len(a)


    def binary_search_grover(self):
        return 0
        



IBMQ.load_account()
       
my_provider = IBMQ.get_provider(hub='ibm-q', group='open', project='main')
        
device = my_provider.get_backend("simulator_mps") 
        
quantum_instance = QuantumInstance(device, shots =1000,skip_qobj_validation = False)
solver = GroverOptimizerFSP(2,2,[[4,2,2,1],[2,3,1,3]],10,quantum_instance,20)
r = solver.solve()
r

# %%



