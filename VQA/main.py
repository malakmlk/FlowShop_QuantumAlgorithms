
from FSPasTSP import FSPasTSP
from timeIndex import FSPTimeIndexform
from VQASolver import VariationalSolver
from qiskit.providers.aer import QasmSimulator
from qiskit.utils import  QuantumInstance
from qiskit.algorithms.optimizers import SPSA
from qiskit import IBMQ
from qiskit_optimization.algorithms import GurobiOptimizer
import sys

from qiskit.providers.aer import QasmSimulator
device = QuantumInstance(QasmSimulator(method='matrix_product_state'), shots=300)
vqa_Solver = VariationalSolver(1,1,'C:\\Users\\malak\\FlowShop_QuantumAlgorithms\\VQA\\instance.txt')
instances = vqa_Solver.read_Data()

#file_path = 'C:\\Users\\malak\\FlowShop_QuantumAlgorithms\\VQA\\analyse.txt'
#sys.stdout = open(file_path, "a+")



options = {
	'backend_name': 'simulator_mps'
}
spsa = SPSA(maxiter=1000)
for key,value in instances.items():
    for i in range(1,3):
        fsp = FSPasTSP(value[0][1],value[1],value[0][0],i)
        op,off = fsp.to_Ising()
        vqa_result = vqa_Solver.VQE_run_time_service(value[0][0],op,options,spsa,"random")  
        vqa_result = vqa_Solver.VQE(op,device)
        print(value)
        print(vqa_result)
        fsp = FSPTimeIndexform(8,value[0][1],value[1],value[0][0])
        op,off = fsp.Ising(i)
        vqa_result = vqa_Solver.VQE_run_time_service(value[0][0],op,options,spsa,"random")  
        vqa_result = vqa_Solver.VQE(op,device)
        print(value)
        print(vqa_result)
        print('-------------------------')


"""
# QUBO VISUALIZATION 
for key,value in instances.items():
    #create time index object
    for i in range(1,4):
        T = sum(value[1][i][j] for i in range(value[0][1]) for j in range(value[0][0]) )
        fsp = FSPTimeIndexform(T,value[0][1],value[1],value[0][0])
        qubo = fsp.to_qubo(i)
        print(qubo)      
        print(GurobiOptimizer().solve(qubo))

"""