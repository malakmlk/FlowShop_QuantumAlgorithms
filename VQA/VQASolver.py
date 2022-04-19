from qiskit.algorithms import VQE
from typing import List
from qiskit import IBMQ
from qiskit.algorithms.optimizers import SPSA
from qiskit.circuit.library import TwoLocal
# useful additional packages
import networkx as nx
from qiskit.utils import algorithm_globals, QuantumInstance
from qiskit import Aer




class VariationalSolver():
    def __init__(self, 
        forme : int,
        method : int , 
        quantumInstance :  QuantumInstance , 
        data : str,            
        )->None:
        
        """
         Different approach to solve QUBO
            Forme : TimeIndex form, Tsp Forme, PositionIndex
            method : Which method to be used VQE,QAOA,XY-QAOA,...
            quantumInstance : Backend to execute the circuit
            data : Instance pathe

        """
        self.forme = forme,
        self.method = method
        self.quantumInstance = quantumInstance
        self.data = data
        
    def read_Data(self) ->dict[int,(int,int,[[]])] :
        Instance,l,M = {},[],[]
        i,p,c = 0,0,''
        file = open(self.data,"r")
        ins = file.read().split()
        while i in range(len(ins)):
            n = int(ins[i])
            m = int(ins[i+1])
            for k in range(i+2,i+m+2):
                for j in range(len(ins[k])) :
                    if ins[k][j] == ',' : 
                        l.append(int(c))
                        c=''
                    else : c += ins[k][j] 
                M.append(l)
                l= []
            Instance[p]=(n,m),M
            M = []
            i = i + m +2 
            p +=1
        return Instance
      
    def VQE(self,operator, device):
        #operator,offeset = problem.to_Ising()
        ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
        spsa = SPSA(maxiter=1000)
        vqe = VQE(ansatz, optimizer=spsa, quantum_instance=device)
        result = vqe.compute_minimum_eigenvalue(operator)
        print(result)
        optimizer_evals = result.optimizer_evals
        return optimizer_evals
    
    def QAOA(self,operator,device):
        
        return 0