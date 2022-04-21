from random import random
from unittest import result
from qiskit.algorithms import VQE
from typing import List
from qiskit import IBMQ
from qiskit.algorithms.optimizers import SPSA
from qiskit.circuit.library import TwoLocal
# useful additional packages
import networkx as nx
from qiskit.utils import algorithm_globals, QuantumInstance
from qiskit import Aer
from qiskit.circuit.library import TwoLocal




class VariationalSolver():
    def __init__(self, 
        forme : int,
        method : int , 
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
        self.data = data
        
    def read_Data(self)  :
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
    
    def VQE_run_time_service(self,n,operator,options,optimizer,init):
        
        ansatz = TwoLocal(n*n,rotation_blocks='ry', entanglement_blocks='cz')
        runtime_inputs = {
            # A parameterized quantum circuit preparing
        # the ansatz wavefunction for the
        # VQE. It is assumed that
        # all qubits are initially in
        # the 0 state.
        'ansatz': ansatz, # object (required)
            
        # A list or dict (with
        # strings as keys) of operators
        # of type PauliSumOp to be
        # evaluated at the final, optimized
        # state.
        'aux_operators': None, # array
            
        # Initial position of virtual qubits
        # on the physical qubits of
        # the quantum device. Default is
        # None.
        'initial_layout': None, # [null,array,object]
            
        # Initial parameters of the ansatz.
        # Can be an array or
        'initial_parameters': init, # [array,string] (required)
            
        # Whether to apply measurement error
        # mitigation in form of a
        # complete measurement fitter to the
        # measurements. 
        'measurement_error_mitigation': None, # boolean
            
        # The Hamiltonian whose smallest eigenvalue
        # we're trying to find. Should
        # be PauliSumOp
        'operator': operator, # object (required)
            
        # The classical optimizer used in
        # to update the parameters in
        # each iteration. 
        'optimizer': optimizer, # object (required)
            
        # The number of shots used
        # for each circuit evaluation. Defaults
        # to 1024.
        'shots': 500 # integer
        }
        provider = IBMQ.enable_account('f7604df03172d5c0401b28b95b54dc35daa6f593c9be6f18c9d5c58deeea94f82daf492cdfb7d39c496e22fd7e6bc979cefad6517f61b22eeea0cdf355b1c0c6')
        IBMQ.load_account()
        provider = IBMQ.get_provider(
            hub='ibm-q',
            group='open',
            project='main'
        )

        job = provider.runtime.run(
            program_id='vqe',
            options=options,
            inputs=runtime_inputs
        )
    
        return result
        
    def QAOA(self,operator,device):
        
        return 0