#!/usr/bin/env python
# coding: utf-8

# In[18]:


#import qiskit
from qiskit import QuantumCircuit,ClassicalRegister,QuantumRegister,execute
from qiskit.utils import QuantumInstance
#import Grover operator parts
from Oracle1Generator import Oracle1
from Oracle2Generator import Oracle2
from Diffuser import diffuser
# import backend in order to run the circuit
from qiskit.providers import Backend, BaseBackend
from qiskit import *
from typing import List
# Visualization
import matplotlib.pyplot as plt
import numpy as np


# In[19]:


def memoryEstimationOracle(q,n,M,cmp):
    N=2**n
    requiredQubits=N*n+2*q+q*(M-1)*(N)+2*q+3+n+cmp+2
    requiredMemorySpace = 2**requiredQubits *32/(1024*1024*1024)
    return [requiredQubits,requiredMemorySpace]


# In[20]:


class GroverFlowSop:
    '''
    Grover search algorithm for solving flow shop
    the oracle part of grover is modified in order to distinguish solution
    the same preparation state and diffuser are used
    '''
    def __init__(
        self,
        num_qubits_job : int,
        num_qubits_cmj : int,
        Pm :List[List[int]],
        num_machine : int,
        upperBound : int,
        iterations,
        quantum_instance : QuantumInstance
                 
    )->None :
        '''
        Args :
           - iterations : in order to specify the number of iteration of grover
           - num_qubits_job : nbr of qubits per job.
           - num_qubits_cimj : nbr of qubits per completion time.
           - num_machine : nbr of machines. 
           - Pm : matrix that contain the processing times of all jobs on machines
           - upper_bound : the upper bound.  
           - quantum_instance : A Quantum Instance or Backend to run the circuits.
        Raises :
        ValueError: if iteration = None
        ValueError:if problem = Non
        '''
        if iterations == None :
            raise ValueError("Pass a value of iteration")
      
        self.iterations = iterations
        self._quantum_instance = None
        self.num_qubits_job = num_qubits_job
        self.num_qubits_cimj =num_qubits_cmj
        self.Pm = Pm
        self.num_machine = num_machine
        self.upperBound = upperBound
        if quantum_instance is not None:
            self.quantum_instance = quantum_instance
    
    def GroverOperator(self)-> QuantumCircuit:
        n=self.num_qubits_job
        N=2**(n)
        q=self.num_qubits_cimj
        x=self.Pm
        M=self.num_machine
        up=self.upperBound
        if n != 1 :
            cmp= 2**(n-1)*(N-1)
            nqubits = N*n+2*q+q*(M-1)*(N)+2*q+3+n+cmp+2
            qRegister = QuantumRegister(nqubits)
            qc = QuantumCircuit(qRegister) 
            qc.append(Oracle1(q,n,M,x,up),[i for i in range(N*n+q*(M-1)*(N)+4*q+2)]+[nqubits-3])
            qc.append(Oracle2(N,n),[i for i in range(N*n)]+[i for i in range(N*n+q*(M-1)*(N)+4*q+2,N*n+q*(M-1)*(N)+4*q+2+n+cmp)]+[nqubits-2])
            qc.ccx(qRegister[nqubits-3],qRegister[nqubits-2],qRegister[nqubits-1])
            #reversibility
            qc.append(Oracle1(q,n,M,x,up).inverse(),[i for i in range(N*n+q*(M-1)*(N)+4*q+2)]+[nqubits-3])
            qc.append(Oracle2(N,n).inverse(),[i for i in range(N*n)]+[i for i in range(N*n+q*(M-1)*(N)+4*q+2,N*n+q*(M-1)*(N)+4*q+2+n+cmp)]+[nqubits-2])
        else :  #when N=2 we can reduce the size of the circuit  
            nqubits = N*n+2*q+q*(M-1)*(N)+2*q+3+2
            qRegister = QuantumRegister(nqubits)
            qc = QuantumCircuit(qRegister) 
            qc.append(Oracle1(q,n,M,x,up),[i for i in range(N*n+q*(M-1)*(N)+4*q+2)]+[nqubits-3])
            qc.cx(0,nqubits-2)
            qc.cx(1,nqubits-2)
            qc.ccx(qRegister[nqubits-3],qRegister[nqubits-2],qRegister[nqubits-1])
            #reversibility
            qc.append(Oracle1(q,n,M,x,up).inverse(),[i for i in range(N*n+q*(M-1)*(N)+4*q+2)]+[nqubits-3])
            qc.cx(0,[nqubits-2])
            qc.cx(1,[nqubits-2])
        qc.append(diffuser(n*N),[i for i in range(n*N)])
        return qc
    
    def ConstructCircuit(self)->QuantumCircuit:
            n=self.num_qubits_job
            N=2**(n)
            q=self.num_qubits_cimj
            x=self.Pm
            M=self.num_machine
            if n == 1 :
                nqubits = N*n+2*q+q*(M-1)*(N)+2*q+3+2
                qRegister = QuantumRegister(nqubits)
                qc = QuantumCircuit(qRegister) 
            else :
                cmp= 2**(n-1)*(N-1)
                nqubits = N*n+2*q+q*(M-1)*(N)+2*q+3+n+cmp+2
                qRegister = QuantumRegister(nqubits)
                qc = QuantumCircuit(qRegister)                 
            #state preparation
            qc.h([i for i in range(n*N)])
            qc.x(-1)
            qc.h(-1)
            #Grover operator
            Grover_op= self.GroverOperator()
            for i in range(self.iterations):
                qc.append(Grover_op,qc.qubits)    
            #Measurement     
            measurement_cr = ClassicalRegister(n*N)
            qc.add_register(measurement_cr)
            qc.measure([i for i in range(n*N)], measurement_cr)    
            return qc
    
    def execute(self):
        qc = self.ConstructCircuit()
        result = self.quantum_instance.execute(qc)   
        state = result.get_counts(qc)
        shots = self.quantum_instance.run_config.shots
        hist = {key[::-1]: val / shots for key, val in sorted(state.items()) if val > 0}
        self._circuit_results = {b: (v / shots) ** 0.5 for (b, v) in state.items()} 
        return hist

    def memoryEstimationOracle1(self):
        M=self.num_machine
        n=self.num_qubits_job
        q=self.num_qubits_cimj
        N=2**n
        requiredQubits=N*n+2*q+q*(M-1)*(N)+2*q+3
        requiredMemorySpace = 2**requiredQubits *32/(1024*1024*1024)
        return [requiredQubits,requiredMemorySpace]
        
    #def memoryEstimationOracle2(self):




        

        

