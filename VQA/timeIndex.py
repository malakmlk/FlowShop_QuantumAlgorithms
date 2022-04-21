#!/usr/bin/env python
# coding: utf-8

# In[28]:


from qiskit_optimization import QuadraticProgram
from qiskit_optimization.translators import from_docplex_mp
from docplex.mp.model import Model
from qiskit_optimization.converters import QuadraticProgramToQubo
from qiskit.utils import algorithm_globals, QuantumInstance
from qiskit_optimization.algorithms import GurobiOptimizer
from typing import List
from VQASolver import VariationalSolver
from itertools import product


# In[30]:


class FSPTimeIndexform():
    """Quantum Optimization for the FSP tIME Index FORME"""
    def __init__(self,timeSpan : int,numberMachine : int,procTime:List[List[int]],numberJobs : int,approach : int)-> None :
        """
        Args : 
        timeSpan : the makespan value
        numberMachine : machine number
        numverJobs : job's number
        """
        self.timespan = timeSpan
        self.numberMachine = numberMachine
        self.numberJobs = numberJobs
        self.procTime = procTime
        self.approach = approach

     
    def to_quadratic_program_validity(self)-> QuadraticProgram :
       
        mdl = Model(name = "FSP_timeIndexOperation")
        N = self.numberJobs
        M = self.numberMachine
        PM = self.procTime 
        T = self.timespan

        # create binary variable xi_t for operation i starts at t
        x = {(i,m,t):mdl.binary_var(name= f"x_{i}_{m}_{t}") for i in range(N) for m in range(M) for t in range(T)}
       
        # constraint 1 : Only one operation starts at t
        for i in range(N): 
            for j in range(M):
                mdl.add_constraint(mdl.sum(x[(i,j,k)] for k in range(T))==1)
 
        #constraint 2 :  the precedence constraint within a job
        for m in range(M-1):
            for i in range(N):
                for t1 in range(T):
                    b= t1+PM[m][i] if t1+PM[m][i] < T else  T 
                    for t2 in range(b):
                        mdl.add_constraint(x[(i,m,t1)]+x[(i,m+1,t2)]<=1)
        #constraint 3 :  No overlapping under the same machine
        for m,i1,i2 in product(range(M),range(N),range(N)):
            if i1 != i2 :
                for t1 in range(T):
                    b= t1+PM[m][i1] if t1+PM[m][i1] < T else T 
                    for t2 in range(t1,b): 
                        mdl.add_constraint(x[(i1,m,t1)]+x[(i2,m,t2)]<=0)

        op = from_docplex_mp(mdl)    
        return op

    def to_quadratic_program_time_threshold(self)->QuadraticProgram :
        mdl = Model(name = "FSP_timeIndexOperation_TTA")
        N = self.numberJobs
        M = self.numberMachine
        PM = self.procTime #the total number of operations
        T = self.timespan

        # create binary variable xi_t for operation i starts at t
        x = {(i,m,t):mdl.binary_var(name= f"x_{i}_{m}_{t}") for i in range(N) for m in range(M) for t in range(T)}
        y = {(i):mdl.binary_var(name = f"y_{i}") for i in range(T)}
        
        # constraints :
        # constraint 1 : Only one operation starts at t
        for i in range(N): 
            for j in range(M):
                mdl.add_constraint(mdl.sum(x[(i,j,k)] for k in range(T))==1)
 
        #constraint 2 :  the precedence constraint within a job
        for m in range(M-1):
            for i in range(N):
                for t1 in range(T):
                    b= t1+PM[m][i] if t1+PM[m][i] < T else  T 
                    for t2 in range(b):
                        mdl.add_constraint(x[(i,m,t1)]+x[(i,m+1,t2)]<=1)
        #constraint 3 :  No overlapping under the same machine
        for m,i1,i2 in product(range(M),range(N),range(N)):
            if i1 != i2 :
                for t1 in range(T):
                    b= t1+PM[m][i1] if t1+PM[m][i1] < T else T 
                    for t2 in range(t1,b): 
                        mdl.add_constraint(x[(i1,m,t1)]+x[(i2,m,t2)]<=1)
            
        # constraint 4 : Only one slot is considred as timespan
        mdl.add_constraint(mdl.sum(y[i] for i in range(T))==1)
        
        # constraint 5 : Each job is completed by the threshold
        for i,j,m,t in product(range(T),range(N),range(M),range(T) ) : 
            if t > i-PM[m][j] :
                mdl.add_constraint(x[(j,m,t)]+y[i]<=1)    
    
        # Objective function :
        mdl.minimize(mdl.sum(y[i]*i for i in range(T)))
        
        op =from_docplex_mp(mdl)
        return op
    
    def to_quadratic_program_approx(self) -> QuadraticProgram: 

        mdl = Model(name = "FSP_timeIndexOperation")
        N = self.numberJobs
        M = self.numberMachine
        PM = self.procTime 
        T = self.timespan

        # create binary variable xi_t for operation i starts at t
        x = {(i,m,t):mdl.binary_var(name= f"x_{i}_{m}_{t}") for i in range(N) for m in range(M) for t in range(T)}
        
        # constraints :
        # constraint 1 : Only one operation starts at t
        for i in range(N): 
            for j in range(M):
                mdl.add_constraint(mdl.sum(x[(i,j,k)] for k in range(T))==1)
 
        #constraint 2 :  the precedence constraint within a job
        for m in range(M-1):
            for i in range(N):
                for t1 in range(T):
                    b= t1+PM[m][i] if t1+PM[m][i] < T else  T 
                    for t2 in range(b):
                        mdl.add_constraint(x[(i,m,t1)]+x[(i,m+1,t2)]<=1)
        #constraint 3 :  No overlapping under the same machine
        for m,i1,i2 in product(range(M),range(N),range(N)):
            if i1 != i2 :
                for t1 in range(T):
                    b= t1+PM[m][i1] if t1+PM[m][i1] < T else T 
                    for t2 in range(t1,b): 
                        mdl.add_constraint(x[(i1,m,t1)]+x[(i2,m,t2)]<=1)
            
        # Objective function :
        mdl.minimize(mdl.sum(
            (t + PM[M-1][i])*x[(i,M-1,t)] 
            for i in range(N) 
            for t in range(T) )
            )

        op = from_docplex_mp(mdl)    
        return op
    
    def to_qubo(self,approach)->QuadraticProgram:  
        conv = QuadraticProgramToQubo()
        return conv.convert(self.to_quadratic_program_time_threshold())

    def Ising(self) -> QuadraticProgram :
         qubitOp, offset = self.to_qubo(1).to_ising()
         return qubitOp, offset 


        



# In[ ]:




