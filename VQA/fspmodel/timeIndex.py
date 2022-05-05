#!/usr/bin/env python
# coding: utf-8

# In[28]:


from qiskit_optimization import QuadraticProgram
from qiskit_optimization.translators import from_docplex_mp
from docplex.mp.model import Model
from qiskit_optimization.converters import QuadraticProgramToQubo
from typing import List
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
       
        # Constraint 2 :  The precedence constraint within a job
        expr = 0
        for i,m in product(range(N),range(M-1)):
            expr = 0
            for t1 in range(T):
                b = t1+PM[m][i]  if t1+PM[m][i] < T else T 
                for t2 in range(b):
                    expr +=  x[(i,m,t1)]*x[(i,m+1,t2)]
            mdl.add_constraint(expr == 0)
               
    
        # Constraint 3 : No overlapping under the same machine
        expr = 0
        for m,i1,i2 in product(range(M),range(N),range(N)):
            expr = 0
            if i1 != i2 :
                for t1 in range(T):
                    b= t1+PM[m][i1] if t1+PM[m][i1] < T else T 
                    for t2 in range(t1,b): 
                        expr += x[(i1,m,t1)] * x[(i2,m,t2)]
            mdl.add_constraint(expr == 0)

        op = from_docplex_mp(mdl)    
        return op

    def to_quadratic_program_time_threshold(self)->QuadraticProgram :

        mdl = Model(name = "FSP_timeIndexOperation_TTA")
        N = self.numberJobs
        M = self.numberMachine  
        PM = self.procTime #the total number of operations
        T = self.timespan

        # Create binary variable xi_t for operation i starts at t
        x = {(i,m,t):mdl.binary_var(name= f"x_{i}_{m}_{t}") for i in range(N) for m in range(M) for t in range(T)}
        y = {(i):mdl.binary_var(name = f"y_{i}") for i in range(T)}
        

        # Constraint 1 : Only one operation starts at t
        for i in range(N): 
            for j in range(M):
                mdl.add_constraint(mdl.sum(x[(i,j,k)] for k in range(T))==1,ctname="cstr1_"+str(i)+str(j))
        
        # Constraint 2 :  The precedence constraint within a job
        expr = 0
        for i,m in product(range(N),range(M-1)):
            expr = 0
            for t1 in range(T):
                b = t1+PM[m][i]  if t1+PM[m][i] < T else T 
                for t2 in range(b):
                    expr +=  x[(i,m,t1)]*x[(i,m+1,t2)]
            mdl.add_constraint(expr == 0)
               
        # Constraint 3 : No overlapping under the same machine
        expr = 0
        for m,i1 in product(range(M),range(N)):
            expr = 0
            for i2 in range(N):
                if i1 != i2 :
                    for t1 in range(T):
                        b= t1+PM[m][i1] if t1+PM[m][i1] < T else T 
                        for t2 in range(t1,b): 
                           expr += x[(i1,m,t1)] * x[(i2,m,t2)]
            mdl.add_constraint(expr == 0)
        
        # Constraint 4 : Only one slot is considred as timespan
        mdl.add_constraint(mdl.sum(y[i] for i in range(T))==1)
        
        # Constraint 5 : Each job is completed by the threshold
        for i in range(T) :
            expr = 0
            for j,m in product(range(N),range(M)):
                for t in range(T):
                    if t > i-PM[m][j] :
                        expr += x[(j,m,t)]*y[i]
            mdl.add_constraint(expr == 0)    
         
        # Objective function :
        mdl.minimize(
            mdl.sum(y[i]*i for i in range(T))
            )
        
        op = from_docplex_mp(mdl)
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
        # Constraint 2 :  The precedence constraint within a job
        expr = 0
        for i,m in product(range(N),range(M-1)):
            expr = 0
            for t1 in range(T):
                b = t1+PM[m][i]  if t1+PM[m][i] < T else T 
                for t2 in range(b):
                    expr +=  x[(i,m,t1)]*x[(i,m+1,t2)]
            mdl.add_constraint(expr == 0)
               
    
        # Constraint 3 : No overlapping under the same machine
        expr = 0
        for m,i1,i2 in product(range(M),range(N),range(N)):
            expr = 0
            if i1 != i2 :
                for t1 in range(T):
                    b= t1+PM[m][i1] if t1+PM[m][i1] < T else T 
                    for t2 in range(t1,b): 
                        expr += x[(i1,m,t1)] * x[(i2,m,t2)]
            mdl.add_constraint(expr == 0)
            
        # Objective function :
        mdl.minimize(mdl.sum(
            (t + PM[M-1][i])*x[(i,M-1,t)] 
            for i in range(N) 
            for t in range(T) )
            )

        op = from_docplex_mp(mdl)    
        return op
    

    def quadratic_pn(self,penality : List[int]) -> QuadraticProgram :

        mdl = Model(name = "FSP PENALITY")
        N = self.numberJobs
        M = self.numberMachine
        PM = self.procTime #the total number of operations
        T = self.timespan

        # create binary variable xi_t for operation i starts at t
        x = {(i,m,t):mdl.binary_var(name= f"x_{i}_{m}_{t}") for i in range(N) for m in range(M) for t in range(T)}
        y = {(i):mdl.binary_var(name = f"y_{i}") for i in range(T)}
        
        # Validity expression 
        constraint = 0

        # expr 1 : penalize solution which had operation with more than one starts time
        expr_1 = 0
        for i in range(N):
            for j in range(M):
                expr_1 += (mdl.sum(x[(i,j,k)] for k in range(T))-1)**2

        # expr 2 : the precedence constraint within a job
        expr_2 = 0
        for i,m in product(range(N),range(M-1)):
            expr = 0
            for t1 in range(T):
                b= t1+PM[m][i] if t1+PM[m][i] < T else  T 
                for t2 in range(b):
                    expr += x[(i,m,t1)] * x[(i,m+1,t2)]
            expr_2 += expr
        
        # expr_3 : No overlapping under the same machine
        expr_3 = 0    
        for m,i1 in product(range(M),range(N)):
            for i2 in range(N):
                expr = 0
                if i1 != i2 :
                    for t1 in range(T):
                        b= t1+PM[m][i1] if t1+PM[m][i1] < T else T 
                        for t2 in range(t1,b): 
                            expr += x[(i1,m,t1)] * x[(i2,m,t2)]
                expr_3 += expr

        # expr_4 : Only one slot is considred as timespan
        expr_4 = (mdl.sum(y[i] for i in range(T)) - 1)**2

        # expr_5 : Each job is completed by the threshold
        expr_5 = 0
        for i,j,m,t in product(range(T),range(N),range(M),range(T) ) : 
            if t > i-PM[m][j] :
                expr_5 += x[(j,m,t)]*y[i]  
    
        constraint =  sum(PM[m][i] for m,i in product(range(M),range(N))) * (expr_1 + expr_2 + expr_3 + expr_4 + expr_5)
 
        # Cost expression 
        cost = 0
        for i in range(T) :
            cost += y[i]*i 
        
        # Objectif function 
        mdl.minimize(
             penality[0]*cost + penality[1]*constraint 
        )

        op = from_docplex_mp(mdl)
        return op


    def to_qubo(self,approach)->QuadraticProgram:  
        conv = QuadraticProgramToQubo(10)
        return conv.convert(self.to_quadratic_program_time_threshold())

    def Ising(self) -> QuadraticProgram :
         qubitOp, offset = self.to_qubo(1).to_ising()
         return qubitOp, offset 



# In[ ]:




