#!/usr/bin/env python
# coding: utf-8

# In[1]:

#!/usr/bin/env python
# coding: utf-8

# In[1]:


from qiskit_optimization import QuadraticProgram
from qiskit_optimization.translators import from_docplex_mp
from docplex.mp.model import Model
from qiskit_optimization.converters import QuadraticProgramToQubo
from typing import List


# In[15]:


class positionIndex : 
    def __init__(self,numberMachine : int,procTime:List[List[int]],numberJobs : int) :
        """
        Args :
        numberMachine : machine number
        numverJobs : job's number
        """
        self.numberMachine = numberMachine
        self.numberJobs = numberJobs
        self.procTime = procTime
        
    
    def to_quadratic_program_approx(self) :
        mdl = Model("Position based model")
        N = self.numberJobs
        M = self.numberMachine
        PM = self.procTime 

        # create binary variable xi_t for operation i starts at t
        x = {(i,j):mdl.binary_var(name= f"x_{i}_{j}") for i in range(N) for j in range(N)}

        #each position containe one job 
        for j in range(N) :
            mdl.add_constraint (mdl.sum( x[(i,j)] for i in range(N))==1)

        #each job is assigned to one position
        for i in range(N) :
            mdl.add_constraint (mdl.sum( x[(i,j)] for j in range(N))==1)

        # Minimize 
        mdl.minimize(
            mdl.sum(
                mdl.sum(
                    mdl.sum(PM[m][i+1]*x[(i+1,j)] for j in range(N)) - mdl.sum(PM[m+1][i]*x[(i,j)] for j in range(N))
                    for i in range(N-1)
                    )
                for m in range(M-1)   
                )
            )
        op = from_docplex_mp(mdl)
        return op

    def to_qubo(self)->QuadraticProgram:  
        conv = QuadraticProgramToQubo()
        return conv.convert(self.to_quadratic_program_approx())

    def to_ising(self) -> QuadraticProgram :
         qubitOp, offset = self.to_qubo().to_ising()
         return qubitOp, offset 


