#!/usr/bin/env python
# coding: utf-8

# In[1]:


#initialization
from quantastica.qiskit_toaster import ToasterBackend
import matplotlib.pyplot as plt
import math 
get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
from math import sqrt
# importing Qiskit
from qiskit import *
from qiskit.circuit.library.arithmetic import DraperQFTAdder
# import basic plot tools
from qiskit.tools.visualization import plot_histogram
from qiskit.utils import QuantumInstance


# In[2]:


# init quantum register with variable x
#number of qubits
#x the init var
def initRegister(n,x):
    xb=format(x,"b")
    var = QuantumRegister(n)
    qc=QuantumCircuit(var)
    for i in range(n-len(xb),n):
        if xb[i-n+len(xb)] == '1':
            qc.x(var[i])
    U_initRegister = qc.to_gate()
    U_initRegister.name = "initRegister"
    return  U_initRegister         


# In[3]:


'''
if var1 < var2 set output to |1>
'''
def lessComparison(n):
    arr =[]
    var1= QuantumRegister(n)
    var2= QuantumRegister(n)
    ancilla = QuantumRegister(n)
    output = QuantumRegister(1)
    qc =  QuantumCircuit(var1,var2,ancilla,output)
    for i in range(n):
        qc.cx(var1[i],ancilla[i])
        qc.cx(var2[i],ancilla[i])
        arr.append(var2[i])
        for j in range(i+1):
            arr.append(ancilla[j])
        qc.mct(arr,output)
        qc.x(ancilla[i])
        arr = []
    #re-init the ancilla qubits in order to reuse them after
    for i in range(n):
        qc.cx(var1[i],ancilla[i])
        qc.cx(var2[i],ancilla[i])
        arr.append(var2[i])
        for j in range(i+1):
            arr.append(ancilla[j])
        qc.x(ancilla[i])
        arr = []
    U_LessCMP = qc.to_gate()
    U_LessCMP.name = "LessCMP"
    return  U_LessCMP  


# In[4]:


def completionTimeInit(n,q,x) :
    position = QuantumRegister(n)
    completionTime= QuantumRegister(q)
    circuit = QuantumCircuit(position,completionTime)
    N = 2**n
    for i in range(2**n):
        #generate control array for job i 
        xb=format(i,'b')
        for j in range(n-len(xb),n):
            if xb[j-n+len(xb)] == '1':
                circuit.x(position[j])
        arr = [position[i] for i in range(n)]        
        #init the completion time to the value of the processing time Pij
        #circuit.append(initRegister(q,x[N-1-i]).control(n),arr)
        pij=format(x[N-1-i],"b")
        for p in range(q-len(pij),q):
            if pij[p-q+len(pij)] == '1':
                circuit.mct(arr,completionTime[p])
        #re-init the control array to 0
        for j in range(n-len(xb),n):
            if xb[j-n+len(xb)] == '1':
                circuit.x(position[j])
    U_initRegister = circuit.to_gate()
    U_initRegister.name = "init CiMj"
    return  U_initRegister             


# In[5]:


def ciMj(q):
    """
    c1 : c1m1
    c2 ; c2m2
    ancilla
    cmp
    c3 
    carry
    """
    arr=[]
    c1 = QuantumRegister(q)
    c2 = QuantumRegister(q)
    c3 = QuantumRegister(q)
    ancilla = QuantumRegister(q)
    cmp = QuantumRegister(1)
    carry = QuantumRegister(1)
    qc = QuantumCircuit(c1,c2,ancilla,cmp,c3,carry)
    qc.append(lessComparison(q),[i for i in range(3*q+1)])
    arr.append(3*q)
    for i in reversed(range(q)):
        arr.append(i)
    for i in reversed(range(3*q+1,4*q+1)):
        arr.append(i)
    arr.append(carry)
    qc.x(cmp[0])
    qc.append(DraperQFTAdder(q, kind='half').to_gate().control(1),arr)
    qc.x(cmp[0])
    for i in range(q):
        arr[i+1]=2*q-i-1
    qc.append(DraperQFTAdder(q, kind='half').to_gate().control(1),arr)    
    qc.append(lessComparison(q).inverse(),[i for i in range(3*q+1)])
    U_ciMj = qc.to_gate()
    U_ciMj.name = "CompletionTime"
    return  U_ciMj


# In[6]:


def Makespan(q,n,M,x):
    N=2**n
    qubits =[]   
    ancillaIndex = [i for i in range(N*n + 2*q + N*q*(M-1),N*n + 2*q + N*q*(M-1)+q+1)]
    completionIndex = [0 for i in range(q)]
    schedule = QuantumRegister(N*n,'position')
    completionTimeMachine1 = QuantumRegister(2*q,'completionTimeMachine1')
    completionTime = QuantumRegister(q*(M-1)*(N),'completionTime')
    ancilla = QuantumRegister(q+2,'ancilla')
    qc = QuantumCircuit(schedule,completionTimeMachine1,completionTime,ancilla) 
    # Initialization of all the completions times
    positionIndex = [idx for idx in range(n)]
    for l in range(q): 
        completionIndex[l]=N*n+q+l
    qc.append(completionTimeInit(n,q,x[0]),positionIndex + completionIndex)
    for i in range(N):
        positionIndex = [i*n + idx for idx in range(n)]
        for j in range(M-1):
            for l in range(q): 
                completionIndex[l] = N*n + 2*q + (i)*(M-1)*q + (j)*q + l
            qc.append(completionTimeInit(n,q,x[j+1]),positionIndex + completionIndex)

    #compute the value of the completion of first job on each machine 
    #machine 1 is a special case because it starts always at 0 while for other machines it depends on the previous machine

    for i in range(1,M):
        for j in reversed(range(q)):
            qubits.append(N*n + q + (i-1)*q + j)
        for j in reversed(range(q,2*q)):
            qubits.append(N*n + q + (i-1)*q + j)
        qubits.append(N*n + 2*q + q*(M-1)*(N)+q)   
        qc.append(DraperQFTAdder(q, kind='half').to_gate(),qubits)
        qubits =[]   

    #For the rest of the completion times we will apply the function ciMj
    #Cij = Cij + max(Ci-1,j,,Ci,j-1) i>1 j>1
    for i in range(1,N):
        positionIndex = [i*n + idx for idx in range(n)]
        for j in range(M):
            if j == 0 :
                #machine 1 we just need to do addition and init the first q qubits
                for l in range(q): 
                    completionIndex[l]=N*n+l
                qc.append(completionTimeInit(n,q,x[0]),positionIndex + completionIndex)  
                completionIndex.reverse()
                qc.append(DraperQFTAdder(q, kind='half').to_gate(),completionIndex + [i for i in reversed(range(N*n +q,N*n+2*q))] + [N*n + 2*q + q*(M-1)*N])
                completionIndex.reverse()
                qc.append(completionTimeInit(n,q,x[0]).inverse(),positionIndex + completionIndex) 

            else :
                #c1,c2,ancilla,cmp,c3,carry
                if j == 1 :
                    #completion 2 depends on completion 1 which is a special case in term of index
                    qc.append(ciMj(q),[l for l in range(N*n + q, N*n + 2*q)] + [l for l in range(N*n+2*q+(i-1)*(M-1)*q,N*n+2*q+(i-1)*(M-1)*q+q)] + ancillaIndex+ [l for l in range(N*n+2*q+i*(M-1)*q,N*n+2*q+i*(M-1)*q+q)] + [N*n + 2*q + N*q*(M-1)+q+1])
                else :
                    qc.append(ciMj(q), [l for l in range(N*n+2*q+(i-1)*(M-1)*q+(j-1)*q,N*n+2*q+(i-1)*(M-1)*q+(j-1)*q+q)] + [l for l in range(N*n+2*q+i*(M-1)*q+(j-2)*q,N*n+2*q+i*(M-1)*q+(j-2)*q+q)]  + ancillaIndex + [l for l in range(N*n+2*q+i*(M-1)*q+(j-1)*q,N*n+2*q+i*(M-1)*q+(j-1)*q+q)]+[N*n + 2*q + N*q*(M-1)+q+1])
        
    U_makespan = qc.to_gate()
    U_makespan.name = "makespan"
    return  U_makespan     


# In[7]:




# In[8]:


qc.draw()


# # Oracle 1
# 

# In[9]:


def Oracle1(q,n,M,x,upBound):
    N=2**n
    schedule = QuantumRegister(n*N)
    completionTime = QuantumRegister(N*(M-1)*q+2*q)
    ancilla = QuantumRegister(q+2)
    upper= QuantumRegister(q)
    validity = QuantumRegister(1)
    qc = QuantumCircuit(schedule,completionTime,ancilla,upper,validity)
    qc.append(Makespan(q,n,M,x),[i for i in range(n*N+N*(M-1)*q+3*q+2)])
    qc.append(initRegister(q,upBound),[i for i in range(n*N+N*(M-1)*q+3*q+2,n*N+N*(M-1)*q+4*q+2)])
    qc.append(lessComparison(q),[i for i in range(n*N+N*(M-1)*q+q,n*N+N*(M-1)*q+2*q)]+[i for i in range(n*N+N*(M-1)*q+3*q+2,n*N+N*(M-1)*q+4*q+2)]+[i for i in range(n*N+N*(M-1)*q+2*q,n*N+N*(M-1)*q+3*q)]+[n*N+N*(M-1)*q+4*q+2])
    U_Oracle1 = qc.to_gate()
    U_Oracle1.name = "Oracle1"
    return  U_Oracle1         


# In[10]:





# In[11]:


"""
array = [i for i in range(N)]+[i for i in range(2+q*N*M-4,2+q*N*M)]
qc.measure(array,cout)
# Use  qasm_simulator
simulator = Aer.get_backend('qasm_simulator')

# Execute the circuit on the qasm simulator
job = execute(qc, simulator, shots=1000)

# Grab results from the job
result = job.result()

# Return counts
counts = result.get_counts(qc)
plot_histogram(counts, figsize=(16, 16),title="flowshop")
"""


# In[12]:





# In[13]:




# # Estimation de l'éspace de mémoire nécessaire 
# 

# In[14]:


def memoryEstimationOracle1(q,n,M):
    N=2**n
    requiredQubits=N*n+2*q+q*(M-1)*(N)+2*q+3
    requiredMemorySpace = 2**requiredQubits *32/(1024*1024*1024)
    return [requiredQubits,requiredMemorySpace]


# In[15]:


n=1
M=2
q=4
memoryEstimationOracle1(q,n,M)


# ---------------------------------------------------------------------------------------------------------------------

# In[16]:


"Oracle 1 compasants"
'''
N=2**n
schedule = QuantumRegister(n*N)
completionTime = QuantumRegister(N*(M-1)*q+2*q)
ancilla = QuantumRegister(q+2)
upper= QuantumRegister(q)
validity = QuantumRegister(1)
qc = QuantumCircuit(schedule,completionTime,ancilla,upper,validity)
qc.append(Makespan(q,n,M,x),[i for i in range(n*N+N*(M-1)*q+3*q+2)])
qc.append(initRegister(q,11),[i for i in range(n*N+N*(M-1)*q+3*q+2,n*N+N*(M-1)*q+4*q+2)])
qc.append(lessComparison(q),[i for i in range(n*N+N*(M-1)*q+q,n*N+N*(M-1)*q+2*q)]+[i for i in range(n*N+N*(M-1)*q+3*q+2,n*N+N*(M-1)*q+4*q+2)]+[i for i in range(n*N+N*(M-1)*q+2*q,n*N+N*(M-1)*q+3*q)]+[n*N+N*(M-1)*q+4*q+2])
qc.draw('mpl')   
'''


# In[17]:


"Completion time"
"""
n = 2
N = 2**n
M = 2
q=4
M1=[2,1,1,3]
M2=[3,1,2,2]
schedule = QuantumRegister(n*N)
c1 = QuantumRegister(q)
c2 = QuantumRegister(q)
c3 = QuantumRegister(q)
c4 = QuantumRegister(q)
ancilla = QuantumRegister(q)
cmp = QuantumRegister(1)
carry = QuantumRegister(1)
cout = ClassicalRegister(n*N+q)
qc = QuantumCircuit(schedule,c1,c2,c3,ancilla,cmp,carry,cout)
qc.h(schedule)

#c1M2
qc.append(completionTimeInit(n,q,M2),[0,1]+[i for i in range(N,N+q)])
#c2M1
qc.append(completionTimeInit(n,q,M1),[2,3]+[i for i in range(N+q,N+2*q)])
#c2M2
qc.append(completionTimeInit(n,q,M2),[2,3]+[i for i in range(N+2*q,N+3*q)])

qc.barrier()


qc.measure([i for i in range(n*N)] + [i for i in range(N+2*q,N+3*q)]  ,cout) 
IBMQ.load_account()
       
my_provider = IBMQ.get_provider(hub='ibm-q', group='open', project='main')
device = my_provider.get_backend("simulator_mps") 
quantum_instance = QuantumInstance(device, shots =1000,skip_qobj_validation = False)
tqc = transpile(qc,device)
result = quantum_instance.execute(tqc)   
state = result.get_counts(tqc)
shots = quantum_instance.run_config.shots
results = {key[::-1]: val / shots for key, val in sorted(state.items()) if val > 0}
circuit_results = {b: (v / shots) ** 0.5 for (b, v) in state.items()} 
results
"""


# In[18]:



'''
n=1
N=2**n
q=4
x1 = [2,3,2,4,4,5,6,7,1,9,10,11,12,13,14,7,15,3,2,3,4,5,6,7,1,9,10,11,12,13,14,7]
x2 = [5,6,2,4,4,5,6,7,1,9,10,11,12,13,14,7,15,3,2,3,4,5,6,7,1,9,10,11,12,13,14,7]
position1 = QuantumRegister(n,'p1')
position2 = QuantumRegister(n,'p2')
completionTime1= QuantumRegister(q,'cp')
completionTime2= QuantumRegister(q,'cp2')
completionTime3= QuantumRegister(q)
completionTime4= QuantumRegister(q)
carry = QuantumRegister(1,'carry')
cout=ClassicalRegister(N+2*q)
circuit = QuantumCircuit(position1,position2,completionTime1,completionTime2,completionTime3,completionTime4,carry,cout)
circuit.h(position1)
circuit.h(position2)
circuit.append(completionTimeInit(n,q,x1),[0]+ [i for i in range(2,2+q)])
circuit.append(completionTimeInit(n,q,x1),[1] + [i for i in range(2+q,2+2*q)])
circuit.append(completionTimeInit(n,q,x2),[0] + [i for i in range(2+2*q,2+3*q)])
circuit.append(completionTimeInit(n,q,x2),[1] + [i for i in range(2+3*q,2+4*q)])
circuit.append(DraperQFTAdder(q, kind='half'), [i for i in reversed(range(2,2+2*q))]+[2+4*q])
#circuit.append(completionTimeInit(n,q,x).inverse(),[i for i in range(n+4*q)])
circuit.draw()
'''


# In[19]:


''''
n=1
N=2**n
q=4
x = [1,2,3,4,4,5,6,7,1,9,10,11,12,13,14,7,15,3,2,3,4,5,6,7,1,9,10,11,12,13,14,7]
position1 = QuantumRegister(n,'p1')
position2 = QuantumRegister(n,'p2')
completionTime1= QuantumRegister(q,'cp')
completionTime2= QuantumRegister(q,'cp2')
cpt= QuantumRegister(q)
carry = QuantumRegister(1,'carry')
cout=ClassicalRegister(N+2*q)
circuit = QuantumCircuit(position1,position2,completionTime1,cpt,completionTime2,carry,cout)
#circuit.h(position1)
#circuit.h(position2)
for i in range(2**n):
        #generate control array for job i 
        xb=format(i,'b')
        for j in range(n-len(xb),n):
            if xb[j-n+len(xb)] == '1':
                circuit.x(position1[j])
        arr = [position1[i] for i in range(n)]
        #init the completion time to the value of the processing time Pij
        #circuit.append(initRegister(q,x[N-1-i]).control(n),arr)  
        pij=format(x[N-1-i],"b")
        for p in range(q-len(pij),q):
            if pij[p-q+len(pij)] == '1':
                circuit.mct(arr,completionTime1[p])
        #re-init the control array to 0
        for j in range(n-len(xb),n):
            if xb[j-n+len(xb)] == '1':
                circuit.x(position1[j])                
        circuit.barrier()                
circuit.barrier()
for i in range(2**n):
        #generate control array for job i 
        xb=format(i,'b')
        for j in range(n-len(xb),n):
            if xb[j-n+len(xb)] == '1':
                circuit.x(position1[j])
        arr = [position1[i] for i in range(n)]
        #init the completion time to the value of the processing time Pij
        #circuit.append(initRegister(q,x[N-1-i]).control(n),arr)  
        pij=format(x[N-1-i],"b")
        for p in range(q-len(pij),q):
            if pij[p-q+len(pij)] == '1':
                circuit.mct(arr,cpt[p])
        #re-init the control array to 0
        for j in range(n-len(xb),n):
            if xb[j-n+len(xb)] == '1':
                circuit.x(position1[j])                
        circuit.barrier()                
circuit.barrier()
for i in range(2**n):
        #generate control array for job i 
        xb=format(i,'b')
        for j in range(n-len(xb),n):
            if xb[j-n+len(xb)] == '1':
                circuit.x(position2[j])
        arr = [position2[i] for i in range(n)]
        #init the completion time to the value of the processing time Pij
        #circuit.append(initRegister(q,x[N-1-i]).control(n),arr)  
        pij=format(x[N-1-i],"b")
        for p in range(q-len(pij),q):
            if pij[p-q+len(pij)] == '1':
                circuit.mct(arr,completionTime2[p])
        #re-init the control array to 0
        for j in range(n-len(xb),n):
            if xb[j-n+len(xb)] == '1':
                circuit.x(position2[j])                
        circuit.barrier()  
#circuit.x(completionTime2[p])
circuit.append(DraperQFTAdder(q, kind='half'), [i for i in range(N,N+2*q)]+[N+2*q]+ )
#circuit.append(completionTimeInit(n,q,x).inverse(),[i for i in range(n+q)])
circuit.draw()
'''


# In[ ]:





# In[20]:


"""
nqubits=N*n+2*q+q*(M-1)*(N)+2*q+2
measurement_rc = ClassicalRegister(n*N + q)
qc.add_register(measurement_rc)
qc.measure([i for i in range(n*N)] + [i for i in range(n*N +q+q*(M-1)*(N) ,n*N+2*q+q*(M-1)*(N))]  ,measurement_rc) 
qc.draw()
IBMQ.load_account()
       
my_provider = IBMQ.get_provider(hub='ibm-q', group='open', project='main')
device = my_provider.get_backend("simulator_mps") 
quantum_instance = QuantumInstance(device, shots =1000,skip_qobj_validation = False)
tqc = transpile(qc,device)
result = quantum_instance.execute(tqc)   
state = result.get_counts(tqc)
shots = quantum_instance.run_config.shots
results = {key[::-1]: val / shots for key, val in sorted(state.items()) if val > 0}
circuit_results = {b: (v / shots) ** 0.5 for (b, v) in state.items()} 
results

"""

