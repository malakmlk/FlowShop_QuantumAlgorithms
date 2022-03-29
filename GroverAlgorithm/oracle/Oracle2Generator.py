#!/usr/bin/env python
# coding: utf-8

# In[1]:


#initialization
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
# importing Qiskit
from qiskit import *
from qiskit.circuit.library import ZGate,OR
# import basic plot tools
from qiskit.tools.visualization import plot_histogram


# In[2]:


#n : the required number inoreder to encode a position
def positions_XOR(n):
    position1=QuantumRegister(n)
    position2=QuantumRegister(n)
    ancilla=QuantumRegister(n)
    v=QuantumRegister(1)
    qc = QuantumCircuit(position1,position2,ancilla,v)
    if n == 1 :
        qc.cx(position1[0],v[0])
        qc.cx(position2[0],v[0])
    else :
        for i in range(n):
            qc.cx(position1[i],ancilla[i])
            qc.cx(position2[i],ancilla[i])    
        qc.append(OR(2,[1,1]).to_gate(),[ancilla[0],ancilla[1],v[0]])
        for i in range(n):
            qc.cx(position1[i],ancilla[i])
            qc.cx(position2[i],ancilla[i]) 
    U_positionsCMP = qc.to_gate()
    U_positionsCMP.name = "CMP_position"
    return  U_positionsCMP 


# In[3]:


#definition de la fonction qui permet de marquer tous les ordonancement possible
# N number of possible positions
# n number of qubits per position
def check_Schedule(N,n):
    cpt=0
    c=0
    #contains the input of the gate position comparison
    arr =[] 
    control = []
    #the required number of comparison
    for i in range(1,N):
        cpt+=i   
    schedule=QuantumRegister(n*N)
    ancilla = QuantumRegister(n)
    check=QuantumRegister(cpt)
    validity=QuantumRegister(1)
    circuit=QuantumCircuit(schedule,ancilla,check,validity)
    for i in range(N) :
        nbr_cpt = N-i-1
        for j in range(nbr_cpt):
            for e in range(n):
                arr.append(schedule[n*i+e])
            for p in range(n):
                arr.append(schedule[n*i+(j+1)*n+p])
            for a in range(n):
                arr.append(ancilla[a])
            arr.append(check[c])
            c+=1
            circuit.append(positions_XOR(n),arr)
            arr =[]
    for i in range(cpt):
        control.append(check[i])
    circuit.mct(control,validity[0])
    U_checkSchedule = circuit.to_gate()
    U_checkSchedule.name = "checkSchedule"
    return  U_checkSchedule        


# In[25]:


n=1
N=2**n
cmp=0
for i in range(1,N):
    cmp+=i
q=QuantumRegister(N*n,'position')
ancilla=QuantumRegister(n,'ancilla')
v=QuantumRegister(cmp,'v')
a=QuantumRegister(1)
cout=ClassicalRegister(N*n+1,name='c') 
qc=QuantumCircuit(q,ancilla,v,a,cout)
qc.h(q)
arr =[]
qc.append(check_Schedule(N,n),qc.qubits)

qc.draw('mpl')


# In[11]:


## measure the address qubits  
qc.measure([0,1,4],cout)


# In[13]:


# Use  qasm_simulator
simulator = Aer.get_backend('qasm_simulator')

# Execute the circuit on the qasm simulator
job = execute(qc, simulator, shots=1000)

# Grab results from the job
result = job.result()

# Return counts
counts = result.get_counts(qc)
plot_histogram(counts, figsize=(60, 60),title="inequality detector")


# In[26]:


def memoryEstimationOracle2(cmp,n):
    N=2**n
    requiredQubits=n*N+n+cmp+1
    requiredMemorySpace = 2**requiredQubits *32/(1024*1024*1024)
    return [requiredQubits,requiredMemorySpace]


# In[27]:


n=1
cmp=0
for i in range(1,N):
    cmp+=i
memoryEstimationOracle2(cmp,n)


# %%
