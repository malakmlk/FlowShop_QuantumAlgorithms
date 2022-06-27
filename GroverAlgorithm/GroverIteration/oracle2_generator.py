#!/usr/bin/env python
# coding: utf-8

# In[1]:


#initialization
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
# importing Qiskit
from qiskit import *
from qiskit.circuit.library import OR
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
def Oracle2(N,n):
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
    U_oracle2 = circuit.to_gate()
    U_oracle2.name = "Oracle2"
    return  U_oracle2       




def memoryEstimationOracle2(n):
    N=2**n
    for i in range(1,N): 
        cmp+=i
    requiredQubits=n*N+n+cmp+1
    requiredMemorySpace = 2**requiredQubits *32/(1024*1024*1024)
    return [requiredQubits,requiredMemorySpace]


# In[27]:


# %%
