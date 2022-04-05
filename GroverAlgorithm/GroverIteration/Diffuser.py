#!/usr/bin/env python
# coding: utf-8

# In[1]:

from qiskit import *
def diffuser(nqubits):
    qc = QuantumCircuit(nqubits)
    # Apply transformation |s> -> |00..0> (H-gates)
    for qubit in range(nqubits):
        qc.h(qubit)
    # Apply transformation |00..0> -> |11..1> (X-gates)
    for qubit in range(nqubits):
        qc.x(qubit)
    # Do multi-controlled-Z gate
    qc.h(nqubits-1)
    qc.mct(list(range(nqubits-1)), nqubits-1)  # multi-controlled-toffoli
    qc.h(nqubits-1)
    # Apply transformation |11..1> -> |00..0>
    for qubit in range(nqubits):
        qc.x(qubit)
    # Apply transformation |00..0> -> |s>
    for qubit in range(nqubits):
        qc.h(qubit)
    # We will return the diffuser as a gate
    U_s = qc.to_gate()
    U_s.name = "U$_s$"
    return U_s

def oracle(q,n,M,x,up) :
    nqubits = N*n+2*q+q*(M-1)*(N)+2*q+3+n+cmp+2
    qRegister = QuantumRegister(nqubits)
    qc = QuantumCircuit(qRegister) 
    qc.append(Oracle1(q,n,M,x,up),[i for i in range(N*n+q*(M-1)*(N)+4*q+2)]+[nqubits-3])
    qc.append(Oracle2(N,n),[i for i in range(N*n)]+[i for i in range(N*n+q*(M-1)*(N)+4*q+2,N*n+q*(M-1)*(N)+4*q+2+n+cmp+1)]+[nqubits-2])
    qc.ccx(qRegister[nqubits-3],qRegister[nqubits-2],qRegister[nqubits-1])
    #reversibility
    qc.append(Oracle1(q,n,M,x,up).inverse(),[i for i in range(N*n+q*(M-1)*(N)+4*q+2)]+[nqubits-3])
    qc.append(Oracle2(N,n).inverse(),[i for i in range(N*n)]+[i for i in range(N*n+q*(M-1)*(N)+4*q+2,N*n+q*(M-1)*(N)+4*q+2+n+cmp+1)]+[nqubits-2])
    U_oracle = qc.to_gate()
    U_oracle.name = "Oracle"
    return  U_oracle 
# In[ ]:




