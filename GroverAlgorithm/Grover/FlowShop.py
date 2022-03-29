from GroverAlgorithm import Oracle1Generator,Oracle2Generator
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

qc.measure([0,1,4],cout)


# Use  qasm_simulator
simulator = Aer.get_backend('qasm_simulator')

# Execute the circuit on the qasm simulator
job = execute(qc, simulator, shots=1000)

# Grab results from the job
result = job.result()

# Return counts
counts = result.get_counts(qc)
plot_histogram(counts, figsize=(60, 60),title="inequality detector",filename="Result.png")
