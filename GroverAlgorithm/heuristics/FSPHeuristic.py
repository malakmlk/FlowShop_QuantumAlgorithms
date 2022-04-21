from time import process_time
from typing import List
import numpy as np
import itertools

class FSPHeuristic:
    """
    Heurstics for solving PFSP problem
    """
    def __init__(self,nb_machine : int ,nb_jobs : int ,proc_time : List[List[int]]) -> None:
        """
        proc_time : matrix that contains the processing times of all jobs under different machine
        nb_machine = number of machines
        nb_jobs : number of jobs
        """
        self.nb_machine =nb_machine
        self.nb_jobs = nb_jobs
        self.proc_time = proc_time
    
    def johnson(self):
        """ 
        Johnson rules for solving PFSP with 2 machines 
        return : sequence,cmax,list jobs
        """
        # Optimal sequence using johnson rules 
        m1_seq = [ j 
                   for j in range(self.nb_jobs) 
                   if self.proc_time[0][j] < self.proc_time[1][j]
                 ]
        m1_seq.sort(key = lambda x : self.proc_time[0][x] )
        m2_seq = [ j 
                   for j in range(self.nb_jobs) 
                   if self.proc_time[1][j] <= self.proc_time[0][j]
                 ]
        m2_seq.sort(key = lambda x : self.proc_time[1][x],reverse = True)
        seq = m1_seq + m2_seq
        jobs_m1,jobs_m2 = [] , []
        job = {
            "id" : seq[0],
            "start_time" : 0,
            "end_time" : self.proc_time[0][seq[0]]
        }
        jobs_m1.append(job)
        job = {
            "id" : seq[0],
            "start_time" : self.proc_time[0][seq[0]],
            "end_time" : self.proc_time[0][seq[0]] + self.proc_time[1][seq[0]]
        }
        print(job)
        jobs_m2.append(job)
        for job_id in seq[1::]:
            start_time = jobs_m1[-1]["end_time"]
            end_time = start_time + self.proc_time[0][job_id]
            job = {
                "id" : job_id,
                "start_time" : start_time,
                "end_time" : end_time
            }
            jobs_m1.append(job)
            start_time = max(jobs_m2[-1]["end_time"],end_time)
            end_time = start_time + self.proc_time[1][job_id]
            job = {
                "id" : job_id,
                "start_time" : start_time,
                "end_time" : end_time
            }
            jobs_m2.append(job)
        optim_cmax = end_time
        return seq,optim_cmax,[jobs_m1,jobs_m2]    

    
    def NEH(self):
        """
        NEH heuristic for solving 
        """
        
    
    def calculate_makespan(self, seq):
        a = np.array(self.proc_time)
        a = np.transpose(a)
        # Order the jobs (rows) in order of the sequence
        a = a[seq]

        b = np.zeros(a.shape)
        jobs = a.shape[0]
        macs = a.shape[1]

        b[0, 0] = a[0, 0]
        # Build first row
        for i in range(1, macs):
            b[0, i] = a[0, i] + b[0, i - 1]
        # Build first column
        for i in range(1, jobs):
            b[i, 0] = a[i, 0] + b[i - 1, 0]
        # Build the rest
        for i in range(1, jobs):
            for j in range(1, macs):
                b[i, j] = a[i, j] + (b[i - 1, j] if b[i - 1, j] > b[i, j - 1] else b[i, j - 1])

        return int(b[-1, -1])


    def brute_force_makespan(self):
        a = np.array(self.proc_time)
        n = self.nb_jobs
        m = self.nb_machine
        perm = [[*row] for row in list(itertools.permutations([i for i in range(n)]))]
        opt_cmax = 1000000
        opt_perm = []
        for per in perm :
            cmax = self.calculate_makespan( a,per)
            if cmax <= opt_cmax : 
                opt_cmax = cmax
                opt_perm = per
                print(cmax,per)
        return opt_cmax,opt_perm


