
import multiprocessing 
import time
import multiprocessing as mp
from affinic_gap import linear_gap_algorytm
from affinic_gap_local import linear_gap_local_algorytm
from multidimensional_N_W import multidimesional_N_W_algoritm
from multidimensional_N_W_local import multidimesional_N_W_algoritm_local
from linear_memory_2 import linear_algorytm
def alignment_with_queue(q ,seqs, go, ge ,s,ret_max=50, linear_memory = False,local = False):
        if(linear_memory and len(seqs)== 2 and go==0 and local==False):
            L = max(((len(seqs[0])+2) * (len(seqs[1])+2) )**(1/2),17 )
            q.put( linear_algorytm(seqs, ge, s, ret_max,L) )
        elif(linear_memory):
            return
        if(local == False):
            if(go == 0):
                q.put( multidimesional_N_W_algoritm(seqs, ge,s ,ret_max) )
            elif(len(seqs)== 2):
                q.put( linear_gap_algorytm(seqs,  go, ge,s ,ret_max) )         
        elif(local == True):
            if(go == 0):
                q.put( multidimesional_N_W_algoritm_local(seqs, ge,s ,ret_max) )
            elif(len(seqs)== 2):
                q.put( linear_gap_local_algorytm(seqs, go, ge,s ,ret_max)  ) 

def alignment_with_no_time_limit(seqs, go, ge ,s,ret_max=50, linear_memory = False,local = False):
        if(linear_memory and len(seqs)== 2 and go==0 and local==False):
            L = max(((len(seqs[0])+2) * (len(seqs[1])+2) )**(1/2),17 )
            return linear_algorytm(seqs, ge, s, ret_max,L)
        elif(linear_memory):
            return
        if(local == False):
            if(go == 0):
                return multidimesional_N_W_algoritm(seqs, ge,s ,ret_max)
            elif(len(seqs)== 2):
                return linear_gap_algorytm(seqs,  go, ge,s ,ret_max)          
        elif(local == True):
            if(go == 0):
                return multidimesional_N_W_algoritm_local(seqs, ge,s ,ret_max)
            elif(len(seqs)== 2):
                return linear_gap_local_algorytm(seqs, go, ge,s ,ret_max) 

def alignment(seqs, go, ge ,s,ret_max=50, linear_memory = False,local = False,time_limit = False):
    if(time_limit):
        ctx = mp.get_context('spawn')
        q = ctx.Queue()
        p = multiprocessing.Process(target=alignment_with_queue, name="Foo", args=(q,seqs, go, ge ,s,ret_max , linear_memory ,local))
        p.start()
        time_left = time_limit
        while(time_left>0):
            time_left = time_left - 0.5
            time.sleep(0.5)
            if(not p.is_alive()):
                break
        if (p.is_alive()):
            p.terminate()
            p.join()
            print("Program is running too long")
            return -1
        else:
            p.join()
            if(not q.empty()):
                return q.get()
    else:
        return alignment_with_no_time_limit(seqs, go, ge ,s,ret_max , linear_memory ,local)
        


    
    
    
    
    