
import multiprocessing 
import time
import multiprocessing as mp
from affinic_gap import linear_gap_algorytm
from affinic_gap_local import linear_gap_local_algorytm
from multidimensional_N_W import multidimesional_N_W_algoritm
from multidimensional_N_W_local import multidimesional_N_W_algoritm_local
from linear_memory_2 import linear_algorytm
def alignment_with_queue(q ,seqs, go, ge ,s,ret_max=50, linear_memory = False,local = False):
        if(linear_memory and len(seqs)== 2 and go==0):
            q.put( linear_algorytm(seqs, s, go, ge, ret_max) )
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
        if(linear_memory and len(seqs)== 2 and go==0):
            return linear_algorytm(seqs, s, go, ge, ret_max)
        if(local == False):
            if(go == 0):
                return multidimesional_N_W_algoritm(seqs, ge,s ,ret_max)
            elif(len(seqs)== 2):
                return linear_gap_algorytm(seqs,  go, ge,s ,ret_max)          
        elif(local == True):
            if(go == 0):
                return multidimesional_N_W_algoritm_local(seqs, ge,s ,ret_max)
            elif(len(seqs)== 2):
                #print("OK")
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
            return -1
        else:
            p.join()
            if(not q.empty()):
                return q.get()
    else:
        return alignment_with_no_time_limit(seqs, go, ge ,s,ret_max , linear_memory ,local)
        

if __name__ == "__main__":
    s1= "AAAA"
    s2= "AA"
    s3= "A"
    tab = [s1, s2 ,s3]
    ge = -0.5
    go = 0
    s =[4,-3]
    max_mathing = 3
    """X = alignment_with_no_time_limit(tab, go, ge , s ,max_mathing ,local = False)
    #print(X)
    s1= "AAAA"
    s2= "AA"
    s3= "A"
    tab = [s1, s2 ,s3]
    ge = -0.5
    go = 0
    s =[4,-3]
    max_mathing = 3
    X = alignment_with_no_time_limit(tab, go, ge , s ,max_mathing ,local = True)
    #print(X)
    tab = [s1, s2 ]
    ge = -0.5
    go = -1
    s =[1,-3]
    max_mathing = 50
    X = alignment_with_no_time_limit(tab, go, ge , s ,max_mathing ,local = False)
    #print(X)
    tab = [s1, s2 ]
    ge = -0.5
    go = -1
    s =[1,-3]
    max_mathing = 50
    X = alignment_with_no_time_limit(tab, go, ge , s ,max_mathing ,local = True)
    #print(X)"""
    tab = [s1, s2 ]
    ge = -0.5
    go = -1
    s =[4,-3]
    max_mathing = 3
    X = alignment(tab, go, ge , s ,max_mathing ,local = False,time_limit=10)
    print(X)
    
    
    
    
    