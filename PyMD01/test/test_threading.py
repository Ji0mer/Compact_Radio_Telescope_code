


import os
import threading
from multiprocessing import Process,Queue
import time

def loop():
    print('thread %s is running...' % threading.current_thread().name)
    n = 0
    while n < 5:
        n = n + 1
        print('thread %s >>> %s' % (threading.current_thread().name, n))
        time.sleep(1)
    print('thread %s ended.' % threading.current_thread().name)

print('thread %s is running...' % threading.current_thread().name)
t = threading.Thread(target=loop, name='LoopThread')
t.start()
for i in range(5):
    print(f" {i} Here is the main threading.")
    time.sleep(0.8)
t.join()
print('thread %s ended.' % threading.current_thread().name)



# 子进程要执行的代码
def run_proc(name,q):
    print('\n\nRun child process %s (%s)...' % (name, os.getpid()))
    for i in range(5):
        print(f"Here still child process. >>> {i}.\n")
        time.sleep(0.5)
        q.put([i,i+1])

if __name__=='__main__':
    q = Queue()
    print('Parent process %s.' % os.getpid())
    p = Process(target=run_proc,args=('test',q))
    print('Child process will start.')
    p.start()
    for i in range(5):
        print(f"Here is parent process. >> {i}")
        temp = q.get()
        print(temp[0],temp[1])
        time.sleep(0.4)
    p.join()
    print('Child process end.')



#%%






