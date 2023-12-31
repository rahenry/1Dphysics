#!/usr/bin/env python3

import multiprocessing as mp

def foo(q):
    q.put('hello')

if __name__ == '__main__':
    q = mp.Queue()
    p = mp.Process(target=foo, args=(q,))
    p.start()
    print(q.get())
    p.join()
