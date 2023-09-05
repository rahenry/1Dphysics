#!/usr/bin/env python3

import logging, os
from args import Q

class nuHandler(logging.NullHandler):
    def __init__(self):
        super().__init__()
        self.things = []

    def filter(self, record):
        return True
    def handle(self, record):
        self.emit(record)
        if Q.log_print == 1:
            print(self.format(record))
    def emit(self, record):
        self.things.append([record, 0])

    def dump(self, f):
        with open(f, 'a+') as w:
            count = 0
            for r in self.things:
                if r[-1] == 0:
                    count += 1
            w.write(f"Dumping {count} records\n")
            for r in self.things:
                if r[-1] == 0:
                    r[-1] = 1
                    # w.write(str(r))
                    w.write(self.format(r[0])+'\n')
