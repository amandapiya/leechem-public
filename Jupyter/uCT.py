#!/usr/bin/python3

import numpy as np
import h5py
import errno
import os

class UCT:
    def __init__(self, h5fn=None):
        '''UCT - Python access to the uCT soma position file
        uct = UCT(h5fn) opens the given "uCT-somata.h5" file.
        uct = UCT() opens the default file in the em170428 directory.
        '''
        if h5fn is None:
            here = os.path.dirname(__file__)
            h5fn = here + '/../data/uCT-somata.h5'
            with h5py.File(h5fn, 'r') as h5:
                self.somapos = h5['somapos']['value'][:] * 2
                self.exitpoint = h5['dpp']['value'][:]
                self.lut = h5['dc']['value'][:]
                self.tails = {}
                tt = h5['tails']['value']
                for k in tt.keys():
                    if k.startswith('_'):
                        k1= int(k[1:])
                        t = tt[k]
                        t = t['value']
                        if len(t)==3:
                            self.tails[k1] = t[:]
                        else:
                            self.tails[k] = [[],[],[]] 
    def somaPosition(self, uctid):
        '''SOMAPOSITION - Position of a soma
        x,y,z = SOMAPOSITION(uctid) returns the position of a soma.
        UCTID must be a UCT ID (for instance, DE3_R is 1; DI1_R is 66).
        X, Y, Z are in UCT pixel coordinates.'''
        return self.somapos[:,uctid-1]
    def exitPointPosition(self, uctid):
        '''EXITPOINTPOSITION - Position of a exit point
        x,y,z = EXITPOINTPOSITION(uctid) returns the position of the point
        where the principal neurite of a cell exits the neuropil.
        UCTID must be a UCT ID (for instance, DE3_R is 1; DI1_R is 66).
        X, Y, Z are in UCT pixel coordinates.'''
        return self.exitpoint[:,uctid-1]
    def tailPosition(self, uctid):
        '''TAILPOSITION - Position of a tail
        x,y,z = TAILPOSITION(uctid) returns the position of the “tail”
        of a cell, i.e., the approximate path of the principal neurite
        between the soma and the “exit point.”
        UCTID must be a UCT ID (for instance, DE3_R is 1; DI1_R is 66).
        X, Y, Z are vectors of UCT pixel coordinates.'''
        return self.tails[uctid]
