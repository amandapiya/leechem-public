#!/usr/bin/python3

import numpy as np
import csv
import errno
import os
import re

class Confidence:
    def __init__(self, csvfn=None):
        '''CONFIDENCE - Python access to the SBEM-uCT-VSD confidence file
        conf = CONFIDENCE(csvfn) opens the given confidence file.
        conf = CONFIDENCE() uses the default file in the em170428 directory.
        Result has several fields that each are dicts with tree IDs as keys:
           uctid (numeric)
           vsdid (letters)
           canoid
           sbemconf (0-100) - confidence of tracing
           vsductconf (0-100) - confidence of vsd to uct mapping
           gmapconf (0-100) - confidence of gmapimg results'''
        self.uctid = {}
        self.vsdid = {}
        self.canoid = {}
        self.sbemconf = {}
        self.vsductconf = {}
        self.gmapconf = {}
        self.vsd2tree = {}

        if csvfn is None:
            here = os.path.dirname(__file__)
            csvfn = here + '/../data/confidence.csv'
        lines = []
        with open(csvfn) as f:
            dl = csv.unix_dialect()
            rdr = csv.reader(f, dl)
            for row in rdr:
                lines.append(row)
        hdr = lines.pop(0)
        r = re.compile('^(\d+)\s*(\((\d+)(-([a-z]+))?\))?$')
        for l in lines:
            ids = l[0]
            cano = l[1]
            sbemc = l[4]
            vumapc = l[5]
            gmapc = l[6]
            m = r.match(ids)
            if m:
                try:
                    tid = m.group(1)
                    uctid = m.group(3)
                    vsdid = m.group(5)
                    if tid is not None:
                        tid = int(tid)
                    if uctid is not None:
                        uctid = int(uctid)
                    self.uctid[tid] = uctid
                    self.vsdid[tid] = vsdid
                    self.vsd2tree[vsdid] = tid
                    self.canoid[tid] = cano
                    if sbemc=='':
                        self.sbemconf[tid] = 100
                        print('Caution: no sbemconf for', ids, '- assuming 100%')
                    else:
                        self.sbemconf[tid] = int(sbemc)
                    if vumapc=='-' or vumapc=='':
                        self.vsductconf[tid] = None
                    else:
                        self.vsductconf[tid] = int(vumapc)
                    if gmapc=='-' or gmapc=='' or gmapc=='?':
                        self.gmapconf[tid] = None
                    else:
                        self.gmapconf[tid] = int(gmapc)
                except:
                    print('Something wrong at', ids)
                    raise
            else:
                pass # print('no match', ids)

    def tracingConfidenceForTree(self, tid):
        '''TRACINGCONFIDENCEFORTREE - Return the (subjective) confidence
        in tracing for the given tree ID as a percentage value. Only
        trees that are presynaptic to DE-3(R) are rated; others return None.'''
        if tid in self.gmapconf:
            return self.gmapconf[tid]
        else:
            return 0

    def roiMappingConfidenceForTree(self, tid):
        '''ROIMAPPINGCONFIDENCEFORTREE - Return the (subjective) confidence
        in mapping the given tree to an ROI in the VSD data as a percentage 
        value. Only trees that are presynaptic to DE-3(R) are rated; others
        return None.'''
        if tid in self.vsductconf:
            return self.vsductconf[tid]
        else:
            return None

    def canonicalMappingConfidenceForROI(self, roiid):
        '''CANONICALMAPPINGCONFIDENCEFORROI - Return the (subjective) 
        confidence in mapping the given ROI to a canonical cell as a
        percentage value. Only ROIs that are associated with an SBEM tree ID
        presynaptic to DE-3(R) are rated. Others return None.'''
        if roiid in self.vsd2tree:
            tid = self.vsd2tree[roiid]
            if tid in self.gmapconf:
                return self.gmapconf[tid]
        return None
            
if __name__=='__main__':
    print('Confidence checker')
    import mapping
    c = Confidence()
    m = mapping.Mapping()
    for k in c.uctid.keys():
        if k not in m.sbem2uct:
            print('confidence has key that is missing from mapping:', k)
        else:
            if c.uctid[k] != m.sbem2uct[k]:
                print('Mismatch for tree', k, ':', c.uctid[k], 'vs', m.sbem2uct[k], 'vsd', c.vsdid[k], 'cano', c.canoid[k])
