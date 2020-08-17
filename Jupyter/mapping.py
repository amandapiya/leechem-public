#!/usr/bin/python3

import numpy as np
import csv
import errno
import os

class Mapping:
    def __init__(self, csvfn=None):
        '''MAPPING - Python access to the SBEM-uCT-VSD mapping
        map = MAPPING(csvfn) opens the given mapping file.
        map = MAPPING() uses the default file in the em170428 directory.
        Result has fields roi2can, roi2uct, roi2sbem that map ROI numbers
        to other IDs; sbem2can, sbem2uct, sbem2roi, sbem2roiid, sbem2tname
        that map SBEM ID to other IDs, and uct2can, uct2roi, uct2sbem that
        map uCT ID numbers to other IDs.'''

        self.roi2can = {}
        self.roi2uct = {}
        self.roi2sbem = {}
        self.sbem2can = {}
        self.sbem2uct = {}
        self.sbem2roi = {}
        self.sbem2roiid = {}
        self.sbem2tname = {}
        self.uct2can = {}
        self.uct2roi = {}
        self.uct2sbem = {}
        self.can2sbem = {}
        self.can2roi = {}
        self.can2uct = {}
        self.roiid2roi = {}
        self.roi2roiid = {}
        
        if csvfn is None:
            here = os.path.dirname(__file__)
            csvfn = here + '/../data/mapping.csv'
        lines = []
        with open(csvfn) as f:
            dl = csv.unix_dialect()
            rdr = csv.reader(f, dl)
            for row in rdr:
                lines.append(row)
        hdr = lines.pop(0)
        for l in lines:
            roi = self.convert_to_number(l[0])
            roiid = l[1]
            can = l[2]
            uct = self.convert_to_number(l[3])
            sbem = self.convert_to_number(l[5])
            tname = l[6]
            if can=='':
                can = None
            if roiid=='':
                roiid = None
                roi = None
            if roi is not None:
                self.roi2can[roi] = can
                self.roi2uct[roi] = uct
                self.roi2sbem[roi] = sbem
                self.roiid2roi[roiid] = roi
                self.roi2roiid[roi] = roiid
            if sbem is not None:
                self.sbem2can[sbem] = can
                self.sbem2uct[sbem] = uct
                self.sbem2roi[sbem] = roi
                self.sbem2roiid[sbem] = roiid
                self.sbem2tname[sbem] = tname
            if uct is not None:
                self.uct2can[uct] = can
                self.uct2sbem[uct] = sbem
                self.uct2roi[uct] = roi
            if can is not None:
                if can in self.can2roi:
                    print(f'Duplicate can: {can}')
                self.can2sbem[can] = sbem
                self.can2uct[can] = uct
                self.can2roi[can] = roi

    def mapTreeIDToUCTID(self, tid):
        '''MAPUCTIDTOUCTID - Return UCT ID associated with 
        given tree ID, or None if there isn't one.'''
        if tid in self.sbem2uct:
            return self.sbem2uct[tid]
        else:
            return None

    def mapTreeIDToROIID(self, tid):
        '''MAPTREEIDTOROIID - Return ROI ID associated with 
        given tree ID, or None if there isn't one.'''
        if tid in self.sbem2roiid:
            return self.sbem2roiid[tid]
        else:
            return None

    def mapTreeIDToCanonicalName(self, tid):
        '''MAPTREEIDTOCANONICALNAME - Return canonical name associated with 
        given tree ID, or None if there isn't one.'''
        if tid in self.sbem2can:
            return self.sbem2can[tid]
        else:
            return None

    def mapROIIDToTreeID(self, roiid):
        '''MAPROIIDTOTREEID - Return tree ID associated with 
        given ROI, or None if there isn't one.'''
        if roiid in self.roiid2roi:
            roi = self.roiid2roi[roiid]
            return self.roi2sbem[roi]
        else:
            return None

    def mapROIIDToUCTID(self, roiid):
        '''MAPROIIDTOUCTID - Return UCT ID associated with 
        given ROI, or None if there isn't one.'''
        if roiid in self.roiid2roi:
            roi = self.roiid2roi[roiid]
            return self.roi2uct[roi]
        else:
            return None

    def mapROIIDToCanonicalName(self, roiid):
        '''MAPROIIDTOCANONICALNAME - Return canonical name associated with 
        given ROI, or None if there isn't one.'''
        if roiid in self.roiid2roi:
            roi = self.roiid2roi[roiid]
            return self.roi2can[roi]
        else:
            return None

    def mapROIIDToNumber(self, roiid):
        '''MAPROIIDTONUMBER - Map a ROI ID to a ROI number.
        ROI IDs are a one- or two-letter string, like "a" or "fr". These
        get mapped to numbers counting from 1.'''
        if roiid in self.roiid2roi:
            return self.roiid2roi[roiid]
        else:
            return None

    def mapCanonicalNameToROIID(self, can):
        '''MAPCANONICALNAMETOROIID - Return ROI ID instantiating given
        canonical name, or None if there isn't one.'''
        if can in self.can2roi:
            roi = self.can2roi[can]
            if roi in self.roi2roiid:
                return self.roi2roiid[roi]
        return None

    def mapCanonicalNameToTreeID(self, can):
        '''MAPCANONICALNAMETOTREEID - Return tree ID associated with
        given canonical name, or None if there isn't one.'''
        if can in self.can2sbem:
            return self.can2sbem[can]
        else:
            return None

    def canonicalNames(self):
        '''CANONICALNAMES - Return a list of all canonical names for which 
        mapping info exists.'''
        return list(self.can2roi.keys())

    def roiIDs(self):
        'ROIIDS - Return a list of all ROI IDs for which mapping info exists.'
        return list(self.roiid2roi.keys())

    def treeIDs(self):
        '''TREEIDS - Return a list of all tree IDs for which mapping info
        exists.'''
        return list(self.sbem2uct.keys())
        
                
    def convert_to_number(self, s):
        try:
            res = int(s)
        except ValueError:
            res = None
        return res
