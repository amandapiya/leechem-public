#!/usr/bin/python

from . import mapping
from . import uCT
from . import sbemdb
import numpy as np

class Geometry:
    '''Class to map uCT coordinates to SBEM coordinates and to retrieve
    soma positions.'''
    def __init__(self):
        mp = mapping.Mapping()
        uct = uCT.UCT()
        db = sbemdb.SBEMDB()
        tids = [ x[0] for x in db.fetch('select tid from trees') ]
        soms = {}
        for tid in tids:
            p = db.nodexyz('tid==%i and typ==1' % tid)
            if len(p[0])>0:
                soms[tid] = [p[0][0], p[1][0], p[2][0]]
        puc = uct.somapos
        pem = puc + np.nan
        for k in range(puc.shape[1]):
            if k+1 in mp.uct2sbem:
                tid = mp.uct2sbem[k+1]
                if tid in soms:
                    pem[0,k] = soms[tid][0]
                    pem[1,k] = soms[tid][1]
                    pem[2,k] = soms[tid][2]

        use = ~np.isnan(pem[0,:])
        N = np.sum(use)
        pem1 = np.transpose(np.vstack((pem[:,use], np.ones((1,N)))))
        puc1 = np.transpose(np.vstack((puc[:,use], np.ones((1,N)))))

        self.A_ = np.linalg.lstsq(puc1, pem1, rcond=-1)[0]

        self.B_ = np.linalg.inv(self.A_)

        pem1 = pem[:,use]
        puc1 = np.matmul(puc1, self.A_)
        puc1 = np.transpose(puc1[:,0:3])

        N = puc.shape[1]
        puc1 = np.transpose(np.vstack((puc, np.ones((1,N)))))
        puc1 = np.matmul(puc1, self.A_)
        puc1 = np.transpose(puc1[:, 0:3])

        for k in range(N):
            if k+1 in mp.uct2sbem:
                tid = mp.uct2sbem[k+1]
                if tid not in soms:
                    soms[tid] = puc1[:,k]

        self.soms = soms
        
    def mapUCTtoSBEM(self, xyz_uct):
        '''MAPUCTTOSBEM - Map uCT pixel coordinates to SBEM space
        xyz_sbem = MAPUCTTOSBEM(xyz_uct), where XYZ_UCT are uCT pixel
        coordinates (either a 3-vector or an Nx3 array), returns the
        corresponding coordinates in SBEM space (again, a vector or Nx3 
        array). The result is expressed in micrometers.'''
        xyz_uct = np.array(xyz_uct) # Convert list/tuple to array
        S = xyz_uct.shape
        if len(S)==1:
            xyz_uct = np.reshape(xyz_uct, [1, 3])
        L = xyz_uct.shape[0]
        xyz_uct = np.append(xyz_uct, np.ones((L,1)), 1)
        xyz_sbem = np.matmul(xyz_uct, self.A_)
        xyz_sbem = xyz_sbem[:,:3]
        if len(S)==1:
            xyz_sbem = np.reshape(xyz_sbem, [3])
        return xyz_sbem

    def mapSBEMtoUCT(self, xyz_sbem):
        '''MAPSBEMTOUCT - Map sbem micron coordinates to UCT pixels
        xyz_uct = MAPSBEMTOUCT(xyz_sbem), where XYZ_SBEM are SBEM micron
        coordinates (either a 3-vector or an Nx3 array), returns the
        corresponding coordinates in UCT space (again, a vector or Nx3 
        array). The result is expressed in pixels.'''
        xyz_sbem = np.array(xyz_sbem) # Convert list/tuple to array
        S = xyz_sbem.shape
        if len(S)==1:
            xyz_sbem = np.reshape(xyz_sbem, [1, 3])
        L = xyz_sbem.shape[0]
        xyz_sbem = np.append(xyz_sbem, np.ones((L,1)), 1)
        xyz_uct = np.matmul(xyz_sbem, self.B_)
        xyz_uct = xyz_uct[:,:3]
        if len(S)==1:
            xyz_uct = np.reshape(xyz_uct, [3])
        return xyz_uct

    def somaLocation(self, tid):
        '''SOMALOCATION - Location of soma in SBEM space
        xyz_sbem = SOMALOCATION(tid), where TID is a SBEM Tree ID, returns
        the location of the soma of that tree in SBEM space (in microns).'''
        if tid in self.soms:
            return self.soms[tid]
        else:
            return None

def lr(xyz_sbem):
    '''LR - Anatomical left-to-right coordinate from SBEM space
    p = LR(xyz_sbem), where XYZ_SBEM is either a 3-vector or an Nx3 array
    of points in SBEM space, returns the left-to-right coordinate. (The
    axis runs from negative values on the left to positive values on the
    right.)'''
    S = xyz_sbem.shape
    if len(S)==2:
        return -xyz_sbem[:,1]
    else:
        return -xyz_sbem[1]

def pa(xyz_sbem):
    '''PA - Anatomical posterior-to-anterior coordinate from SBEM space
    p = PA(xyz_sbem), where XYZ_SBEM is either a 3-vector or an Nx3 array
    of points in SBEM space, returns the posterior-to-anterior coordinate.
    (The axis runs from negative values on the posterior end to positive
    values on the anterior end.)'''
    S = xyz_sbem.shape
    if len(S)==2:
        return -xyz_sbem[:,2]
    else:
        return -xyz_sbem[2]

def vd(xyz_sbem):
    '''VD - Anatomical ventral-to-dorsal coordinate from SBEM svdce
    p = VD(xyz_sbem), where XYZ_SBEM is either a 3-vector or an Nx3 array
    of points in SBEM space, returns the ventral-to-dorsal coordinate.
    (The axis runs from negative values on the ventral side to positive
    values on the dorsal side.)'''
    S = xyz_sbem.shape
    if len(S)==2:
        return xyz_sbem[:,0]
    else:
        return xyz_sbem[0]
    
