#!/usr/bin/python3

import numpy as np
import sqlite3
import errno
import os

class Node:
    def __init__(self, nid, tid, typ, x, y, z, cdate, uid):
        '''NODE - A tree node from the NODES table.
        Contains nid, tid, typ, x, y, z, cdate, uid fields from the table.'''
        self.nid = nid
        self.tid = tid
        self.typ = typ
        self.x = x
        self.y = y
        self.z = z
        self.cdate = cdate
        self.uid = uid

class SBEMDB:
    def __init__(self, dbfn=None):
        '''SBEMDB - Pythonic access to SBEMDB
        db = SBEMDB(dbfn) opens the given database file. db = SBEMDB() opens
        the database file in the em170428 directory.'''
        if dbfn is None:
            here = os.path.dirname(__file__)
            dbfn = here + '/../data/170428_pub.sbemdb'
        if not os.path.isfile(dbfn):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                    dbfn)
        self.db = sqlite3.connect(dbfn)
        self.dbfn = dbfn

    def fetch(self, query):
        c = self.db.cursor()
        c.execute(query)
        return c.fetchall()
        
    def nodetypes(self):
        '''NODETYPES - Get nodetype enum
        NODETYPES() returns a mapping of node type numbers to descriptions.'''
        res = {}
        for row in self.fetch('select typ, descr from nodetypes'):
            res[row[0]] = row[1]
        return res
    
    def nodexyz(self, where=''):
        '''NODEXYZ - Get a map of node location
        (x, y, z, nid) = NODEXYZ(whereclause) performs a SELECT on the NODES 
        with the given where clause. For instance,
        nnn = db.nodes('tid==444 and typ==6') finds
        all the postsynaptic terminals on tree #444.
        (TYP is useful for selecting node types more broadly. Values are:
          1: Soma
          2: Exit point from neuropil
          3: Regular tree node
          5: Presynaptic site
          6: Postsynaptic site.)
        Results X, Y, Z are expressed in um. 
        Positive X is to anatomical left.
        Positive Y is to posterior.
        Positive Z is to ventral side.
        ID is node ID'''
        if where != '':
            where = f'where {where}'
        rows = self.fetch(f'select nid, x, y, z from nodes {where}')
        N = len(rows)
        x = np.zeros(N)
        y = np.zeros(N)
        z = np.zeros(N)
        nid = np.zeros(N, dtype=int)
        for n in range(N):
            x[n] = rows[n][1]
            y[n] = rows[n][2]
            z[n] = rows[n][3]
            nid[n] = rows[n][0]
        x = self.pixtoum(x)
        y = self.pixtoum(y)
        z = self.slicetoum(z)
        return (x, y, z, nid)

    def onenodexyz(self, nid):
        '''ONENODEXYZ - Get location of a single node
        (x, y, z) = ONENODEXYZ(nid) returns the location of the given node.
        The result is in microns'''
        x,y,z,nid = self.nodexyz(f'nid={nid}')
        if len(x):
            return x[0], y[0], z[0]
        else:
            raise ValueError(f'No such node: {nid}')

    def somaxyz(self, tid):
        '''SOMAXYZ - Get soma location and node ID
        (x, y, z, nid) = SOMAXYZ(tid) returns the soma location and node ID
        for the given tree.'''
        x,y,z,nid = self.nodexyz(f'tid={tid} and typ=1')
        if len(x):
            return x[0], y[0], z[0], nid[0]
        else:
            return None, None, None, None

    def nodeDetails(self, where=''):
        '''NODEDETAILS - Get a map of node details
        nnn = NODEDETAILS(whereclause) performs a SELECT on the NODES with the
        given WHERE clause. For instance, nnn = db.nodes('tid==444 and typ==6')
        finds all the postsynaptic terminals on tree #444.
        Result is a map from nid to class Node.
        If WHERE is not given, returns all nodes.'''
        if where!='':
            where = f'where {where}'
        res = {}
        query = f'select nid, tid, typ, x, y, z, cdate, uid from nodes {where}'
        for row in self.fetch(query):
            res[row[0]] = Node(row[0], row[1], row[2],
                               row[3], row[4], row[5],
                               row[6], row[7])
        return res
        
    def pixtoum(self, x, a=0):
        return x * .0055 * 2.**a

    def umtopix(self, x, a=0):
        return (x / 2.**a) / .0055

    def slicetoum(self, z):
        return z * .050

    def umtoslice(self, z):
        return z / .050

    def presyntrees(self, tid=444, where=''):
        '''PRESYNTREES - Get list of trees presynaptic to a given other tree
        pst = PRESYNTREES(tid) returns information on all the trees that are
        presynaptic to the given tree. If TID is not given, it defaults to
        444, which is the tree ID for cell DE-3.
        The result is a dict from TID (tree IDs), to tuples (SYNCNT, TNAME)
        where SYNCNT is a count of synapses from each tree to the target
        and TNAME is the name of the tree.
        PRESYNTREES(tid, clause) where CLAUSE is a SQL where clause in terms
        of, e.g., trees.tname adds additional constraints. '''
        clause = f'b.tid={tid}'
        if where != '':
            clause += f' and ({where})'

        query = f'''select a.tid, count(1) as cnt, tname
        from nodes as a
        inner join syncons as sca on a.nid==sca.nid
        inner join synapses as s on sca.sid==s.sid 
        inner join syncons as scb on scb.sid==s.sid
        inner join nodes as b on scb.nid==b.nid
        inner join trees on a.tid==trees.tid
        where ({clause}) and a.typ==5 and b.typ==6
        group by a.tid'''
        res = {}
        for row in self.fetch(query):
            res[row[0]] = (row[1], row[2])
        return res

    def simplesegments(self, where):
        '''SIMPLESEGMENTS - Get coordinates of a tree
        segs = SIMPLESEGMENTS(tid) returns all segments of a tree. TID must
        be a numeric tree ID.
        segs = SIMPLESEGMENTS(clause) returns a subset of segments.
        CLAUSE should be a WHERE clause that specifies a tree, for
        instance, 'nodes.tid==444'.
        SEGS will be a Nx6 matrix of (x1,y1,z1,x2,y2,z2) coordinates.
        Coordinates are returned in microns.'''
        if not (type(where)==str):
            clause = f'a.tid={where}'
        else:
            clause = f'({where})'
        query = f'''select a.x as ax, a.y as ay, a.z as az,
            b.x as bx,b.y as by,b.z as bz
            from nodes as a
            inner join nodecons as nc on a.nid==nc.nid1
            inner join nodes as b on nc.nid2==b.nid
            where a.nid<b.nid and ({clause})'''
        rows = self.fetch(query)
        K = len(rows)
        res = np.zeros((K, 6))
        for k in range(K):
            res[k,:] = np.array(rows[k])
        for c in [0, 1, 3, 4]:
            res[:,c] = self.pixtoum(res[:,c])
        for c in [2, 5]:
            res[:,c] = self.slicetoum(res[:,c])

        return res
    
    def segments(self, where):
        '''SEGMENTS - Get coordinates of a tree
        (xx, yy, zz) = SEGMENTS(tid) returns all segments of a tree. TID must
        be a numeric tree ID.
        (xx, yy, zz) = SEGMENTS(clause) returns a subset of segments.
        CLAUSE should be a WHERE clause that specifies a tree, for
        instance, 'nodes.tid==444'.
        Results XX, YY, ZZ can be used directly for plotting.
        Coordinates are returned in microns;
        breaks in the data are marked by NANs.'''
        if not (type(where)==str):
            clause = f'nodes.tid={where}'
        else:
            clause = f'({where})'
		
        c = self.db.cursor()
        wht = '''nodes.nid as aid,
           nodes.x as ax, nodes.y as ay, nodes.z as az,
           b.nid as bid,b.x as bx,b.y as by,b.z as bz'''
        frm = ''' from nodes inner join nodecons as c on nodes.nid==c.nid1
           inner join nodes as b on c.nid2==b.nid'''
        query = f'''select {wht} {frm} 
            where aid<bid and ({clause}) order by aid'''
        xx = []
        yy = []
        zz = []
        lastid = None
        for row in self.fetch(query):
            aid = row[0]
            ax = row[1]
            ay = row[2]
            az = row[3]
            bid = row[4]
            bx = row[5]
            by = row[6]
            bz = row[7]
            if aid != lastid:
                if lastid is not None:
                    xx.append(np.nan)
                    yy.append(np.nan)
                    zz.append(np.nan)
                xx.append(ax)
                yy.append(ay)
                zz.append(az)
            xx.append(bx)
            yy.append(by)
            zz.append(bz)
            lastid = bid
        xx = self.pixtoum(np.array(xx))
        yy = self.pixtoum(np.array(yy))
        zz = self.slicetoum(np.array(zz))
        return (xx, yy, zz)      

    def synapses(self, where='', extended=False):
        '''SYNAPSES - Find position of synapses
        (xx, yy, zz, pretid, posttid) = SYNAPSES(where) returns the 
        coordinates and pre- and postsynaptic tree IDs of synapses.
        WHERE may specify something about the presynaptic tree or the 
        postsynaptic tree. E.g., 'post.tid==444'.
        If WHERE is not given, then all synapses connections will be retrieved.
        Results XX, YY, ZZ, expressed in microns.
        The result is the average of the position of the presynaptic node and
        of the postsynaptic node.
        Only finds synapses that have at least one presynaptic terminal and
        at least one postsynaptic terminal.
        EXTENDED tells if synapse ids, presynaptic and postsynaptic node ids
        should be added to output. In that case, the return value is
        (xx, yy, zz, pretid, posttid, sid, prenid, postnid).'''
        clause = f'a.typ==5 and b.typ==6'
        if where != '': 
            clause += f' and ({where})'

        query = f'''select a.x as ax, a.y as ay, a.z as az,
                b.x as bx, b.y as by, b.z as bz,
                pre.tid as atid, post.tid as btid,
                s.sid as sid, a.nid as prenid, b.nid as postnid
                from trees as pre
                inner join nodes as a on pre.tid==a.tid
                inner join syncons as sca on a.nid==sca.nid
                inner join synapses as s on sca.sid==s.sid
                inner join syncons as scb on scb.sid==s.sid
                inner join nodes as b on scb.nid==b.nid
                inner join trees as post on b.tid==post.tid
                where {clause}'''
        
        rows = self.fetch(query)
        N = len(rows)
        xx = np.zeros(N)
        yy = np.zeros(N)
        zz = np.zeros(N)
        pretid = np.zeros(N, dtype=int)
        posttid = np.zeros(N, dtype=int)
        synid = np.zeros(N, dtype=int)
        prenid = np.zeros(N, dtype=int)
        postnid = np.zeros(N, dtype=int)
        for n in range(N):
            row = rows[n]
            xx[n] = (row[0]+row[3])/2.0
            yy[n] = (row[1]+row[4])/2.0
            zz[n] = (row[2]+row[5])/2.0
            pretid[n] = row[6]
            posttid[n] = row[7]
            synid[n] = row[8]
            prenid[n] = row[9]
            postnid[n] = row[10]
        xx = self.pixtoum(xx)
        yy = self.pixtoum(yy)
        zz = self.slicetoum(zz)
        output = (xx, yy, zz, pretid, posttid)
        if extended:
            output += (synid, prenid, postnid)
        return output
    
    def distanceAlongTree(self, nid):
        '''DISTANCEALONGTREE - Distance along tree between nodes
        dd = DISTANCEALONGTREE(nid) returns the distance along the tree from
        the given node as a dict mapping nodes to distance in microns.
        See example of use in demo.py.'''
            
        res = self.fetch(f'select tid from nodes where nid={nid}')
        if len(res) != 1:
            raise KeyError(f'Node {nid} not found')
        tid = res[0][0]
        res = self.fetch(f'''select
           n1.x, n1.y, n1.z, n1.nid, n2.x, n2.y, n2.z, n2.nid
           from nodes as n1
           inner join nodecons as s on n1.nid=s.nid1
           inner join nodes as n2 on n2.nid=s.nid2
           where n1.nid<n2.nid and n1.tid={tid}''')
        neighbors = {}
        # NEIGHBORS will be a map from nid to a dict of (nid: distance) pairs
        facxy = self.pixtoum(1)**2
        facz = self.slicetoum(1)**2
        for x1, y1, z1, n1, x2, y2, z2, n2 in res:
            dst = np.sqrt(facxy*(x1-x2)**2 + facxy*(y1-y2)**2 + facz*(z1-z2)**2)
            if n1 not in neighbors:
                neighbors[n1] = {}
            neighbors[n1][n2] = dst
            if n2 not in neighbors:
                neighbors[n2] = {}
            neighbors[n2][n1] = dst

        distances = { nid: 0 }
        frontier = [ nid ]
        while len(frontier):
            newfrontier = []
            for k in frontier:
                dst = distances[k]
                for n,d in neighbors[k].items():
                    if n not in distances:
                        newfrontier.append(n)
                        distances[n] = dst + d
            frontier = newfrontier
        return distances

    def nodeToEdgeDistance(self, nid, nid1, nid2):
        '''NODETOEDGEDISTANCE - Distance between node and an edge 
        NODETOEDGEDISTANCE(nid, nid1, nid2) returns the distance between 
        the node defined by NID and the EDGE defined by the two nodes 
        NID1 and NID2. The result is expressed in microns.'''
        p = self.onenodexyz(nid)
        p1 = self.onenodexyz(nid1)
        p2 = self.onenodexyz(nid2)
        def dif(pa, pb):
            return [pa[k] - pb[k] for k in range(3)]
        def inner(a, b):
            return sum([a[k]*b[k] for k in range(3)])
        def addmul(p, dp, t):
            return [p[k] + t*dp[k] for k in range(3)]
        def length(dp):
            return np.sqrt(inner(dp,dp))
        # Let P, P1, P2 be the points (x,y,z), (x1,y1,z1) and (x2,y2,z2)
        # The line segment is defined by
        # P' = P1 + T (P2-P1)
        # We are seeking the value of T that minimizes D := ||P' - P||
        # constrained by 0 <= T <= 1.
        # Find that by demanding that d(D²)/dT = 0:
        #    D² = sum_i [ P1_i + T (P2_i - P1_i) - P_i ]²,
        # where i iterates over x, y, z. Thus:
        #   0 = sum_i [ P1_i + T (P2_i - P1_i) - P_i ] [P2_i - P1_i].
        # Separting terms with and without T is easy.
        #   0 = sum_i [P1_i - P_i][P2_i - P1_i] + T sum_i[P2_i - P1_i]^2
        dif1 = dif(p1, p)
        dif2 = dif(p2, p1)
        t = -inner(dif1, dif2) / inner(dif2, dif2)
        if t<0:
            t = 0
        elif t>1:
            t = 1
        pp = addmul(p1, dif2, t)
        return length(dif(p, pp))
    
