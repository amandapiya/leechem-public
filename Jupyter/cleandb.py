import os;
import numpy as np;
from countrows import count_rows;
from shutil import copyfile;
from sbemdb import SBEMDB;
from mapping import Mapping;

def clean_db(db):
    src = db.dbfn
    dst = os.path.join('..','data', '170428_pub_modified.sbemdb')
    copyfile(src, dst)
    db = SBEMDB(dst)
        
    # CONSTRAINTS
    # get trees which have roi
    roi_tree_ids = []
    mapper = Mapping()
    for roi in mapper.roi2sbem:
        if roi != '' and mapper.roi2sbem[roi] != '' and mapper.roi2sbem[roi] is not None:
            roi_tree_ids.append(mapper.roi2sbem[roi])
            
    # get trees, which have at least 1 connection
    cons_tree_ids = []
    c = db.db.cursor()
    query = 'select pre.tid, post.tid'\
           + ' from trees as pre' \
           + ' inner join nodes as a on pre.tid==a.tid' \
           + ' inner join syncons as sca on a.nid==sca.nid' \
           + ' inner join synapses as s on sca.sid==s.sid' \
           + ' inner join syncons as scb on scb.sid==s.sid' \
           + ' inner join nodes as b on scb.nid==b.nid' \
           + ' inner join trees as post on b.tid==post.tid' \
           + ' where a.typ==5 and b.typ==6'
    c.execute(query)
    for row in c.fetchall():
        if row[1] != 444: continue # post tree should be 444
        if row[0] == 444: continue # pre tree should not be 444
        if row[0] not in cons_tree_ids: cons_tree_ids.append(row[0])
        if row[1] not in cons_tree_ids: cons_tree_ids.append(row[1])

    tree_ids_keep = np.intersect1d(roi_tree_ids, cons_tree_ids)
    tree_ids_sql = ','.join(map(str, tree_ids_keep))
    
    # TREES TABLE
    count_rows('trees', before=True, db=db)
    c = db.db.cursor()
    c.execute(f'delete from trees where tid NOT IN ({tree_ids_sql})')
    count_rows('trees', before=False, db=db)
    
    # NODES TABLE
    count_rows('nodes', before=True, db=db)
    c = db.db.cursor()
    c.execute(f'delete from nodes where tid NOT IN ({tree_ids_sql})')
    count_rows('nodes', before=False, db=db)
    
    # NODECONS TABLE
    count_rows('nodecons', before=True, db=db)
    c = db.db.cursor()
    c.execute(f'select nid from nodes where tid IN ({tree_ids_sql})')
    data = c.fetchall()
    remaining_nids = []
    for row in data:
        remaining_nids.append(row[0])
    remaining_nids_sql = ','.join(map(str, remaining_nids))  
    c = db.db.cursor()
    c.execute(f'delete from nodecons where nid1 NOT IN ({remaining_nids_sql}) or nid2 NOT IN ({remaining_nids_sql})')
    count_rows('nodecons', before=False, db=db)
    
    # SYNCONS TABLE
    count_rows('syncons', before=True, db=db)
    c = db.db.cursor()
    c.execute(f'delete from syncons where nid NOT IN ({remaining_nids_sql})')
    count_rows('syncons', before=False, db=db)
    
    # SYNAPSES TABLE
    count_rows('synapses', before=True, db=db)
    c = db.db.cursor()
    c.execute(f'select distinct sid from syncons')
    data = c.fetchall()
    remaining_sids = []
    for row in data:
        remaining_sids.append(row[0])
    remaining_sids_sql = ','.join(map(str, remaining_sids))
    c = db.db.cursor()
    # sid 52: 444->444 | sid 486: 444->546
    c.execute(f'delete from synapses where sid NOT IN ({remaining_sids_sql}) or sid=486')
    count_rows('synapses', before=False, db=db)
    
    # COMMIT
    db.db.commit()
    return db

def clean_db_uct(db):
    src = db.dbfn
    dst = os.path.join('..','data', '170428_ganglion10_modified_uct.sbemdb')
    copyfile(src, dst)
    db = SBEMDB(dst)
        
    # CONSTRAINTS
    # get trees which have roi
    uct_tree_ids = []
    mapper = Mapping()
    for uct in mapper.uct2sbem:
        if uct != '' and mapper.uct2sbem[uct] != '' and mapper.uct2sbem[uct] is not None:
            uct_tree_ids.append(mapper.uct2sbem[uct])
            
    # get trees, which have at least 1 connection
    cons_tree_ids = []
    c = db.db.cursor()
    query = 'select pre.tid, post.tid'\
           + ' from trees as pre' \
           + ' inner join nodes as a on pre.tid==a.tid' \
           + ' inner join syncons as sca on a.nid==sca.nid' \
           + ' inner join synapses as s on sca.sid==s.sid' \
           + ' inner join syncons as scb on scb.sid==s.sid' \
           + ' inner join nodes as b on scb.nid==b.nid' \
           + ' inner join trees as post on b.tid==post.tid' \
           + ' where a.typ==5 and b.typ==6'
    c.execute(query)
    for row in c.fetchall():
        if row[1] != 444: continue # post tree should be 444
        if row[0] == 444: continue # pre tree should not be 444
        if row[0] not in cons_tree_ids: cons_tree_ids.append(row[0])
        if row[1] not in cons_tree_ids: cons_tree_ids.append(row[1])

    tree_ids_keep = np.intersect1d(uct_tree_ids, cons_tree_ids)
    tree_ids_sql = ','.join(map(str, tree_ids_keep))
    
    # TREES TABLE
    count_rows('trees', before=True, db=db)
    c = db.db.cursor()
    c.execute(f'delete from trees where tid NOT IN ({tree_ids_sql})')
    count_rows('trees', before=False, db=db)
    
    # NODES TABLE
    count_rows('nodes', before=True, db=db)
    c = db.db.cursor()
    c.execute(f'delete from nodes where tid NOT IN ({tree_ids_sql})')
    count_rows('nodes', before=False, db=db)
    
    # NODECONS TABLE
    count_rows('nodecons', before=True, db=db)
    c = db.db.cursor()
    c.execute(f'select nid from nodes where tid IN ({tree_ids_sql})')
    data = c.fetchall()
    remaining_nids = []
    for row in data:
        remaining_nids.append(row[0])
    remaining_nids_sql = ','.join(map(str, remaining_nids))  
    c = db.db.cursor()
    c.execute(f'delete from nodecons where nid1 NOT IN ({remaining_nids_sql}) or nid2 NOT IN ({remaining_nids_sql})')
    count_rows('nodecons', before=False, db=db)
    
    # SYNCONS TABLE
    count_rows('syncons', before=True, db=db)
    c = db.db.cursor()
    c.execute(f'delete from syncons where nid NOT IN ({remaining_nids_sql})')
    count_rows('syncons', before=False, db=db)
    
    # SYNAPSES TABLE
    count_rows('synapses', before=True, db=db)
    c = db.db.cursor()
    c.execute(f'select distinct sid from syncons')
    data = c.fetchall()
    remaining_sids = []
    for row in data:
        remaining_sids.append(row[0])
    remaining_sids_sql = ','.join(map(str, remaining_sids))
    c = db.db.cursor()
    # sid 52: 444->444 | sid 486: 444->546
    c.execute(f'delete from synapses where sid NOT IN ({remaining_sids_sql}) or sid=486')
    count_rows('synapses', before=False, db=db)
    
    # COMMIT
    db.db.commit()
    return db
