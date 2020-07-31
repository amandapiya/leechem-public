from distance import distance;
import warnings;
import networkx as nx;

def find_path(db, nid1, nid2):
    '''
    Finds DISTANCE and PATH between 2 nodes: NID1 and NID2.
    If the nodes belong to different trees, displays 'Nodes are on different trees'
    If one/both nodes are pre/post-synaptic add distance from the node(s) to the middle of the synapse(s) they belong to
        to the total distance
    '''
    if nid1 == nid2:
        raise ValueError('Node ids have to be different')
        
    nodes = db.nodeDetails(f'nid={nid1} or nid={nid2}')
    if nid1 not in nodes or nid2 not in nodes:
        raise ValueError('One/both nodes are not in the database')
        
    tid1 = nodes[nid1].tid
    tid2 = nodes[nid2].tid
    if tid1 != tid2:
        raise ValueError('Nodes are on different trees')
        
    typ1 = nodes[nid1].typ
    typ2 = nodes[nid2].typ
        
    x, y, z, nids = db.nodexyz(f'tid={tid1}')
    # nodes' coordinates map
    node_coords = {}
    for i, nid in enumerate(nids):
        node_coords[nid] = [x[i], y[i], z[i]]
        
    # add distance to the center of synapse
    add_to_length = 0
    if typ1 in [5, 6]: # is synapse (5 - pre, 6-post synaptic)
        (xx, yy, zz, pretid, posttid) = db.synapses(f'a.nid={nid1} or b.nid={nid1}')
        if len(xx) == 0: warnings.warn(f'Synapse node {nid1} is connected to a tree, which is not on the list of available trees.')
        else: add_to_length += distance(xx[0], yy[0], zz[0], *(node_coords[nid1]))
    if typ2 in [5, 6]: # is synapse (5 - pre, 6-post synaptic)
        (xx, yy, zz, pretid, posttid) = db.synapses(f'a.nid={nid2} or b.nid={nid2}')
        if len(xx) == 0: warnings.warn(f'Synapse node {nid2} is connected to a tree, which is not on the list of available trees.')
        else: add_to_length += distance(xx[0], yy[0], zz[0], *(node_coords[nid2]))
    
    c = db.db.cursor()
    c.execute(f'select c.nid1, c.nid2' +
              f' from nodes inner join nodecons as c on nodes.nid==c.nid1 inner join nodes as b on c.nid2==b.nid' +
              f' where nodes.tid={tid1} and b.tid={tid1}')
    # create graph
    G = nx.DiGraph()
    for row in c.fetchall():
        G_nid1, G_nid2 = row[0], row[1]
        G.add_edge(G_nid1, G_nid2, distance=distance(*(node_coords[G_nid1] + node_coords[G_nid2])))
        
    length = nx.dijkstra_path_length(G, nid1, nid2, 'distance')
    length += add_to_length
    path = nx.dijkstra_path(G, nid1, nid2, 'distance')
    return length, path
    
def find_path_ext(db, pre_tid=None, post_tid=None):
    '''
    Takes as arguments presynaptic tree id and postsynaptic tree id.
    Returns maps in the following forms:
        1) {(pre_tid, post_tid): {synapse_id: distance}} # distance is from synapse to post_tree's soma
        2) {(pre_tid, post_tid): {synapse_id: [path]}}   # path are nodes
    '''
    where1 = ''
    where2 = ''
    where3 = ''
    where4 = 'typ=1'
    if pre_tid is not None:
        where1 += f'pre.tid={pre_tid}'
        where2 += f'tid={pre_tid}'
        where3 += f'where nodes.tid={pre_tid} or b.tid={pre_tid}'
    if post_tid is not None:
        if where1 != '':
            where1 += f' and post.tid={post_tid}'
        else:
            where1 += f'post.tid={post_tid}'

        if where2 != '':
            where2 += f' or tid={post_tid}'
        else:
            where2 += f'tid={post_tid}'

        if where3 != '':
            where3 += f' or nodes.tid={post_tid} or b.tid={post_tid}'
        else:
            where3 += f'where nodes.tid={post_tid} or b.tid={post_tid}'

        if where4 != '':
            where4 += f' and tid={post_tid}'

    x, y, z, nids = db.nodexyz(where2)
    # nodes' coordinates map
    node_coords = {}
    for i, nid in enumerate(nids):
        node_coords[nid] = [x[i], y[i], z[i]]


    select_subject = 'c.nid1, c.nid2'
    from_tables = 'nodes inner join nodecons as c on nodes.nid==c.nid1 inner join nodes as b on c.nid2==b.nid'
    c = db.db.cursor()
    c.execute(f'select {select_subject}'
              f' from {from_tables} {where3}')
    # create graph
    G = nx.DiGraph()
    for row in c.fetchall():
        G_nid1, G_nid2 = row[0], row[1]                
        G.add_edge(G_nid1, G_nid2, distance=distance(*(node_coords[G_nid1] + node_coords[G_nid2])))

    # get targets(somas) nids
    treeid2nid = {}
    data = db.nodeDetails(where4)
    for nid in data:
        tid = data[nid].tid
        treeid2nid[tid] = nid

    # get synapses 
    xx, yy, zz, pretid, posttid, synid, prenid, postnid = db.synapses(where1, extended=True)

    # calculations
    results_length = {} # {(source_tree_id, target_tree_id): {synapse_id: distance}}
    results_path = {} # {(source_tree_id, target_tree_id): {synapse_id: path}}
    for i in range(len(pretid)):
        pre_tid = pretid[i]
        post_tid = posttid[i]
        post_nid = postnid[i]
        syn_id = synid[i]
        if (pre_tid, post_tid) not in results_length: 
            results_length[(pre_tid, post_tid)] = {}
            results_path[(pre_tid, post_tid)] = {}
        if post_tid in treeid2nid: # if has soma
            length = nx.dijkstra_path_length(G, post_nid, treeid2nid[post_tid], 'distance')
            length += distance(xx[i], yy[i], zz[i], *node_coords[post_nid])
            path = nx.dijkstra_path(G, post_nid, treeid2nid[post_tid], 'distance')
            results_length[(pre_tid, post_tid)][syn_id] = length
            results_path[(pre_tid, post_tid)][syn_id] = path
        else:
            results_length[(pre_tid, post_tid)][syn_id] = 'Target tree has no soma'
            results_path[(pre_tid, post_tid)][syn_id] = 'Target tree has no soma'

    return results_length, results_path