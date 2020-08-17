import math
import numpy as np
import itertools
import networkx


class Synapse:
    def __init__(self, sid, x, y, z, pre_nid, post_nid):
        self.sid = sid

        self.x = x
        self.y = y
        self.z = z

        self.pre_nid = pre_nid
        self.post_nid = post_nid

    def __repr__(self):
        return f"<Synapse ID: {self.sid}>"

class MatrixSynapseDistance:
    """ Distances between synapses returned as matrix
    The 1st row and 1st column of matrix contain synapse id
    """
    def __init__(self, matrix):
        self._matrix = matrix
    
    def distance(self, a, b):
        """Calculates distances based on the given distance matrix"""
        syn_1_ind = np.where(self._matrix[0] == a.sid)[0][0]
        syn_2_ind = np.where(self._matrix[0] == b.sid)[0][0]
        return self._matrix[syn_1_ind, syn_2_ind]

class EuclideanSynapseDistance:
    """Euclidean distance calculator for Synapse objects.

    Instantiation takes no variables.
    """

    @staticmethod
    def distance(a, b):
        """Calculates the euclidean distance between 2 Synapse object positions."""

        return math.sqrt(
            (a.x - b.x) ** 2
            + (a.y - b.y) ** 2
            + (a.z - b.z) ** 2
        )


class GraphSynapseDistance:
    """Graph distance calculator for Synapse objects.

    Instantiation requires networkx graph, with 'distance' field for each edge.
    """

    def __init__(self, graph):
        self._graph = graph

    def distance(self, a, b):
        """Calculates the graph distance between the post nodes of 2 Synapse objects."""

        return networkx.dijkstra_path_length(self._graph, a.post_nid, b.post_nid, 'distance')


class ConstraintChaining:
    """Constrains the hierarchical clustering process by minimal nearest neighbor
    distance within the cluster.
    """

    def __init__(self, distance_calculator, distance_limit):
        self._distance_calculator = distance_calculator
        self._distance_limit = distance_limit

    def check(self, cluster):
        for i, node_i in enumerate(cluster):
            for j, node_j in enumerate(cluster):
                if i == j:
                    continue
                if self._distance_calculator.distance(node_i, node_j) < self._distance_limit:
                    break
            else:
                return False
        return True


class ConstraintDiameter:
    """Constrains the cluster by its diameter (largest distance between any 2 nodes)."""

    def __init__(self, distance_calculator, distance_limit):
        self._distance_calculator = distance_calculator
        self._distance_limit = distance_limit

    def check(self, cluster):
        for node_i, node_j in itertools.combinations(cluster, 2):
            if self._distance_calculator.distance(node_i, node_j) > self._distance_limit:
                return False
        return True


class ClosestLinkage:

    def __init__(self, distance_calculator):
        self._distance_calculator = distance_calculator

    def distance(self, nodes_1, nodes_2) -> float:
        minimal = None
        for node_i, node_j in itertools.product(nodes_1, nodes_2):
            distance = self._distance_calculator.distance(node_i, node_j)

            if distance == 0.0:
                print(f"0 distance {node_i.nid}, {node_j.nid}")

            if minimal is None:
                minimal = distance
            elif distance < minimal:
                minimal = distance
        assert minimal is not None, "Minimal distance calculation failed, check input!"

        return minimal


class Validator:
    def __init__(self, conditions):
        self._validity_conditions = conditions

    def check_validity(self, cluster):
        for condition in self._validity_conditions:
            valid = condition.check(cluster)
            if not valid:
                return False
        return True


def hierarchical_clustering(linkage, validator, nodes):
    clusters = [[i] for i in nodes]

    while len(clusters) > 1:
        minimal = None
        ambiguous = True
        merge = None

        for cluster_1, cluster_2 in itertools.combinations(clusters, 2):
            merged_cluster = cluster_1 + cluster_2
            if validator.check_validity(merged_cluster):
                distance = linkage.distance(cluster_1, cluster_2)
                if minimal is None:
                    ambiguous = False
                    minimal = distance
                    merge = cluster_1, cluster_2
                elif distance < minimal:
                    ambiguous = False
                    minimal = distance
                    merge = cluster_1, cluster_2
                elif distance == minimal:
                    ambiguous = True

        if minimal is not None:
            if ambiguous:
                print(f"WARNING potentially ambiguous situation, {clusters}, {merge}")
            clusters.remove(merge[0])
            clusters.remove(merge[1])
            clusters.append(merge[0] + merge[1])
        else:
            break

    return clusters
