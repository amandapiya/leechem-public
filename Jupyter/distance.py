import numpy as np;

def distance(x1, y1, z1, x2, y2, z2):
    return np.sqrt(np.power(x1 - x2, 2) + np.power(y1 - y2, 2) + np.power(z1 - z2, 2));

