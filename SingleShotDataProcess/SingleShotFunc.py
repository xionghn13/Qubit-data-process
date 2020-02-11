import numpy as np
from sklearn.cluster import MiniBatchKMeans, KMeans


def getBlobCenters(x, y, num_blob):
    num_points = len(x)
    matrix = np.array([x, y]).transpose()
    if num_points > 10e3:
        kmeans = MiniBatchKMeans(n_clusters=num_blob).fit(matrix)
    else:
        kmeans = KMeans(n_clusters=num_blob).fit(matrix)
    centers = kmeans.cluster_centers_
    return centers


def getCenterHeights(X, Y, H, centers):

    return centers