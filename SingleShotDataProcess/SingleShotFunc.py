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
    num_batch = min(10000, num_points)
    labels = kmeans.labels_[:num_batch]
    x_batch = x[:num_batch]
    y_batch = y[:num_batch]
    sigmas = np.zeros_like(centers)
    for i in range(len(centers)):
        sigmas[i, 0] = np.std(x_batch[labels == i])
        sigmas[i, 1] = np.std(y_batch[labels == i])
    return centers, sigmas


def getCenterHeights(X, Y, H, centers):
    num_blob = len(centers)
    heights = np.zeros((num_blob,))
    for i in range(num_blob):
        x_c = centers[i, 0]
        y_c = centers[i, 1]
        ind_x = np.argmin((np.unique(X) - x_c) ** 2)
        ind_y = np.argmin((np.unique(Y) - y_c) ** 2)
        heights[i] = H[ind_y, ind_x]

    return heights