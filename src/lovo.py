import numpy as np
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

def getPartialRMSD(moving, ref, phi):
    """

    :param moving:
    :param ref:
    :param phi:
    :return: rmsd, minimal_movement_indices
    """

    # translate
    moved = moving - (centroid(moving) - centroid(ref))
    # calculate rotational matrix
    rot = quaternion_rotate(moved, ref)
    # calculate aligned moved position
    moved = np.dot(moved, rot)

    _msd = get_msd(moved, ref)
    core_indices = np.argsort(_msd)[0:int(len(moved) * phi)]

    while True:
        rot = quaternion_rotate(moved[core_indices], ref[core_indices])
        moved = np.dot(moved, rot)
        _msd = get_msd(moved, ref)
        new_core_indices = np.argsort(_msd)[0:int(len(moved) * phi)]
        if np.array_equal(new_core_indices, core_indices):
            # Check if converge
            break
        else:
            core_indices = new_core_indices

    return rmsd(moved[core_indices], ref[core_indices]), core_indices


def normalize_msd(msd):
    rmsd = np.mean(msd)
    return 1-np.exp(-np.square(msd / rmsd)/2)


def alovo_traj(coordset, ref, n_clusters=3, similarity_function=normalize_msd):
    """

    :param coordset:
    :param ref_index:
    :param n_clusters:
    :param similarity_function:
    :return: section_indices, coord_msd
    """
    if len(coordset.shape) != 3 or coordset.shape[2] != 3:
        raise ValueError('Invalid Coordset dimensions, must be (X,Y,3), you provided {}'.format(coordset.shape))


    #  Calculate initial msd after iteration
    coord_msd = np.zeros(coordset.shape[:2])
    for i in np.arange(coordset.shape[0]):
        coordset[i] = coordset[i] - (centroid(coordset[i]) - centroid(ref))
        rot = quaternion_rotate(coordset[i], ref)
        coordset[i] = np.dot(coordset[i], rot)
        coord_msd[i] = get_msd(coordset[i], ref)

    # Calculate the b-factor / rmsf



    # Calculate the mean msd across the trajectory
    mean_normalized_coord_msd = similarity_function(np.mean(coord_msd, axis=0))

    kmean = KMeans(n_clusters=n_clusters, random_state=0)

    kmean.fit(np.stack([mean_normalized_coord_msd, np.zeros(len(mean_normalized_coord_msd))]).transpose())
    section_indices = np.array([kmean.labels_ == k for k in np.arange(kmean.n_clusters)])
    section_indices = section_indices[np.argsort([np.mean(mean_normalized_coord_msd[_]) for _ in section_indices])]

    k = 0
    while k < 100:
        for i in np.arange(coordset.shape[0]):
            rot = quaternion_rotate(coordset[i][section_indices[0]], ref[section_indices[0]])
            coordset[i] = np.dot(coordset[i], rot)
            coord_msd[i] = get_msd(coordset[i], ref)
            # normalized_coord_msd[i] = similarity_function(coord_msd[i])
        mean_normalized_coord_msd = similarity_function(np.mean(coord_msd, axis=0))
        kmean.fit(np.stack([mean_normalized_coord_msd, np.zeros_like(mean_normalized_coord_msd)]).transpose())
        new_section_indices = np.array([kmean.labels_ == k for k in np.arange(kmean.n_clusters)])
        new_section_indices = new_section_indices[
            np.argsort([np.mean(mean_normalized_coord_msd[_]) for _ in new_section_indices])]

        if np.array_equal(section_indices[0], new_section_indices[0]):
            section_indices = new_section_indices
            break
        else:
            section_indices = new_section_indices
            k += 1

    # plt.scatter(mean_normalized_coord_msd,)

    return section_indices, np.mean(coord_msd,axis=0)


def rmsf(coordset):
    """
    Calculate the rmsf or the b-factor, which is average of distance to centroid across trajectory per atom, return
    rmsf array with length of number ofs atoms
    :param coordset:
    :return:
    """
    for i in np.arange(1, coordset.shape[0]):
        coordset[i] = coordset[i] - (centroid(coordset[i]) - centroid(coordset[0]))
        rot = quaternion_rotate(coordset[i], coordset[0])
        coordset[i] = np.dot(coordset[i], rot)
    average_coordset = np.mean(coordset,axis=0)
    rmsf = np.average(np.sqrt(np.sum(np.square((coordset - average_coordset)),axis=2)),axis=0)
    return rmsf



def quaternion_transform(r):
    """
    Get optimal rotation
    note: translation will be zero when the centroids of each molecule are the
    same
    """
    Wt_r = makeW(*r).T
    Q_r = makeQ(*r)
    rot = Wt_r.dot(Q_r)[:3, :3]
    return rot


def makeW(r1, r2, r3, r4=0):
    """
    matrix involved in quaternion rotation
    """
    W = np.asarray([
        [r4, r3, -r2, r1],
        [-r3, r4, r1, r2],
        [r2, -r1, r4, r3],
        [-r1, -r2, -r3, r4]])
    return W


def makeQ(r1, r2, r3, r4=0):
    """
    matrix involved in quaternion rotation
    """
    Q = np.asarray([
        [r4, -r3, r2, r1],
        [r3, r4, -r1, r2],
        [-r2, r1, r4, r3],
        [-r1, -r2, -r3, r4]])
    return Q


def quaternion_rotate(X, Y):
    """
    Calculate the rotation

    Parameters
    ----------
    X : ndarray
        (N,D) matrix, where N is points and D is dimension.
    Y:  ndarray
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rot : matrix
        Rotation matrix (D,D)
    """
    N = X.shape[0]
    W = np.asarray([makeW(*Y[k]) for k in range(N)])
    Q = np.asarray([makeQ(*X[k]) for k in range(N)])
    Qt_dot_W = np.asarray([np.dot(Q[k].T, W[k]) for k in range(N)])
    W_minus_Q = np.asarray([W[k] - Q[k] for k in range(N)])
    A = np.sum(Qt_dot_W, axis=0)
    eigen = np.linalg.eigh(A)
    r = eigen[1][:, eigen[0].argmax()]
    rot = quaternion_transform(r)
    return rot


def centroid(X):
    """
    Calculate the centroid from a vectorset X.

    https://en.wikipedia.org/wiki/Centroid
    Centroid is the mean position of all the points in all of the coordinate
    directions.

    C = sum(X)/len(X)

    Parameters
    ----------
    X : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    C : float
        centeroid

    """
    C = np.mean(X, axis=0)
    return C


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------
    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
        Root-mean-square deviation

    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i] - w[i]) ** 2.0 for i in range(D)])
    return np.sqrt(rmsd / N)


def get_msd(V, W) -> np.ndarray:
    """
    Calculate mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------
    V : ndarray
        (N,D) matrix, where N is points and D is dimension.
    W : ndarray
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
        Root-mean-square deviation

    """
    return np.sqrt(np.sum(np.square(V - W), axis=1))
