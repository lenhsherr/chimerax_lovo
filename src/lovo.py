import numpy as np


def getPartialRMSD(moving, ref, phi):
    """

    :param moving:
    :param ref:
    :param phi:
    :return: rmsd, minimal_movement_indices
    """

    moved =  moving - (centroid(moving) - centroid(ref))
    msd_series = [np.sum(msd(moved,ref))]

    rot = quaternion_rotate(moved, ref)

    moved = np.dot(moved, rot)



    _msd = msd(moved, ref)
    truncated_sorted_indices = np.arange(len(moved))
    while not ifConverged(msd_series,cutoff=0.01):

        truncated_sorted_indices = np.argsort(_msd)[0:int(len(moved) * phi)]
        rot = quaternion_rotate(moved[truncated_sorted_indices],ref[truncated_sorted_indices])

        moved = np.dot(moved,rot)

        _msd = msd(moved, ref)

        _msd_l = np.sum(_msd[truncated_sorted_indices])

        msd_series.append(np.sum(_msd_l))

    return rmsd(moved[truncated_sorted_indices],ref[truncated_sorted_indices]), truncated_sorted_indices

def ifConverged(series,cutoff = 0.01):

    if len(series) < 3:
        return False

    if np.max(series)  == np.min(series):
        return True

    else:
        difference =  np.abs((series[-1]-series[-2])/(np.max(series) - np.min(series)))
        return difference < cutoff

def quaternion_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD

    based on doi:10.1016/1049-9660(91)90036-O

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    P : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
    """
    rot = quaternion_rotate(P, Q)
    P = np.dot(P, rot)
    return rmsd(P, Q)

def quaternion_msd(P, Q):
    """
    Rotate matrix P unto Q and calculate the per point RSD

    based on doi:10.1016/1049-9660(91)90036-O

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    P : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
    """
    rot = quaternion_rotate(P, Q)
    P = np.dot(P, rot)
    return msd(P, Q)


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
    X : array
        (N,D) matrix, where N is points and D is dimension.
    Y: array
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
    C = X.mean(axis=0)
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
        rmsd += sum([(v[i] - w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)

def msd(V,W) -> np.array:
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
    return np.sum(np.square(V-W),axis=1)