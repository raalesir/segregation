"""
module for  keeping simulation routines
"""

from regrow import regrow, regrow_biased
from scipy.special import comb
import  numpy as np
import logging

def URW(n, n_steps):
    """
    Uniform random walks (URW) in the configuration space of the ring lattice polymer

    :param n: number o beads (monomers)
    :type n: int
    :param n_steps: number of conformations
    :type n_steps: int
    :return: Numpy array of intersection numbers [0,n_max] and corresponding counts
    :rtype: tuple
    """
    intersections = []
    for i in range(n_steps):
        if i % 10000 == 0:
            logging.info("passed %3.1f %%"  %(i/n_steps * 100))
        coords = regrow(n, 0, 0, 0, [])
        coords = np.array(coords).T.astype(float)
        #         print(coords, coords.shape, "\n")
        #         sys.exit()
        #     scatter1.x = lines1.x = coords[0,:]
        #     scatter1.y = lines1.y = coords[1,:]
        #     scatter1.z = lines1.z = coords[2,:]
        #     sleep(.2)
        #         if not np.sum(np.abs(np.diff(coords, axis=1))) == coords.shape[1]-1:
        #             print('failed  to build')
        #             print(coords, "\n")
        #             sys.exit()
        #     print(coords.T)
        u, c = np.unique(coords.T, axis=0, return_counts=True)
        #     print(u,c)
        intersections.append(sum([comb(v, 2) for v in c]))

    bins, counts_ = np.unique(intersections, return_counts=True)
    counts_ = [c / sum(counts_) for c in counts_]

    return bins, counts_



def WL(n, max_overlaps, min_overlaps=0, grain=1, exclude=(), alpha=0, sweep_length=1000, ds_min=0.0000001,
       flatness=0.3):
    """
    WL  procedure for calculation of overlapping entropy


    """
    if min_overlaps > 0:
        #         max_overlaps += 1
        counts = np.zeros(max_overlaps + 2)
        s = np.zeros(max_overlaps + 2)
        indexes = list(set([ind - ind % grain for ind in range(min_overlaps, max_overlaps + 1) if ind not in exclude]))
        indexes += [max_overlaps + 1]  # adding index for max_overlaps+ intersections


    else:
        counts = np.zeros(max_overlaps + 1)
        s = np.zeros(max_overlaps + 1)
        indexes = list(set([ind - ind % grain for ind in range(min_overlaps, max_overlaps + 1) if ind not in exclude]))

    print('s size is: ', s.shape)
    logging.info('indexes: %s'%indexes)
    logging.info("max overlap %i, min_overlaps  %i" %(max_overlaps, min_overlaps))

    ds = .1

    o_ = min_overlaps
    w_o = 1
    k_o = 0
    collect_s = []
    sweep_number = 0

    #     indexes = [ind for ind in range(max_overlaps) if ind not in exclude]


    while ds > ds_min:
        sweep_number += 1
        for i in range(sweep_length):

            #             coords_n = regrow(n, 0,0,0, [])
            coords_n, w_n, k_n = regrow_biased(n, 0, 0, 0, [], w=1, alpha=alpha, k=0)

            coords_n = np.array(coords_n)
            u, c = np.unique(coords_n, axis=0, return_counts=True)
            #             c = Counter(coords_n).values()
            n_ = int(sum([comb(v, 2) for v in c]))
            #             n_ = sum([el for el in c if  el>1])
            #             n_ = n_ - n_%grain

            if (n_ > max_overlaps) and (min_overlaps > 0):
                n_ = max_overlaps + 1
                if np.random.rand() < (w_n / w_o) * np.exp(-alpha * (k_o - k_n)) * np.exp(s[o_] - s[n_]):
                    w_o = w_n
                    k_o = k_n
                    o_ = n_
            elif (n_ <= max_overlaps) and (n_ >= min_overlaps):
                n_ = n_ - n_ % grain
                if np.random.rand() < (w_n / w_o) * np.exp(-alpha * (k_o - k_n)) * np.exp(s[o_] - s[n_]):
                    w_o = w_n
                    k_o = k_n
                    o_ = n_

            counts[o_] += 1
            s[o_] += ds
        # print(o_)

        t = counts[indexes]
        #         print(counts)

        #         t = counts[counts>0]
        #     print(t)
        mean = sum(t) / len(t)
        print('sweep number',sweep_number, 'mean=',round(mean, 2), 'max=',max(t), 'min=',min(t), end='\r')

        #     print(max(t)/mean -1, 1- min(t)/mean )
        if (max(t) / mean - 1 < flatness) & (1 - min(t) / mean < flatness):
            # print(counts)
            counts = 0 * counts
            # print(repr(s))
            collect_s.append(s.copy())
            ds = ds / 1.5
            # print('')
            logging.info("ds=%e, ds_min=%f, sweep number=%i"%(ds, ds_min, sweep_number))

    return collect_s, sweep_number
