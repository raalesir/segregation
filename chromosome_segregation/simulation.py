"""
module for  keeping simulation routines
"""
try:
    from regrow import regrow, regrow_biased, regrow_saw
    import  consts
    import  aux
except ImportError:
    from chromosome_segregation.regrow import regrow_biased, regrow, regrow_saw
    # import  aux


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
        # coords = regrow(n, 0, 0, 0, [])
        coords, _, _ = regrow_biased(n, 0, 0, 0, [], w=1, alpha=0.0, k=0)

        coords = np.array(coords).astype(float)
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
        u, c = np.unique(coords, axis=0, return_counts=True)
        #     print(u,c)
        intersections.append(sum([comb(v, 2) for v in c]))



    bins, counts_ = np.unique(intersections, return_counts=True)
    counts_ = [c / sum(counts_) for c in counts_]

    return bins, counts_


def URW_saw(n, n_steps, box):
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
    failed_to_grow = 0
    minx = miny = minz = maxx = maxy = maxz = 0

    for i in range(1,n_steps//50):
        if i % 1000 == 0:
            print("passed %3.1f %%, failed to grow: %3.1f%%"  %(i/n_steps * 100, failed_to_grow/i*100))
        # coords = regrow(n, 0, 0, 0, [])
        coords, _, _ = regrow_saw(n, 0, 0, 0, [], w=1, alpha=0.0, k=0)
        if len(coords) < n:
            failed_to_grow +=1
            # print('skipping, failed to grow')
        else:
            coords = np.array(coords).astype(float)
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
            if minx > min(coords[:, 0]): minx = min(coords[:, 0])
            if miny > min(coords[:, 1]): miny = min(coords[:, 1])
            if minz > min(coords[:, 2]): minz = min(coords[:, 2])

            if maxx < max(coords[:, 0]): maxx = max(coords[:, 0])
            if maxy < max(coords[:, 1]): maxy = max(coords[:, 1])
            if maxz < max(coords[:, 2]): maxz = max(coords[:, 2])

            u, c = np.unique(coords, axis=0, return_counts=True)
            if max(c) >1: print('error')
            # print(box)

        # #     print(u,c)
        # intersections.append(sum([comb(v, 2) for v in c]))

    print(minx, maxx, max(abs(minx), maxx))
    print(miny, maxy, max(abs(miny), maxy))
    print(minz, maxz, max(abs(minz), maxz))


    maxx = max(abs(minx), maxx)
    maxy = max(abs(miny), maxy)
    maxz = max(abs(minz), maxz)
    print('finished estimating limits:', maxx, maxy, maxz)
    max_size = int(max(maxx, maxy, maxz))
    # dx = int(maxx/n_bins); dy = int(maxy/n_bins); dz = int(maxz/n_bins)
    # print("dx=",dx,"dy=", dy,"dz=", dz)
    failed_to_grow = 0

    a, b,c = box[0], box[1], box[2]
    # scale = int(max_size/ min(box)) +1
    # print("box=%s, scale=%i"%(box, scale))
    size_boxes = []

    sx = list(range(a-max(box), a + max_size+4 ))#scale+1))
    sy = list(range(b-max(box), b+max_size+4)) #scale+1))
    sz = list(range(c-max(box), c + max_size+4 ))#scale+1))
    #
    # sx = list(range(0, a + max_size+2))
    # sy = list(range(b-a, len(sx)-b+a))
    # sz = list(range(c-a, len(sx)-c+a))

    print(sx)
    print(sy)
    print(sz)


    for i in range(1, n_steps):
        if i % 10000 == 0:
            print("passed %3.1f %%, failed to grow: %3.1f%%" % (i / n_steps * 100, failed_to_grow / i * 100))
        # coords = regrow(n, 0, 0, 0, [])
        coords, _, _ = regrow_saw(n, 0, 0, 0, [], w=1, alpha=0.0, k=0)
        if len(coords) < n:
            failed_to_grow += 1
            # print('skipping, failed to grow')
        else:
            coords = np.array(coords)#.astype(float)

            # extreme_x = max(abs(min(coords[:, 0])), max(coords[:, 0]))
            extreme_x = (max(coords[:, 0]) - min(coords[:, 0]) +1)
            # extreme_y = max(abs(min(coords[:, 1])), max(coords[:, 1]))
            extreme_y = (max(coords[:, 1]) - min(coords[:, 1]) +1)
            # extreme_z = max(abs(min(coords[:, 2])), max(coords[:, 2]))
            extreme_z = (max(coords[:, 2]) - min(coords[:, 2]) +1)

            box_x = sx.index(extreme_x) # extreme_x - extreme_x % dx) / dx
            box_y = sy.index(extreme_y) #(extreme_y - extreme_y % dy) / dy
            box_z = sz.index(extreme_z)  #(extreme_z - extreme_z % dz) / dz
            box = max(box_x, box_y, box_z)

            size_boxes.append(int(box))

    # bins, counts_ = np.unique(intersections, return_counts=True)
    # counts_ = [c / sum(counts_) for c in counts_]

    # array = aux.list_to_arr(size_boxes+)

    return [el + sx[0] for el in size_boxes] #shifting so first box is always zero....



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

    logging.info('s size is: %s' %(s.shape))
    logging.info('indexes: %s'%indexes)
    logging.info("max overlap %i, min_overlaps  %i" %(max_overlaps, min_overlaps))

    ds = .1

    o_ = min_overlaps
    w_o = 1
    k_o = 0
    collect_s = []
    sweep_number = 0


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
                if np.random.random() < (w_n / w_o) * np.exp(-alpha * (k_o - k_n)) * np.exp(s[o_] - s[n_]):
                    w_o = w_n
                    k_o = k_n
                    o_ = n_
            elif (n_ <= max_overlaps) and (n_ >= min_overlaps):
                n_ = n_ - n_ % grain
                if np.random.random() < (w_n / w_o) * np.exp(-alpha * (k_o - k_n)) * np.exp(s[o_] - s[n_]):
                    w_o = w_n
                    k_o = k_n
                    o_ = n_

            counts[o_] += 1
            s[o_] += ds
        # print(o_)

        t = counts[indexes]
        mean = sum(t) / len(t)
        print('sweep number',sweep_number, 'mean=',round(mean, 2), 'max=',max(t), 'min=',min(t), end='\r')
        # print(t, end='\r')

        if (max(t) / mean - 1 < flatness) & (1 - min(t) / mean < flatness):
            counts = 0 * counts
            collect_s.append(s.copy())
            ds = ds / 1.5
            logging.info("ds=%e, ds_min=%f, sweep number=%i"%(ds, ds_min, sweep_number))

    return collect_s, sweep_number



def WL_saw(n, max_overlaps, min_overlaps=0, grain=1, exclude=(), alpha=0, sweep_length=1000, ds_min=0.0000001,
       flatness=0.3):
    """
    WL  procedure for calculation of SAW distribution based on size....

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

    logging.info('s size is: %s' % (s.shape))
    logging.info('indexes: %s' % indexes)
    logging.info("max overlap %i, min_overlaps  %i" % (max_overlaps, min_overlaps))

    ds = .1

    o_ = min_overlaps
    w_o = 1
    k_o = 0
    collect_s = []
    sweep_number = 0

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
                if np.random.random() < (w_n / w_o) * np.exp(-alpha * (k_o - k_n)) * np.exp(s[o_] - s[n_]):
                    w_o = w_n
                    k_o = k_n
                    o_ = n_
            elif (n_ <= max_overlaps) and (n_ >= min_overlaps):
                n_ = n_ - n_ % grain
                if np.random.random() < (w_n / w_o) * np.exp(-alpha * (k_o - k_n)) * np.exp(s[o_] - s[n_]):
                    w_o = w_n
                    k_o = k_n
                    o_ = n_

            counts[o_] += 1
            s[o_] += ds
        # print(o_)

        t = counts[indexes]
        mean = sum(t) / len(t)
        print('sweep number', sweep_number, 'mean=', round(mean, 2), 'max=', max(t), 'min=', min(t), end='\r')
        # print(t, end='\r')

        if (max(t) / mean - 1 < flatness) & (1 - min(t) / mean < flatness):
            counts = 0 * counts
            collect_s.append(s.copy())
            ds = ds / 1.5
            logging.info("ds=%e, ds_min=%f, sweep number=%i" % (ds, ds_min, sweep_number))

    return collect_s, sweep_number
