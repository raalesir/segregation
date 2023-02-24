"""
module for  keeping simulation routines
"""
try:
    from regrow import regrow, regrow_biased, regrow_saw
    import  consts
    import  aux
except ImportError:
    from chromosome_segregation.regrow import regrow_biased, regrow, regrow_saw
    from  chromosome_segregation import  aux


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

    sx = list(range(a-max(box), a + max_size+10 ))#scale+1))
    sy = list(range(b-max(box), b+max_size+10)) #scale+1))
    sz = list(range(c-max(box), c + max_size+10 ))#scale+1))
    #
    # sx = list(range(0, a + max_size+2))
    # sy = list(range(b-a, len(sx)-b+a))
    # sz = list(range(c-a, len(sx)-c+a))

    logging.info('x-array is: %s' %sx)
    logging.info('y-arrya is: %s' %sy)
    logging.info('z-array is: %s' %sz)
    print('x-array is: %s' %sx)
    print('y-array is: %s' %sy)
    print('z-array is: %s' %sz)



    for i in range(1, n_steps):
        if i % 10000 == 0:
            logging.info("passed %3.1f %%, failed to grow: %3.1f%%" % (i / n_steps * 100, failed_to_grow / i * 100))
        # coords = regrow(n, 0, 0, 0, [])
        coords, _, _ = regrow_saw(n, 0, 0, 0, [], w=1, alpha=0.0, k=0)
        if len(coords) < n:
            failed_to_grow += 1
            # print('skipping, failed to grow')
        else:
            coords = np.array(coords)#.astype(float)

            # extreme_x = max(abs(min(coords[:, 0])), max(coords[:, 0]))
            extreme_x = (max(coords[:, 0]) - min(coords[:, 0]) +1)  # otherwise the volume could be  zero!
            # extreme_y = max(abs(min(coords[:, 1])), max(coords[:, 1]))
            extreme_y = (max(coords[:, 1]) - min(coords[:, 1]) +1)
            # extreme_z = max(abs(min(coords[:, 2])), max(coords[:, 2]))
            extreme_z = (max(coords[:, 2]) - min(coords[:, 2]) +1)

            # box_x = sx.index(extreme_x) # extreme_x - extreme_x % dx) / dx
            # box_y = sy.index(extreme_y) #(extreme_y - extreme_y % dy) / dy
            # box_z = sz.index(extreme_z)  #(extreme_z - extreme_z % dz) / dz
            # box = max(box_x, box_y, box_z)

            box = aux.get_box(sx=sx, sy=sy, sz=sz, l=(extreme_x, extreme_y, extreme_z))
            size_boxes.append(int(box))

    # bins, counts_ = np.unique(intersections, return_counts=True)
    # counts_ = [c / sum(counts_) for c in counts_]

    # array = aux.list_to_arr(size_boxes+)

    return [el + sx[0] for el in size_boxes] #shifting so first box is always zero....



def WL_saw(n, indexes, sweep_length=1000, ds_min=0.0000001, flatness=0.3,
        decrease = 2.0, scale_alpha=2.0, shift_alpha=0.0):
    """
    WL  procedure for calculation of SAW enropies inside matroshkas
    """

    s = np.zeros(len(indexes[0]))
    counts = np.zeros(s.shape)
    indexes_ = list(range(len(indexes[0])))
    wn = np.zeros(s.shape)

    logging.info('s size is: %s' % (s.shape))
    logging.info('counts size is: %s' % (counts.shape))

    logging.info('indexes: %s' % indexes_)
    #     logging.info("max overlap %i, min_overlaps  %i" %(max_overlaps, min_overlaps))


    ds = .1
    box_o = -1
    w_o = 1
    k_o = 0
    collect_s = []
    sweep_number = 0
    alpha = 0.0
    frozen =0
    logging.info('scale alpha is: %f'%scale_alpha)
    logging.info('shift alpha is: %f' %shift_alpha)

    while ds > ds_min:
        sweep_number += 1
        failed_to_grow = 0
        out_of_range = 0
        alpha_o = 0
        for i in range(sweep_length):
            alpha = (np.random.rand()-shift_alpha)*scale_alpha #choice((0,1,2))

            coords_n, w_n, k_n, prob = regrow_saw(n, 0, 0, 0, [], w=1, alpha=alpha, k=0, prob=1)
            if len(coords_n) < n:
                failed_to_grow += 1
            # print('skipping, failed to grow')
            else:
                coords_n = np.array(coords_n)  # .astype(float)

                # u, c = np.unique(coords_n, axis=0, return_counts=True)
                # if coords_n.shape != u.shape:
                #     print('wrong SAW')

                extreme_x = (max(coords_n[:, 0]) - min(coords_n[:, 0]) ) # otherwise the volume could be  zero!
                # if extreme_x == 0:   extreme_x =1
                extreme_y = (max(coords_n[:, 1]) - min(coords_n[:, 1]))
                # if extreme_y == 0:  extreme_y =1
                extreme_z = (max(coords_n[:, 2]) - min(coords_n[:, 2]))
                # if extreme_z == 0:   extreme_z =1

                box_n = aux.get_box(sx=indexes[0], sy=indexes[1], sz=indexes[2], l=(extreme_x, extreme_y, extreme_z))

                if (box_n == None):
                    box_n = indexes_[-1]
                    out_of_range += 1
                # print((w_n / w_o) * np.exp(-alpha * (k_o - k_n)), k_n, w_n)
                # if np.random.random() < np.exp(s[box_o] - s[box_n]):
                if frozen > sweep_length/5:
                    print('resetting frozen')
                if (np.random.random() < (w_n / w_o) * np.exp(-alpha_o*k_o
                    +alpha* k_n) * np.exp(s[box_o] - s[box_n])) or (frozen > sweep_length/5):
                    
                    alpha_o = alpha
                    box_o = box_n
                    #if np.exp(-alpha*k_n)/w_n != wn[box_o]:
                    #    wn[box_o] = (wn[box_o] + np.exp(-alpha*k_n)/w_n) / 2

                    frozen =0
                    w_o=w_n
                    k_o=k_n
                else:
                    if box_o == box_n:
                       pass
#                        frozen +=1
                # print(extreme_x, extreme_y, extreme_z, box_o)
                counts[box_o] += 1
                s[box_o] += ds
                # if box_o == 1:
                #   print(w_o, wn[box_o])

        # print(o_)

        t = counts[indexes_]
        mean = sum(t) / len(t)
        # print('sweep number', sweep_number, 'mean=', round(mean, 2), 'max=', max(t), 'min=', min(t), end='\r')
        # print(t, end='\r')
        logging.info("visits %s" %t)
        logging.info("S %s" %s)
        #print(t, s)
        # print(wn)
        logging.info('failed to grow rate: %2.0f%%, out of range: %2.0f%% '%(100.0*failed_to_grow/sweep_length, 100.0*out_of_range/sweep_length))
        if (max(t) / mean - 1 < flatness) & (1 - min(t) / mean < flatness):
            counts = 0 * counts
            collect_s.append(s.copy())
            ds = ds / decrease
            logging.info("ds=%e, ds_min=%f, sweep number=%i" % (ds, ds_min, sweep_number))

    return collect_s, sweep_number



def WL_saw_mpi(n, indexes, sweep_length=1000, ds_min=0.0000001, flatness=0.3,
        decrease = 2.0, scale_alpha=2.0, shift_alpha=0.0):
    """
    WL  procedure for calculation of SAW enropies inside matroshkas
    """

    s = np.zeros(len(indexes[0]))
    counts = np.zeros(s.shape)
    indexes_ = list(range(len(indexes[0])))
    wn = np.zeros(s.shape)

    logging.info('s size is: %s' % (s.shape))
    logging.info('counts size is: %s' % (counts.shape))

    logging.info('indexes: %s' % indexes_)
    #     logging.info("max overlap %i, min_overlaps  %i" %(max_overlaps, min_overlaps))


    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    ds = 0.1
    box_o = -1
    w_o = 1
    k_o = 0
    probability_old = 1
    probability_native_old = 1
    collect_s = []
    sweep_number = 0
    alpha = 0.0
    frozen =0
    logging.info('scale alpha is: %f'%scale_alpha)
    logging.info('shift alpha is: %f' %shift_alpha)
    alpha_o = 0
    
    begin = True
    #alphas = np.linspace(scale_alpha, 0, size)
    # nonuniform setting for alpha
    #t  = np.log(np.linspace(1,20 ,size))[::-1]
    #t = t - np.min(t)
    #alphas = t/np.max(t)*scale_alpha

    coords_n, w_n, k_n, probability, probability_native_old  = regrow_saw(n, 0, 0, 0, [], w=1, alpha=alpha, k=0, prob=1, prob_native=1)
    
    while ds > ds_min:
        sweep_number += 1
        failed_to_grow = 0
        out_of_range = 0
        zero_box = 0
    #    if rank == 0:
    #      s = s/size
    #      counts = counts/size
 
        s = comm.bcast(s, root=0)
        counts = comm.bcast(counts, root=0)

        for i in range(sweep_length):
            
            #if size < 1:
            
            #alpha = (np.random.rand()-shift_alpha)*scale_alpha 
            alpha = scale_alpha #np.random.choice((0,scale_alpha))
           # else:
            #    alpha = alphas[rank]

            coords_n, w_n, k_n, probability, probability_native  = regrow_saw(n, 0, 0, 0, [], w=1, alpha=alpha, k=0, prob=1, prob_native=1)
            if len(coords_n) < n:
                failed_to_grow += 1
            # print('skipping, failed to grow')
            else:
                #if probability != probability_old:
                #   logging.error('probabilities are:  %e, %e' %(probability,probability_old)) 
 
                coords_n = np.array(coords_n)  # .astype(float)

                # u, c = np.unique(coords_n, axis=0, return_counts=True)
                # if coords_n.shape != u.shape:
                #     print('wrong SAW')

                extreme_x = (max(coords_n[:, 0]) - min(coords_n[:, 0]) ) # otherwise the volume could be  zero!
                # if extreme_x == 0:   extreme_x =1
                extreme_y = (max(coords_n[:, 1]) - min(coords_n[:, 1]))
                # if extreme_y == 0:  extreme_y =1
                extreme_z = (max(coords_n[:, 2]) - min(coords_n[:, 2]))
                # if extreme_z == 0:   extreme_z =1

                box_n = aux.get_box(sx=indexes[0], sy=indexes[1], sz=indexes[2], l=(extreme_x, extreme_y, extreme_z))
                #print(rank, box_n)
                if (box_n == None):
                    box_n = indexes_[-1]
                    out_of_range += 1
                elif box_n == 0:
                    zero_box +=1
                    #logging.info('configuration probability for 0th box for walker %i is: %e' % (rank, probability))
                # print((w_n / w_o) * np.exp(-alpha * (k_o - k_n)), k_n, w_n)
                # if np.random.random() < np.exp(s[box_o] - s[box_n]):
                if frozen > sweep_length/2:
                    logging.info("resetting frozen. Walker %i could not get out of box=%i" %(rank, box_o))
               #if (np.random.random() < (w_n / w_o) * np.exp(-alpha_o*k_o
               #     +alpha* k_n) * np.exp(s[box_o] - s[box_n])) or (frozen > sweep_length/5):
                #if (np.random.random() <  (probability_old/probability) * np.exp(s[box_o] - s[box_n])) or (frozen > sweep_length/2):
                if (np.random.random() < (probability_native/probability_native_old)) or (frozen > sweep_length/2):
                #if (np.random.random() < 1)or (frozen > sweep_length/10):
                #if (np.random.random() < np.exp(s[box_o] - s[box_n]))or (frozen > sweep_length/10):
                    alpha_o = alpha
                    box_o = box_n
                    probability_old = probability
                    probability_native_old = probability_native
                    #if np.exp(-alpha*k_n)/w_n != wn[box_o]:
                    #    wn[box_o] = (wn[box_o] + np.exp(-alpha*k_n)/w_n) / 2

                    frozen =0
                    w_o=w_n
                    k_o=k_n
                    
                    #logging.info('probabilities should be equal: %e and %e'%(probability, np.exp(-alpha_o*k_o)/w_o))
                else:
    #                if box_o == box_n:
                  #pass
                  frozen +=1
                # print(extreme_x, extreme_y, extreme_z, box_o)
                counts[box_o] += 1
                s[box_o] += ds
              
                # if box_o == 1:
                #   print(w_o, wn[box_o])

        # print(o_)
        #counts = comm.reduce(counts, op=MPI.SUM, root=0)
        #s = comm.reduce(s, op=MPI.SUM, root=0)
        tmp_counts = comm.gather(counts, root=0)
        tmp_s = comm.gather(s, root=0)

        logging.info("fraction of generated confs for 0th box: %f, absolute number %i, rank %i" % (1.0*zero_box/(sweep_length-failed_to_grow), zero_box, rank))
        if rank ==0:
          counts = np.max(np.array(tmp_counts), axis=0)
          s = np.max(np.array(tmp_s), axis=0)
          #logging.info(tmp_counts)
          t = counts[indexes_]
          mean = sum(t) / len(t)
        # print('sweep number', sweep_number, 'mean=', round(mean, 2), 'max=', max(t), 'min=', min(t), end='\r')
        # print(t, end='\r')
          logging.info("visits %s" %t)
          logging.info("S %s" %s)
        #print(t, s)
        # print(wn)
          logging.info('failed to grow rate: %2.0f%%, out of range: %2.0f%% '%(100.0*failed_to_grow/sweep_length, 100.0*out_of_range/sweep_length))
          if (max(t) / mean - 1 < flatness) & (1 - min(t) / mean < flatness):
              counts = 0 * counts
              collect_s.append(s.copy())
              ds = ds / decrease
              logging.info("ds=%e, ds_min=%f, sweep number=%i" % (ds, ds_min, sweep_number))
        ds = comm.bcast(ds, root=0)

    return collect_s, sweep_number

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




def WL_CM_polymers(n, max_overlaps, min_overlaps=0, grain=1, exclude=(), alpha=0, sweep_length=1000, ds_min=0.0000001,
       flatness=0.3):
    """
    Calculates S(R_cm) -- entropy of a center of mass distance between two polymers, confined in the tube.
    We fix two polymers at  (0,0,0). First we generate SAW for  one polymer, secondly -- the second, such that it
    does not overlap with the first one. Both are inside the tube.
    Calculate cm-distance and accumulate S(R_cm).
    Regrow the first chain, so it does not  overlap with the second etc.

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
            coords_n, w_n, k_n = regrow_biased(n, 0, 0, 0, [], w=1, alpha=0, k=0)

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
