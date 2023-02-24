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
import math

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
    def add_key(k, ncs):
      if k in ncs.keys():
        ncs[k]+=1
      else:
        ncs[k] = 1
      return ncs

    n_boxes = int(n*1.2)
    
    #bins = np.linspace(1.15,  0.2,  n_boxes)
    #bins = np.linspace(0.95,  0.05,  n_boxes)
    #bins = np.linspace(0.4, 0.02,  n_boxes)
    #bins = np.linspace(1.0, 0.3,  n_boxes)
    #bins = np.linspace(1.15, 0.02,  n_boxes)
    s = np.zeros(n_boxes)
    counts = np.zeros(n_boxes)    
    logging.info('number of contact boxes is: %i' %(n_boxes))
    #s = np.zeros(len(indexes[0]))
    #counts = np.zeros(s.shape)
    indexes_ = list(range(20,len(s)-35))
    #indexes_ = list(range(10))
    wn = np.zeros(s.shape)

    logging.info('s size is: %s' % (s.shape))
    logging.info('counts size is: %s' % (counts.shape))

    logging.info('indexes: %s' % indexes_)
    #     logging.info("max overlap %i, min_overlaps  %i" %(max_overlaps, min_overlaps))


    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank ==0:
      number_of_contacts_global = {}
      list_probs_global = {}
      for i in range(n_boxes):
        list_probs_global[i] = []
    ds = 0.1
    ds0 = ds
    box_o = -1
    w_o = 1
    k_o = 0
    probability_old = 1
    probability_nativei_old = 1
    collect_s = []
    sweep_number = 0
    alpha = 0.0
    frozen =0
    logging.info('scale alpha is: %f'%scale_alpha)
    #logging.info('shift alpha is: %f' %shift_alpha)
    alpha_o = 0
    
    begin = True
    #alphas = np.linspace(scale_alpha, 0, size)
    # nonuniform setting for alpha
    #t  = np.log(np.linspace(1,20 ,size))[::-1]
    #t = t - np.min(t)
    #alphas = t/np.max(t)*scale_alpha
    initial_alpha = 2.0
    box_o = len(s) +1
    nc_bin_o  = n_boxes +1
    while nc_bin_o >  n_boxes: #indexes_[-1]:
      coords_n, w_n, k_n, probability_old, probability_native_old  = regrow_saw(n, 0, 0, 0, [], w=1, alpha=scale_alpha, k=0, prob=1, prob_native=1,
box = indexes_[-1], indexes=indexes)
    
      coords_n = np.array(coords_n)  # .astype(float)
      nc = aux.calculate_number_of_contacts(coords_n)
      nc_bin_o = nc#np.digitize(1.0*nc/n, bins)

      extreme_x = (max(coords_n[:, 0]) - min(coords_n[:, 0]) ) # otherwise the volume could be  zero!
      extreme_y = (max(coords_n[:, 1]) - min(coords_n[:, 1]))
      extreme_z = (max(coords_n[:, 2]) - min(coords_n[:, 2]))
      box_o = aux.get_box(sx=indexes[0], sy=indexes[1], sz=indexes[2], l=(extreme_x, extreme_y, extreme_z))
      if not box_o: box_o = len(s)
    logging.info("first contact box, %i, first space box is %i" %(nc_bin_o, box_o))
    aver_dist = - np.log(probability_old)/scale_alpha/n
    #probability_old =  0.0# np.exp(-scale_alpha * aver_dist)
    aver_probs = np.zeros((2,len(s)))

    bin_count = 0  # keeping statistics on last bin  visit  dynamics

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
        aver_probs = comm.bcast(aver_probs, root=0)
        
        list_probs = {}
        for i in range(n_boxes):
          list_probs[i] = []

        number_of_contacts = {}

        for i in range(sweep_length):
            
            #if size < 1:
            
            #alpha = (np.random.rand()-shift_alpha)*scale_alpha 
            alpha = scale_alpha #np.random.choice((0,scale_alpha))
            #alpha = np.random.choice((-4.5 ,0.0, -2.0, 2.5))
           # else:
            #    alpha = alphas[rank]

            coords_n, w_n, k_n, probability, probability_native  = regrow_saw(n, 0, 0, 0, [], w=1, alpha=alpha, k=0, prob=1, prob_native=1, box = indexes_[-1], indexes=indexes)
            if len(coords_n) < n:
                failed_to_grow += 1
            # print('skipping, failed to grow')
            else:
                #print(probability)
                #if probability <  probability_old:
                #logging.error('probability:  %e, probability_native  %e' %(probability,probability_native)) 
 
                coords_n = np.array(coords_n)  # .astype(float)

                #u, c = np.unique(coords_n, axis=0, return_counts=True)
                if False:
                #if coords_n.shape != u.shape:
                   print('wrong SAW')
                   #pass
                else:
                #aver_dist = - np.log(probability)/alpha/n
                #probability =/ n#  np.exp(-2*alpha * aver_dist)
                  nc = k_n# aux.calculate_number_of_contacts_new(coords_n)
                  nc_bin_n = nc#np.digitize(1.0*nc/n, bins)
                  if nc != k_n:
                     logging.error("numbers of contacts do not match! From regrow: %i, from whole chain %i"%(k_n, nc))
                  # accumulating probabilities for each box
                  #list_probs 
                #print(rank, box_n)
                  if (nc_bin_n  > indexes_[-1]):
                      nc_bin_n  = indexes_[-1]
                      out_of_range += 1
                  elif nc_bin_n < indexes_[0]:
                    nc_bin_n = indexes_[0]
                    #extreme_x = (max(coords_n[:, 0]) - min(coords_n[:, 0]) ) # otherwise the volume could be  zero!
                    #extreme_y = (max(coords_n[:, 1]) - min(coords_n[:, 1]))
                    #extreme_z = (max(coords_n[:, 2]) - min(coords_n[:, 2]))
                    #box_n = aux.get_box(sx=indexes[0], sy=indexes[1], sz=indexes[2], l=(extreme_x, extreme_y, extreme_z))
                    #pass
                  elif nc_bin_n == 0:
                      zero_box +=1
                      #nc = aux.calculate_number_of_contacts(coords_n)
                      logging.info('spacial box number is  %i for contact bin %i '%(box_n, nc_bin_n)) 
                      #logging.info('configuration probability for 0th box for walker %i is: %e, aver_dist is: %f' % (rank, probability, aver_dist))
                  #elif box_n == indexes_[-2]:
                  #    pass # print(box_n, extreme_x, extreme_y, extreme_z)
                  elif nc_bin_n == indexes_[-1]:
                      logging.info('contact bin %i, nc is %i, prob %e '%(nc_bin_n,nc, probability)) #if np.random.random() < .01:
                  #    pass #rint(box_n, extreme_x, extreme_y, extreme_z)
                  elif nc_bin_n == 1:
                      #pass
                      #nc = aux.calculate_number_of_contacts(coords_n)
                      logging.info('contact bin %i, nc is %i, prob %e '%(nc_bin_n,nc, probability)) #if np.random.random() < .01:
                      #logging.info('spacial box number is  %i for contact bin %i '%(box_n, nc_bin_n)) 
                      #logging.info('configuration probability for 1th box for walker %i is: %e, aver_dist is: %f' % (rank, probability, aver_dist))
                    #print(repr(coords_n))
                    #sys.exit()
                  #elif box_n == 10:
                  #    pass
                      #if np.random.random() < .1:
                      #   logging.info('configuration probability for 10th box for walker %i is: %e, aver_dist is: %f' % (rank, probability, aver_dist))
                  elif nc_bin_n == 5:
                      pass
                      #nc = aux.calculate_number_of_contacts(coords_n)
                      #logging.info('number of contacts for box %i is %i'%(box_n, nc)) 
                      #logging.info('contact bin %i, nc is %i, prob %e '%(nc_bin_n,nc, probability)) #if np.random.random() < .01:
                      #   logging.info('configuration probability for 5th box for walker %i is: %e, aver_dist is: %f' % (rank, probability, aver_dist))
                  elif nc_bin_n == 6:
                      #logging.info('contact bin %i, nc is %i, prob %e '%(nc_bin_n,nc, probability)) #if np.random.random() < .01:
                      pass
                  #elif box_n > [indexes_[-1]:
                  #else:
                      # accumulating probabilities for each box
                  #    print(box_n, 'wrong box')     #print(repr(coords_n))
                      #sys.exit()
                # print((w_n / w_o) * np.exp(-alpha * (k_o - k_n)), k_n, w_n)
                # if np.random.random() < np.exp(s[box_o] - s[box_n]):
                #if frozen > sweep_length/100:
                #    logging.info("resetting frozen. Walker %i could not get out of box=%i" %(rank, box_o))
               #if (np.random.random() < (w_n / w_o) * np.exp(-alpha_o*k_o
               #     +alpha* k_n) * np.exp(s[box_o] - s[box_n])) or (frozen > sweep_length/5):
                  #if (np.random.random() <  (probability_old/probability) * np.exp(s[box_o] - s[box_n])) or (frozen > sweep_length/2):
                #if (np.random.random() < (probability_native/probability_native_old)*(probability_old/probability)) or (frozen > sweep_length/10):
                  #if (np.random.random() < 1):
                #if (np.random.random() < (probability_old/probability))or (frozen > sweep_length/100):
                #if (np.random.random() < (probability_native/probability_native_old))or (frozen > sweep_length/100):
                  #print(s, box_o, box_n,  probability_old, probability)
                  #print(np.exp(s[box_o] - s[box_n]))
                  if aver_probs[0, nc_bin_n] == 0:
                       res = 1
                       #logging.info('res =1 for box=%i'%nc_bin_n)
                  else:
                       res = aver_probs[0, nc_bin_o]/aver_probs[0, nc_bin_n]
                  res = 1
                  #if (np.random.random() < np.exp(s[box_o] - s[box_n]) ) and  (box_n <= indexes_[-1]):
                  #if (np.random.random() < np.exp(s[nc_bin_o] - s[nc_bin_n])*probability_old/probability):
                  rnd = np.random.random()
                  if (rnd < np.exp(s[nc_bin_o] - s[nc_bin_n])*probability_old/probability):
                      if s[nc_bin_o] - s[nc_bin_n] >  5:
                        logging.info('rare bin visited rnd=%f, from %5.2f to %5.2f, from box %i to box %i'%(rnd, s[nc_bin_o] , s[nc_bin_n], nc_bin_o,     nc_bin_n))
                      if s[nc_bin_o] - s[nc_bin_n] < - 5:
                        logging.info('low prob event rnd=%f, from %5.2f to %5.2f, from box %i to box %i'%(rnd, s[nc_bin_o] , s[nc_bin_n], nc_bin_o, nc_bin_n))
                      if nc_bin_n ==  indexes_[-1]:
                         logging.debug("entered the last bin %f with %f from bin %f having %f"%(nc_bin_n, s[nc_bin_n], nc_bin_o, s[nc_bin_o]))
                         bin_count  = counts[nc_bin_n]
                      if nc_bin_o == indexes_[-1]:
                         bin_count = counts[nc_bin_o] - bin_count                         
                         logging.debug("left the last bin %f with %f to  bin %f having %f after %i steps"%(nc_bin_o, s[nc_bin_o], nc_bin_n, s[nc_bin_n], bin_count))
                         
                  #if (np.random.random() < 2 ):#np.exp(s[nc_bin_o] - s[nc_bin_n]))or (frozen > sweep_length/10):
                  #if (box_n <= indexes_[-1]):
                  #if (np.random.random() < (res) )or (frozen > sweep_length/10):
                      alpha_o = alpha
                      nc_bin_o = nc_bin_n
                      probability_old = probability
                      probability_native_old = probability_native

                      frozen =0
                      w_o=w_n
                      k_o=k_n
                      
                      #en = aux.get_energy(coords_n)
                      #list_probs[nc_bin_n].append(probability) 
                      #number_of_contacts = add_key(nc, number_of_contacts)

                      #if aver_probs[1,nc_bin_o] < 100:
                      #  aver_probs[0, nc_bin_o]  = (aver_probs[0,nc_bin_o]*aver_probs[1,nc_bin_o] + probability_old)/(aver_probs[1,nc_bin_o] +1)
                      #  aver_probs[1, nc_bin_o] += 1
                      #else:
   
                      #  if (((probability_old/aver_probs[0,nc_bin_o]) >0.1) and ((probability_old/aver_probs[0,nc_bin_o]) <10) ): 
                      #    aver_probs[0, nc_bin_o]  = (aver_probs[0,nc_bin_o]*aver_probs[1,nc_bin_o] + probability_old)/(aver_probs[1,nc_bin_o] +1)
                      #    aver_probs[1, nc_bin_o] += 1
                      #  else:
                      #    aver_probs[0, nc_bin_o]  = (aver_probs[0,nc_bin_o]*aver_probs[1,nc_bin_o] + aver_probs[0,nc_bin_o])/(aver_probs[1,nc_bin_o] +1)
                      #    aver_probs[1, nc_bin_o] += 1
 
                        
                      #counts[box_o] += 1
                      #s[box_o] += ds
                    #logging.info('probabilities should be equal: %e and %e'%(probability, np.exp(-alpha_o*k_o)/w_o))
                  else:
                    pass
                    #if nc_bin_n == nc_bin_o:
                    #  print(10*'***')
                    #  print(nc_bin_n, s[nc_bin_n])#pass
                    #frozen +=1
                # print(extreme_x, extreme_y, extreme_z, box_o)
                  counts[nc_bin_o] += 1
                  s[nc_bin_o] += ds
                  
                  number_of_contacts = add_key(nc, number_of_contacts)
                  #list_probs = add_key(probability, list_probs)
                  list_probs[nc_bin_n].append(probability) 
                # if box_o == 1:
                #   print(w_o, wn[box_o])

        # print(o_)
        #print(aver_probs, rank)
        #sys.exit()
        #counts = comm.reduce(counts, op=MPI.SUM, root=0)
        #s = comm.reduce(s, op=MPI.SUM, root=0)
        tmp_aver_probs = comm.gather(aver_probs, root=0)
        tmp_counts = comm.gather(counts, root=0)
        tmp_s = comm.gather(s, root=0)
        tmp_list_probs = comm.gather(list_probs, root =0)
        tmp_number_of_contacts = comm.gather(number_of_contacts, root =0)

        if failed_to_grow <  sweep_length:
          logging.info("fraction of generated confs for 0th box: %f, absolute number %i, rank %i" % (1.0*zero_box/(sweep_length-failed_to_grow), zero_box, rank))
        else:
          logging.info('all confs failed to grow at rank %i' %rank)
        if rank ==0:
          counts = np.mean(np.array(tmp_counts), axis=0)
          #counts[indexes_] = counts[indexes_]- np.min(counts[indexes_])
          s = np.mean(np.array(tmp_s), axis=0)
          #s[indexes_] = s[indexes_]- np.min(s[indexes_])
          #s = s - np.min(s)
          t = np.array(tmp_aver_probs)
          t1 = np.sum(t, axis=0, dtype=float)
          t2 = np.sum(t!=0, axis=0, dtype=float)
          aver_probs = np.divide(t1, t2, out=np.zeros_like(t1), where=t2!=0)
          #logging.info("aver_probs is: %s "%str(t))
          #logging.info("shape of aver_probs is: %s "%str(aver_probs.shape))
          #logging.info(tmp_counts)
           
          for d in tmp_number_of_contacts:
            for k in d.keys():
              if k  in  number_of_contacts_global:
                number_of_contacts_global[k] += d[k]
              else:
                number_of_contacts_global[k] = d[k]
          logging.info("contacts dict for rank=%i: %s"%(rank, number_of_contacts_global))



          for d in tmp_list_probs:
             for k in d.keys():
               list_probs_global[k].append(d[k])  

          for k in list_probs_global.keys():
             if len(list_probs_global[k][0] ) < 100000:  
               list_probs_global[k] = [[item for sublist  in list_probs_global[k] for item in sublist ]]
          logging.info('lengths of list_probs are: %s' %[len(list_probs_global[k][0]) for k in list_probs_global.keys() ])
          logging.info('list_probs are: %s' %[np.median(aux.reject_outliers(np.array(list_probs_global[k][0]))) for k in list_probs_global.keys() ])

          t = counts[indexes_]
          mean = sum(t) / len(t)
        # print('sweep number', sweep_number, 'mean=', round(mean, 2), 'max=', max(t), 'min=', min(t), end='\r')
        # print(t, end='\r')
          logging.info("visits %s" %t)
          logging.info("S %s" %repr(s))
          #logging.info("average probs for states %s" %(repr(aver_probs[0,:])))
        #print(t, s)
        # print(wn)
          if failed_to_grow <  sweep_length:
            logging.info('failed to grow rate: %2.0f%%, out of range in %%  %2.0f%%, and abs number %i  '%(100.0*failed_to_grow/sweep_length, 100.0*out_of_range/(sweep_length-failed_to_grow), out_of_range))
          
          if (max(t) / mean - 1 < flatness) & (1 - min(t) / mean < flatness):
              counts = 0 * counts
              collect_s.append(s.copy())
              ds = ds / decrease
              logging.info("ds=%e, ds_min=%f, sweep number=%i" % (ds, ds_min, sweep_number))

        ds = comm.bcast(ds, root=0)

    return collect_s, sweep_number




def WL_saw_mpi_fast(n, indexes, sweep_length=1000, ds_min=0.0000001, flatness=0.3,
        decrease = 2.0, scale_alpha=2.0, shift_alpha=0.0):
    """
    WL  procedure for calculation of SAW enropies inside matroshkas
    """
    def add_key(k, ncs):
      if k in ncs.keys():
        ncs[k]+=1
      else:
        ncs[k] = 1
      return ncs

    n_boxes = 10#int(n*1.2)
    
    s = np.zeros(n_boxes)
    counts = np.zeros(n_boxes)    
    logging.info('number of contact boxes is: %i' %(n_boxes))
    exclude = []
    if n == 8:
       exclude.append(7)
    indexes_ = [el for el in list(range(0,len(s))) if el not in exclude]
    wn = np.zeros(s.shape)
    nc_bin_o = int((indexes_[-1] + indexes_[0]) /2)
    logging.info('first bin is: %i' %nc_bin_o)
 
    logging.info('s size is: %s' % (s.shape))
    logging.info('counts size is: %s' % (counts.shape))

    logging.info('indexes: %s' % indexes_)

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank ==0:
      number_of_contacts_global = {}
      list_probs_global = {}
      list_number_of_candidates_global = {}
      for i in range(n_boxes):
        list_probs_global[i] = []
        list_number_of_candidates_global[i] = []

    ds = 0.1
    ds0 = ds
    box_o = -1
    w_o = 1
    k_o = 0
    probability_old = 1
    probability_nativei_old = 1
    collect_s = []
    sweep_number = 0
    alpha = 0.0
    frozen =0
    logging.info('scale alpha is: %f'%scale_alpha)
    #logging.info('shift alpha is: %f' %shift_alpha)
    alpha_o = 0
    
    begin = True
    #alphas = np.linspace(scale_alpha, 0, size)
    # nonuniform setting for alpha
    #t  = np.log(np.linspace(1,20 ,size))[::-1]
    #t = t - np.min(t)
    #alphas = t/np.max(t)*scale_alpha
    bin_count = 0  # keeping statistics on last bin  visit  dynamics

    while ds > ds_min:
        sweep_number += 1
        failed_to_grow = 0
        out_of_range = 0
        zero_box = 0
 
        list_probs = {}
        list_number_of_candidates = {}
        for i in range(n_boxes):
          list_probs[i] = []
          list_number_of_candidates[i] = []
 
        number_of_contacts = {}
        collect_contacts = []
        collect_number_of_candidates = []

        for i in range(sweep_length):
            
           #alpha = (np.random.rand()-shift_alpha)*scale_alpha 
            alpha = scale_alpha #np.random.choice((0,scale_alpha))
            alpha = np.random.uniform(-1.8,  .9)
           # else:
            #    alpha = alphas[rank]

            coords_n, w_n, k_n, probability, probability_native  = regrow_saw(n, 0, 0, 0, [], w=1, alpha=alpha, k=0, prob=1, prob_native=0, box = indexes_[-1], indexes=indexes)
            if len(coords_n) < n:
                failed_to_grow += 1
            # print('skipping, failed to grow')
            else:
                coords_n = np.array(coords_n)  # .astype(float)
                
                if math.isnan(probability):
                   logging.error("probability is undefined") 
                
                #u, c = np.unique(coords_n, axis=0, return_counts=True)
                if False:
                #if coords_n.shape != u.shape:
                   print('wrong SAW')
                   #pass
                else:
                  # collecting number of contacts for a configuration                
                  collect_contacts.append(k_n)
                  if k_n in list_probs:
                    list_number_of_candidates[k_n].append(probability_native)
                    list_probs[k_n].append(probability) 


        tmp_collect_contacts = comm.gather(collect_contacts, root=0)
        tmp_list_probs = comm.gather(list_probs, root =0)
        tmp_list_number_of_candidates = comm.gather(list_number_of_candidates, root =0)

        if rank == 0:
           for d in tmp_list_probs:
             for k in d.keys():
               list_probs_global[k].append(d[k])  
           for d in tmp_list_number_of_candidates:
             for k in d.keys():
               list_number_of_candidates_global[k].append(d[k])  


           for k in list_probs_global.keys():
             if len(list_probs_global[k][0] ) < 100000:  
               list_probs_global[k] = [[item for sublist  in list_probs_global[k] for item in sublist ]]
              # if (len(list_probs_global[k][0]) > 0) and (np.median(aux.reject_outliers(np.array(list_probs_global[k][0]))) == np.nan):
              #    logging.error('weird probs detected %s'%(list_probs_global[k][0]))
               #logging.debug('keys %s lens %i median %e entries %s' %(k, len(list_probs_global[k][0]), np.median(aux.reject_outliers(np.array(list_probs_global[k][0]))), list_probs_global[k][0]) )

               list_number_of_candidates_global[k] = [[item for sublist  in list_number_of_candidates_global[k] for item in sublist ]]
           
           logging.info('lengths of list_probs are: %s' %[len(list_probs_global[k][0]) for k in list_probs_global.keys() ])
           logging.info('means of number of candidates are: %s' %[np.mean(list_number_of_candidates_global[k][0])/n for k in list_number_of_candidates_global.keys() ])
           logging.info('list_probs are: %s' %[np.mean(aux.reject_outliers(np.array(list_probs_global[k][0]))) for k in list_probs_global.keys() ])
           #logging.info('list_probs are: %s' %[list_probs_global[k][0] for k in list_probs_global.keys() ])

           global_contacts = [el for  sub in tmp_collect_contacts for el in sub]
           logging.info("the length of global contact list is: %i" %len(global_contacts))

           for nc_bin_n in global_contacts:
    
                  if (nc_bin_n  > indexes_[-1]):
                      nc_bin_n  = indexes_[-1]
                      out_of_range += 1
                  elif nc_bin_n < indexes_[0]:
                    nc_bin_n = indexes_[0]
                  elif nc_bin_n == 0:
                      pass#zero_box +=1
                      #nc = aux.calculate_number_of_contacts(coords_n)
                      #logging.info('spacial box number is  %i for contact bin %i '%(box_n, nc_bin_n)) 
                      #logging.info('configuration probability for 0th box for walker %i is: %e, aver_dist is: %f' % (rank, probability, aver_dist))
                  #elif box_n == indexes_[-2]:
                  #    pass # print(box_n, extreme_x, extreme_y, extreme_z)
                  elif nc_bin_n == indexes_[-2]:
                      #logging.info('contact bin %i,  probs  are %s '%(nc_bin_n, list_probs_global[nc_bin_n][0])) #if np.random.random() < .01:
                      pass #rint(box_n, extreme_x, extreme_y, extreme_z)
                  elif nc_bin_n == int(indexes_[-2] + indexes_[0])/2:
                      pass
                      #logging.info('contact bin %i, prob %s '%(nc_bin_n,  list_probs_global[nc_bin_n][0])) #if np.random.random() < .01:
                  

                  rnd = np.random.random()
                  if (s[nc_bin_o] > s[nc_bin_n]) or (rnd < np.exp(s[nc_bin_o] - s[nc_bin_n])):
                      if s[nc_bin_o] - s[nc_bin_n] >  5:
                        logging.info('rare bin visited rnd=%f, from %5.2f to %5.2f, from box %i to box %i'%(rnd, s[nc_bin_o] , s[nc_bin_n], nc_bin_o,     nc_bin_n))
                      if s[nc_bin_o] - s[nc_bin_n] < - 5:
                        logging.info('low prob event rnd=%f, from %5.2f to %5.2f, from box %i to box %i'%(rnd, s[nc_bin_o] , s[nc_bin_n], nc_bin_o, nc_bin_n))
                      #if nc_bin_n ==  indexes_[-1]:
                      #   logging.debug("entered the last bin %f with %f from bin %f having %f"%(nc_bin_n, s[nc_bin_n], nc_bin_o, s[nc_bin_o]))
                      #   bin_count  = counts[nc_bin_n]
                      #if nc_bin_o == indexes_[-1]:
                      #   bin_count = counts[nc_bin_o] - bin_count                         
                      #   logging.debug("left the last bin %f with %f to  bin %f having %f after %i steps"%(nc_bin_o, s[nc_bin_o], nc_bin_n, s[nc_bin_n], bin_count))
                         
                      alpha_o = alpha
                      nc_bin_o = nc_bin_n
                      probability_old = probability
                      probability_native_old = probability_native

                      frozen =0
                      w_o=w_n
                      k_o=k_n
                      
                  else:
                    pass
                  counts[nc_bin_o] += 1
                  s[nc_bin_o] += ds
                  
                  #number_of_contacts = add_key(nc, number_of_contacts)
                  #list_probs = add_key(probability, list_probs)
                  #list_probs[nc_bin_n].append(probability) 

           t = counts[indexes_]
           mean = sum(t) / len(t)
        # print('sweep number', sweep_number, 'mean=', round(mean, 2), 'max=', max(t), 'min=', min(t), end='\r')
        # print(t, end='\r')
           logging.info("visits %s" %t)
           logging.info("S %s" %repr(s))
          #logging.info("average probs for states %s" %(repr(aver_probs[0,:])))
        #print(t, s)
        # print(wn)
           if failed_to_grow <  sweep_length:
            logging.info('failed to grow rate: %2.0f%%, out of range in %%  %2.0f%%, and abs number %i  '%(100.0*failed_to_grow/sweep_length, 100.0*out_of_range/(sweep_length-failed_to_grow), out_of_range))
          
           if (max(t) / mean - 1 < flatness) & (1 - min(t) / mean < flatness):
              counts = 0 * counts
              collect_s.append(s.copy())
              ds = ds / decrease
              logging.info("ds=%e, ds_min=%f, sweep number=%i" % (ds, ds_min, sweep_number))
              #logging.info('list_probs are: %s' %[list_probs_global[k][0] for k in list_probs_global.keys() ])

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
