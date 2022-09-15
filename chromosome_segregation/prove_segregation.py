###
# Proving that two mixes ring polymers will segregate
###

LOG_FILE="segregation_prove.log"

import logging
import  os
import numpy as np

try:
    import consts, regrow, aux
except:
    from chromosome_segregation import consts, regrow, aux

RESULTS_FOLDER = os.path.join(os.path.dirname(consts.__file__), consts.RESULTS_SEGREGATION_PROVE)



def save_results(results, area, n):
    """

    :param results:
    :type results:
    :return:
    :rtype:
    """

    OUT_FOLDER = os.path.join(RESULTS_FOLDER, 'area_'+str(area))

    if not os.path.exists(OUT_FOLDER):
        os.makedirs(OUT_FOLDER)

    OUT_FOLDER = os.path.join(OUT_FOLDER, 'n_'+str(n))

    if not os.path.exists(OUT_FOLDER):
        os.makedirs(OUT_FOLDER)

    experiment_folder = aux.get_subfolder(OUT_FOLDER)

    OUT_FOLDER = os.path.join(OUT_FOLDER, experiment_folder)

    data = np.array(results)

    np.savetxt(os.path.join(OUT_FOLDER, 'data.csv'),data, fmt='%3.2f %3.2f'  )




def get_mixing_degree(coords1, coords2):
    """
    calculates the mixing degree between two polymers by calculating the number of the same x-coordinates
    The larger the overlap the larger the mixing degree.

    degree is in [0, 1] interval.


    :param coords1:
    :type coords1:
    :param coords2:
    :type coords2:
    :return:
    :rtype:
    """

    arr1 = np.array(coords1)[:,0]
    arr2 = np.array(coords2)[:,0]

    overlaps1 = np.array([el for el in arr1 if el in arr2]).shape[0]
    overlaps2 = np.array([el for el in arr2 if el in arr1]).shape[0]

    degree =  (overlaps1+overlaps2) / 2.0 / arr1.shape[0]

    if degree > 1:
        print(degree, arr1, arr2, overlaps1, overlaps2, arr1.shape[0] )
        import sys
        sys.exit()

    return degree




def get_cm_distance(coords1, coords2):
    """
    distance between centers of mass of polymers

    :param coords1: coordinates of the first polymer
    :type coords1: list of lists [ [1,2,3], [4,5,6] ]
    :param coords2:
    :type coords2:
    :return:  cms' distance
    :rtype: float
    """

    cm1 = np.average(np.array(coords1), axis=0)
    cm2 = np.average(np.array(coords2), axis=0)
    diff = cm1 - cm2

    return abs(diff[0])# np.dot(diff,diff)


def run_urw(n, iter_max=10000):

    logging.info('getting cached OR calculating from the scratch..')
    consts.caches = aux.get_grow_caches(fname=consts.GROW_CACHES_FOLDER,
                                    params=(n + 1, 25, 25, 25))
    logging.info('done calculating (n,dx,dy,dz) array')

    cm_distribution = []

    max_n_fails = 5
    for i in range(iter_max):
        if i%1000 == 0:
            logging.info("%3.1f %%"%(1.0*i/iter_max*100))
        coords1 = coords2 = []
        n_fails = 0
        while (len(coords1) < n) and (n_fails < max_n_fails):
            coords1, w1, k = regrow.regrow_saw_segregation_prove(n, 0, 0, 0, [], w=1, alpha=0.0, k=0, coords=[])
            n_fails +=1
        # print('coords1', w1, k, coords1, n_fails)
        n_fails = 0
        while (len(coords2) < n) and (n_fails < max_n_fails):
            coords2, w2, k = regrow.regrow_saw_segregation_prove(n, 0, 0, 0, [], w=1, alpha=0.0, k=0, coords=coords1)
            n_fails +=1
        # print('coords2', w2, k, coords2)
        if (n_fails < max_n_fails) and (len(coords2) == n) and (len(coords1) == n):
            cm_distance = get_cm_distance(coords1, coords2)
            mixing_degree = get_mixing_degree(coords1, coords2)


            # logging.info("cm distance is: %3.1f" %cm_distance)
            cm_distribution.append((cm_distance, mixing_degree))
        # count = consts.caches[n-2, 1, 0, 0]
        # count = aux.n_conf(n,0,0,0)
        # print(1/count, consts.caches[n-1,0,0,0],  consts.caches[n-2,1,0,0],consts.caches[n-3,0,0,0])

        # tmp1 = set(tuple(el) for el in coords1)
        # tmp2 = set(tuple(el) for el in coords2)

    return cm_distribution #tmp1.intersection( tmp2)



def run_urw_simultaneously(n, iter_max=10000):

    logging.info('getting cached OR calculating from the scratch..')
    consts.caches = aux.get_grow_caches(fname=consts.GROW_CACHES_FOLDER,
                                    params=(n + 1, 25, 25, 25))
    logging.info('done calculating (n,dx,dy,dz) array')

    cm_distribution = []

    max_n_fails = 5
    for i in range(iter_max):
        if i%1000 == 0:
            logging.info("%3.1f %%"%(1.0*i/iter_max*100))
        coords1 = coords2 = []
        n_fails = 0
        while (len(coords1) < n) and (n_fails < max_n_fails):
            coords1, coords2 = regrow.regrow_saw_segregation_prove_two_chains(n,0,0,0, n, 0, 0, 0, res1=[], res2=[])
            # coords1, w1, k = regrow.regrow_saw_segregation_prove(n, 0, 0, 0, [], w=1, alpha=0.0, k=0, coords=[])
            n_fails +=1
        logging.debug("coords1 length is: %i; coords2 length is: %i "%(len(coords1), len(coords2)))
        if (n_fails < max_n_fails) and (len(coords2) == n) and (len(coords1) == n):
            cm_distance = get_cm_distance(coords1, coords2)
            mixing_degree = get_mixing_degree(coords1, coords2)


            logging.debug("cm distance is: %3.1f" %cm_distance)
            cm_distribution.append((cm_distance, mixing_degree))
        # count = consts.caches[n-2, 1, 0, 0]
        # count = aux.n_conf(n,0,0,0)
        # print(1/count, consts.caches[n-1,0,0,0],  consts.caches[n-2,1,0,0],consts.caches[n-3,0,0,0])

        # tmp1 = set(tuple(el) for el in coords1)
        # tmp2 = set(tuple(el) for el in coords2)
        # print(repr(coords1))
        # print(repr(coords2))
    return cm_distribution #tmp1.intersection( tmp2)


def run_wl(n, ds_min=0.01):

    logging.info("running WL")

    logging.info('getting cached OR calculating from the scratch..')
    consts.caches = aux.get_grow_caches(fname=consts.GROW_CACHES_FOLDER,
                                        params=(n + 1, 25, 25, 25))
    logging.info('done calculating (n,dx,dy,dz) array')

    scale = 4
    max_dist = 5*scale
    s = np.zeros(max_dist)
    visits = np.zeros(s.shape)


    ds = 1
    max_n_fails = 5
    o_ = 1
    sweep_number = 0
    flatness = 0.3

    while ds > ds_min:
        sweep_number += 1

        for i in range(1000):
            coords1 = coords2 = []
            n_fails = 0
            while (len(coords1) < n) and (n_fails < max_n_fails):
                coords1, w1, k = regrow.regrow_saw_segregation_prove(n, 0, 0, 0, [], w=1, alpha=0.0, k=0, coords=[])
                n_fails += 1
            # print('coords1', w1, k, coords1, n_fails)
            n_fails = 0
            while (len(coords2) < n) and (n_fails < max_n_fails):
                coords2, w2, k = regrow.regrow_saw_segregation_prove(n, 0, 0, 0, [], w=1, alpha=0.0, k=0,
                                                                     coords=coords1)
                n_fails += 1
            if (n_fails < max_n_fails) and (len(coords2) == n) and (len(coords1) == n):
                n_ = int(np.round(get_cm_distance(coords1, coords2) * scale))
                if n_> max_dist: logging.info("cm distance is: %3.1f" % n_)
                if  n_ < max_dist:
                    if (np.random.random() <  np.exp(s[o_] - s[n_])) :
                        o_ = n_

            visits[o_] += 1
            s[o_] += ds

        mean = np.mean(visits) #(t) / len(t)
        print('sweep number', sweep_number, 'mean=', round(mean, 2), 'max=', max(visits), 'min=', min(visits), end='\r')
        # print(t, end='\r')

        if (max(visits) / mean - 1 < flatness) & (1 - min(visits) / mean < flatness):
            visits[:] = 0
            # counts = 0 * counts
            # collect_s.append(s.copy())
            ds = ds / 1.5
            logging.info("ds=%e, ds_min=%f, sweep number=%i" % (ds, ds_min, sweep_number))

    return s



if __name__ == "__main__":

    logging.basicConfig(
                level=logging.INFO,
                format="%(asctime)s [%(levelname)s] %(module)s.%(funcName)s:%(lineno)d %(message)s",
                handlers=[
                    logging.FileHandler(LOG_FILE),
                    logging.StreamHandler()
                ]
            )

    logging.info("starting segregation prove routine")

    # print(run_urw(n=10))
    # run_wl(n=10)
    ns = [10]
    for n in ns:
        for i in range(3):
            iter_max = 1000 * n
            logging.info("running URW for n=%i, number of iterations = %i" %(n, iter_max))
            # results = run_urw(n=n, iter_max=iter_max)
            results = run_urw_simultaneously(n=n, iter_max=iter_max)

            area = consts.size_z  * consts.size_y
            save_results(results, area, n)



