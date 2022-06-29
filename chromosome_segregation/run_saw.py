import  numpy as np


try:
    from  chromosome_segregation.simulation import URW_saw
    from chromosome_segregation.aux import cache_n_conf, get_grow_caches

    from chromosome_segregation  import overlaps

    from  chromosome_segregation import  consts
except:
    from simulation import URW_saw
    from aux import cache_n_conf, get_grow_caches
    import overlaps
    import consts

import  matplotlib.pyplot as plt

from statistics import Counter


import  logging

def get_n(box, density):
    t = density*box[0]*box[1]*box[2]
    return  int(t - t%2)


def list_to_arr(boxes):
    items = Counter(boxes).items()
#     keys = Counter(boxes).keys()
    arr = np.zeros(30)
    for item in Counter(boxes).items():
#     print(item)
        arr[item[0]] = item[1]
    return arr/np.sum(arr)


def f(x, a,b,c):
    return a*x*np.log(x) +  b*x +c


def plot(total_results, thicknesses):
    plt.figure(figsize=(12, 8))
    # thickness=3


    for result, thickness in zip(total_results, thicknesses):
        x = []
        y = []
        boxes = [[i, thickness, thickness] for i in range(n_boxes, 0, -1)]
        for i in range(len(boxes)):
            n = get_n(boxes[i], density)
            x.append(1 - n / get_n(boxes[0], density))
            #         x.append(n)
            specific_free_energy = f(1 / n, a=0.375, b=.1347, c=-.2459)
            n_conformations_total = overlaps.Overlap(n).n_conformations
            saw_fraction = np.exp(specific_free_energy * n)
            n_saws = saw_fraction * n_conformations_total

            print(i, boxes[i], n, result[i][boxes[i][0]], result[i][boxes[i][0]] * n_saws)  # , list_to_arr(results[i]))
            y.append(-np.log(result[i][boxes[i][0]] * n_saws) / (n))

        plt.plot(x, y, linestyle='--', marker='o', markersize=10, label='thickness=%i' % thickness)
        thickness -= 1
    plt.grid()
    plt.xlabel("$1-N/N_0$")
    plt.ylabel("$ F/N$")
    plt.title('density %3.1f; For SAW' % (density))
    plt.legend()
    plt.savefig('thickness.png')


def run(density, n_boxes, thicknesses):

        logging.info("calculating maximal number of monomers...")
        max_n = get_n([n_boxes, thicknesses[-1], thicknesses[-1]], density)
        logging.info('max_n=%i' % (max_n))


        #consts.caches = cache_n_conf(N_=max_n + 1, dx=25, dy=25, dz=25)
        logging.info('getting cached OR calculating from the scratch..')
        consts.caches = get_grow_caches(fname='grow_caches.txt.gz',
                params=(max_n+1, 25,25,25))
        logging.info('done calculating n,dx,dy,dz array')
        total_results = []
        for thickness in thicknesses:
            boxes = [[i, thickness, thickness] for i in range(n_boxes, 0, -1)]
            logging.info('boxes: %s for thickness %i' % (boxes, thickness))

            nsteps = np.linspace(20000 * thickness, 200000, len(boxes))
            nsteps = [int(el - el % 1000) for el in nsteps]
            nsteps = nsteps[::-1]
            print('number of steps: %s' % nsteps)

            results = []
            for i in range(len(boxes)):
                print(
                    "density: %f, box: %s, #monomers: %i, #steps: %i" % (density, boxes[i], get_n(boxes[i],density), nsteps[i]))
                boxes_ = URW_saw(get_n(boxes[i], density), nsteps[i], box=boxes[i])
                results.append(list_to_arr(boxes_))

            total_results.append(results)

        plot(total_results, thicknesses)


if __name__ == "__main__":

    density = 0.5
    n_boxes = 7

    thicknesses = list(range(2, 4))
    logging.info("running SAWs with the parameters: density=%3.1f, n_boxes=%i, thicknesses=%s" %(density, n_boxes, thicknesses))

    run(density, n_boxes, thicknesses)

