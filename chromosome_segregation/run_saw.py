import  numpy as np

import  os

LOG_FILE="saw.log"

try:
    from  chromosome_segregation.simulation import URW_saw
    from chromosome_segregation.aux import cache_n_conf, get_grow_caches

    from chromosome_segregation  import overlaps

    from  chromosome_segregation import  consts
    from chromosome_segregation import aux
except:
    from simulation import URW_saw
    from aux import cache_n_conf, get_grow_caches
    import overlaps
    import consts
    import  aux

import  matplotlib.pyplot as plt

from statistics import Counter


import  logging


RESULTS_FOLDER = os.path.join(os.path.dirname(consts.__file__), consts.RESULTS_SAW)
GROW_CACHES_FOLDER = os.path.join(os.path.dirname(consts.__file__), 'grow_caches.txt.gz')


def get_n(box, density):
    t = density*box[0]*box[1]*box[2]
    return  int(t - t%2)


def list_to_arr(boxes):
    # items = Counter(boxes).items()
#     keys = Counter(boxes).keys()
    arr = np.zeros(30)
    for item in Counter(boxes).items():
#     print(item)
        arr[item[0]] = item[1]
    return arr/np.sum(arr)


def f(x, a,b,c):
    return a*x*np.log(x) +  b*x +c


def plot(total_results, thicknesses,n_boxes,density):
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



def save_results(results, density):
    """

    :param results:
    :type results:
    :return:
    :rtype:
    """

    OUT_FOLDER = os.path.join(RESULTS_FOLDER, 'density_'+str(density))

    if not os.path.exists(OUT_FOLDER):
        os.makedirs(OUT_FOLDER)

    experiment_folder = aux.get_subfolder(OUT_FOLDER)

    OUT_FOLDER = os.path.join(OUT_FOLDER, experiment_folder)

    data = np.array(results)

    np.savetxt(os.path.join(OUT_FOLDER, 'data.csv'),data, fmt='%3.1f %3.1f %3.1f %.4f %i'  )



def process_result(distribution, box, density):
    """
    given normalized distribution and specific box calculates free energy for the box
    :return:
    :rtype:
    """

    n = get_n(box, density)
    specific_free_energy = f(1 / n, a=0.375, b=.1347, c=-.2459)
    n_conformations_total = overlaps.Overlap(n).n_conformations
    saw_fraction = np.exp(specific_free_energy * n)
    n_saws = saw_fraction * n_conformations_total

    return -np.log(distribution[box[0]] * n_saws) / n




def run(density, n_boxes, thicknesses):

        logging.info("calculating maximal number of monomers...")
        max_n = get_n([n_boxes, thicknesses[-1], thicknesses[-1]], density)
        logging.info('max_n=%i' % (max_n))


        #consts.caches = cache_n_conf(N_=max_n + 1, dx=25, dy=25, dz=25)
        logging.info('getting cached OR calculating from the scratch..')
        consts.caches = get_grow_caches(fname = GROW_CACHES_FOLDER ,
                params=(max_n+1, 25,25,25))
        logging.info('done calculating (n,dx,dy,dz) array')
        total_results = []
        total_results1 = []

        for thickness in thicknesses:
            boxes = [[i, thickness, thickness] for i in range(n_boxes, 0, -1)]
            logging.info('boxes: %s for thickness %i' % (boxes, thickness))

            nsteps = np.linspace(20000 * thickness, 200000, n_boxes)
            nsteps = [int(el - el % 1000) for el in nsteps]
            nsteps = nsteps[::-1]
            print('number of steps: %s' % nsteps)

            results = []
            for i in range(n_boxes):
                n = get_n(boxes[i],density)
                print(
                    "density: %f, box: %s, #monomers: %i, #steps: %i" % (density, boxes[i], n, nsteps[i]))
                all_boxes = URW_saw(n, nsteps[i], box=boxes[i])
                logging.info('making distribution of boxes')

                specific_free_energy = process_result(distribution=list_to_arr(all_boxes), box=boxes[i], density=density)
                total_results1.append( (*boxes[i], specific_free_energy, n)  )
                results.append(list_to_arr(all_boxes))

            total_results.append(results)

        save_results(total_results1, density)


        plot(total_results, thicknesses, n_boxes, density)


if __name__ == "__main__":

    logging.basicConfig(
                level=logging.INFO,
                format="%(asctime)s [%(levelname)s] %(module)s.%(funcName)s:%(lineno)d %(message)s",
                handlers=[
                    logging.FileHandler(LOG_FILE),
                    logging.StreamHandler()
                ]
            )

    density = 0.5
    n_boxes = 8

    thicknesses = list(range(2, 6))
    logging.info("running SAWs with the parameters: density=%3.1f, n_boxes=%i, thicknesses=%s" %(density, n_boxes, thicknesses))

    run(density, n_boxes, thicknesses)

