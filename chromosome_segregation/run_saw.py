import  numpy as np

import  os
import pandas as pd

LOG_FILE="saw.log"

try:
    from  chromosome_segregation.simulation import URW_saw,WL_saw
    from chromosome_segregation.aux import cache_n_conf, get_grow_caches

    from chromosome_segregation  import overlaps

    from  chromosome_segregation import  consts
    from chromosome_segregation import aux
except:
    from simulation import URW_saw, WL_saw
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
    arr = np.zeros(60)
    for item in Counter(boxes).items():
#     print(item)
        arr[item[0]] = item[1]
    return arr/np.sum(arr)


def f(x, a,b,c):
    return a*x*np.log(x) +  b*x +c


def plot(total_results, thicknesses_x, thicknesses_y, n_boxes, density):
    plt.figure(figsize=(12, 8))
    # thickness=3


    for result, thickness_x, thickness_y in zip(total_results, thicknesses_x, thicknesses_y):
        x = []
        y = []
        boxes = [[i, thickness_x, thickness_y] for i in range(n_boxes, 0, -1)]
        for i in range(len(boxes)):
            n = get_n(boxes[i], density)
            if n>0:
                x.append(1 - n / get_n(boxes[0], density))
            #         x.append(n)
                specific_free_energy = f(1 / n, a=0.375, b=.1347, c=-.2459)
                n_conformations_total = overlaps.Overlap(n).n_conformations
                saw_fraction = np.exp(specific_free_energy * n)
                n_saws = saw_fraction * n_conformations_total

                print(i, boxes[i], n, result[i][boxes[i][0]], result[i][boxes[i][0]] * n_saws)  # , list_to_arr(results[i]))
                y.append(-np.log(result[i][boxes[i][0]] * n_saws) / (n))

        plt.plot(x, y, linestyle='--', marker='o', markersize=10, label='thickness=%i,%i' % (thickness_x, thickness_y))


        # thickness -= 1
    plt.grid()
    plt.xlabel("$1-N/N_0$")
    plt.ylabel("$ F/N$")
    plt.title('density %3.1f; For SAW' % (density))
    plt.legend()
    plt.savefig('thickness.png')


def load_and_process_results():

    data_per_density = [os.path.join(RESULTS_FOLDER, el) for el in os.listdir(RESULTS_FOLDER) if el.startswith('density')]
    logging.info('results  folder per density: %s'%data_per_density)

    density_data = {}
    for d in data_per_density:
        tmp = []
        for root, dirs, files in os.walk(d):
            if files[0].startswith('data'):
                tmp.append(os.path.join(root, files[0]))
        density_data[d.split('/')[-1]] = tmp
        logging.info('density and corresponding data files: %s' %density_data[d.split('/')[-1]])


    for k, files in density_data.items():
        logging.info('working with density: %s' %k)
        dfs = []
        for file in files:
            logging.info('reading file: %s'% file)
            dfs.append(pd.read_csv(file, sep=' ', header=None))
        df = pd.concat(dfs)
        df.columns = ['x', 'y', 'z', 'en', 'n']
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        #     print(df.head(n=50))
        #     aggregated = df.groupby(['x','y','z']).agg({'en':['mean', 'std'], 'n':['mean']})
        aggregated = df.groupby(['x', 'y', 'z']).agg(en_mean=('en', 'mean'), en_std=('en', 'std'), n=('n', 'mean'))

        aggregated.reset_index(inplace=True)
        aggregated['area'] = aggregated['y'] * aggregated['z']

        make_energy_plot(aggregated, save_to=os.path.join(RESULTS_FOLDER, k, 'f_n.png'))
        make_energy_plot_length(aggregated, save_to=os.path.join(RESULTS_FOLDER, k, 'f_length_s.png'))
        make_energy_plot_area(aggregated, save_to=os.path.join(RESULTS_FOLDER, k, 'f_s_length.png'))



def make_energy_plot_area(data, save_to):
    plt.rc('font', size=16)

    plt.figure(figsize=(14, 8))
    for k, group in data.groupby('x'):
        print(k)
        plt.scatter(group['area'], group['en_mean'], label='length=' + str(k))
        plt.errorbar(group['area'], group['en_mean'], yerr=group['en_std'], capsize=10)
    plt.grid()
    plt.legend()
    plt.title('$F/N$ as a function of $S$ for different lengths; contentration=0.5')
    plt.xlabel('area, $S$')
    plt.ylabel('$F/N$')

    plt.savefig(save_to)



def make_energy_plot_length(data, save_to):
    plt.rc('font', size=16)

    plt.figure(figsize=(14, 8))
    for k, group in data.groupby('area'):
        # print(k)
        plt.plot(group['x'], group['en_mean'], marker='p', label='area=' + str(k))
    # plt.errorbar(group['area'], group['en_mean'], yerr=group['en_std'],capsize=10)
    plt.grid()
    plt.legend()
    plt.title('$F/N$ as a function of $x$ for different S; concentration=0.5')
    plt.xlabel('length, $x$')
    plt.ylabel('$F/N$')

    plt.savefig(save_to)



def make_energy_plot(data, save_to):
    plt.rc('font', size=16)
    plt.figure(figsize=(16, 10))

    for k, group in data.groupby('area'):
        plt.plot(1 / group['n'], group['en_mean'], linestyle='-', marker='p', markersize=10, label='$S=$' + str(int(k)))
        plt.errorbar(1/group['n'], group['en_mean'], yerr=group['en_std'], capsize=4)

    plt.grid()
    plt.legend()
    plt.xlabel('$1/N$')
    plt.ylabel('$F/N$')
    plt.title('$F/N$ for different areas; concentration=0.5')

    plt.savefig(save_to)


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
    os.rename(LOG_FILE, os.path.join(OUT_FOLDER, LOG_FILE) )



def process_result(distribution, box, start_from, density):
    """
    given normalized distribution and specific box calculates free energy for the box
    :return:
    :rtype:
    """

    n = get_n(box, density)
    specific_excess_entropy = f(1 / n, a=0.375, b=.1347, c=-.2459)
    logging.info('specific_excess_entropy from the analytical curve  for n=%i is: %4.3f' %(n, specific_excess_entropy))
    n_conformations_total = overlaps.Overlap(n).n_conformations
    logging.info('the total number of phantom conformation for n=%i is %e' %(n, n_conformations_total))
    saw_fraction = np.exp(specific_excess_entropy * n)
    n_saws = saw_fraction * n_conformations_total
    logging.info('number of SAWs for n=%i is %e' %(n, n_saws))

    logging.info('number of SAWs for n=%i inside the box=%s is %e'  %(n, box, np.sum(distribution[:start_from+1]) * n_saws))
    specific_free_energy_for_box = -np.log(np.sum(distribution[:start_from+1]) * n_saws) /n

    logging.info('specific free energy for n=%i and box=%s is: %5.3f' %(n, box, specific_free_energy_for_box))
    return specific_free_energy_for_box




def run(density, n_boxes, thicknesses_x, thicknesses_y):

        logging.info("calculating maximal number of monomers...")
        max_n = get_n([n_boxes[0]+1, thicknesses_x[-1], thicknesses_y[-1]], density)
        logging.info('max_n=%i' % (max_n))


        #consts.caches = cache_n_conf(N_=max_n + 1, dx=25, dy=25, dz=25)
        logging.info('getting cached OR calculating from the scratch..')
        consts.caches = get_grow_caches(fname = GROW_CACHES_FOLDER ,
                params=(max_n+1, 25,25,25))
        logging.info('done calculating (n,dx,dy,dz) array')
        total_results = []
        total_results1 = []

        for thickness_x, thickness_y in zip(thicknesses_x, thicknesses_y):
            #boxes = [[i, thickness_x, thickness_y] for i in range(n_boxes[0],n_boxes[1]-1, -1)]
            boxes = [[i, thickness_x, thickness_y] for i in range(n_boxes, n_boxes-1, -1)]
            logging.info('boxes: %s for thickness_x=%i, thickness_y=%i' % (boxes, thickness_x, thickness_y))

            nsteps = np.linspace(1000000 * max(thickness_x, thickness_y), 9000000, n_boxes[0]-n_boxes[1]+1)
            nsteps = [int(el - el % 1000) for el in nsteps]
            nsteps = nsteps[::-1]
            print('number of steps: %s' % nsteps)

            results = []
<<<<<<< HEAD
            for i in range(n_boxes[0]-n_boxes[1]+1):
=======
            for i in range(len(boxes)):
>>>>>>> ec59c3d4a4b6a17641b19d8327e40300dd92d546
                n = get_n(boxes[i],density)
                # print(
                #     "density: %f, box: %s, #monomers: %i, #steps: %i" % (density, boxes[i], n, nsteps[i]))

                logging.info('running WL_saw with n=%i, nsteps=%i, box=%s' %(n,nsteps[i], boxes[i]))
                if n>0:
                    # all_boxes = URW_saw(n, nsteps[i], box=boxes[i])

                    extend_to_left = 1
                    indexes = aux.get_indexes(boxes[i], extend_to_left=extend_to_left, extend_to_right=3, length=30)
                    logging.info(indexes)
                    s, sweep_number =  WL_saw(n, indexes, sweep_length=10000, ds_min=0.000001, flatness=0.3)
                    logging.info('s: %s', s)
                    # logging.info('making distribution for box %s' %boxes[i])
                    # logging.info(list_to_arr(all_boxes))

                    specific_free_energy = process_result(distribution=s[-1]/sum(s[-1]), box = boxes[i], start_from=extend_to_left, density=density)
                    # specific_free_energy = process_result(distribution=list_to_arr(all_boxes), box=boxes[i], density=density)

                    total_results1.append( (*boxes[i], specific_free_energy, n)  )
                    # results.append(list_to_arr(all_boxes))

                    # total_results.append(results)

        save_results(total_results1, density)


        plot(total_results, thicknesses_x, thicknesses_y, n_boxes, density)



if __name__ == "__main__":

    logging.basicConfig(
                level=logging.INFO,
                format="%(asctime)s [%(levelname)s] %(module)s.%(funcName)s:%(lineno)d %(message)s",
                handlers=[
                    logging.FileHandler(LOG_FILE),
                    logging.StreamHandler()
                ]
            )

    density = 0.4
<<<<<<< HEAD
    n_boxes = (27,23)
=======
    n_boxes = 15
>>>>>>> ec59c3d4a4b6a17641b19d8327e40300dd92d546

    thicknesses_x = list(range(3, 4))
    thicknesses_y = [el  for el in thicknesses_x]

    logging.info("running SAWs with the parameters: density=%3.1f, n_boxes=%i, thicknesses=(%s,%s)" %
                 (density, n_boxes[0]-n_boxes[1]+1, thicknesses_x, thicknesses_y))

    run(density, n_boxes, thicknesses_x, thicknesses_y)

