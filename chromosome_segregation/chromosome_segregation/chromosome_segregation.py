"""Main module."""

from  chromosome_segregation.chromosome_segregation.aux import  plot, glue_s, rescale_s, save_results, load_data

from  chromosome_segregation.chromosome_segregation.simulation import WL, URW

import  logging
import  time
import  numpy as np


def run(n, n_steps, bins=None, counts=None):
    """
    A single experiment.
    The experiment includes 3 independent runs. 1) URW 2)WL for left part of distribution 3)WL for right part
    of distribution

    :param n: number of monomers
    :type n: int
    :param n_steps: number of conformation changes
    :type n_steps: int
    :param bins: Numpy array with number of overlaps [0,1,...n_max]
    :type bins: numpy array
    :param counts: visit counts for URW for ``bins``
    :type counts: numpy array
    :return: ``bins``,  ``counts``, ``s_left`` -- left wing of entropy, ``s_right`` -- right wing of entropy, metrics -- simulation parameters
    :rtype: tuple
    """

    metrics = {}

    #########
    # URW
    #########
    metrics['n'] = n
    start_time = time.time()
    logging.info("starting URW for %i beads, %i steps"%(n, n_steps))
    bins, counts = URW(n=n, n_steps=n_steps)
    metrics['URW_run_time'] = time.time() - start_time
    logging.info("Running URW took %4.0f seconds" %metrics['URW_run_time'] )

    # bins_bias_100, counts_bias_100, out = URW_biased(n=n, n_steps=n_steps, alpha=70.)
    logging.info("normalizing to 1 counts...")
    counts = counts / np.nansum(counts)
    # print(counts)
    metrics['n_steps'] = n_steps

    # counts = np.array(counts)
    # print('counts', counts, )

    sweep_length = 1000
    metrics['sweep_length'] = sweep_length

    #########
    # WL  left
    #########
    if n > 30:
        max_overlaps = np.argmax(np.array(counts) > 10 ** (-2))
        alpha = 1.5  # 2.0#1.7
        ds_min = 1 * 10 ** -7  # 0.0000001
        flatness = 0.3  # 0.1
    else:
        max_overlaps = np.argmax(np.array(counts)) + 5
        alpha = 0.0
        ds_min = 1 * 10 ** -8  # 0.0000001
        flatness = 0.05

    logging.info("running WL-left with max_overlaps=%i, ds_min=%f, alpha=%f, flatness=%3.2f"
                 %(int(max_overlaps),ds_min, alpha, flatness))
    # print('max overlap = ', max_overlaps)
    metrics['max_overlaps_left'] = int(max_overlaps)
    metrics['flatness_left'] = flatness

    metrics['alpha_left'] = alpha

    metrics['ds_min_left'] = ds_min
    start_time = time.time()
    s_left, n_sweeps_left = WL(n, max_overlaps, exclude=(), alpha=alpha, sweep_length=sweep_length, ds_min=ds_min,
                               flatness=flatness)

    #     s_left, n_sweeps_left = None, None
    metrics['WL_left_run_time'] = time.time() - start_time
    logging.info("Running WL-left took %4.0f seconds and %i sweeps" %(metrics['WL_left_run_time'], n_sweeps_left ))

    metrics['n_sweeps_left'] = n_sweeps_left

    #########
    # WL  right
    #########

    if n > 30:
        alpha = -0.5
        grain = 5
    elif (n <= 30) and (n > 20):
        alpha = -0.5
        grain = 3
    else:
        alpha = -0.7
        grain = 1

    exclude = ()
    if (n == 12): exclude = (23, 24, 26, 27, 28, 29)
    if (n == 10): exclude = (15, 17, 18, 19)

    metrics['grain_right'] = grain
    metrics['alpha_right'] = alpha

    min_overlaps = np.nanargmax(counts) + np.nanargmax(counts[np.nanargmax(counts):] < 10 ** -2)

    addition = np.nanargmax(counts[np.nanargmax(counts):] <= 1 / n_steps)
    if addition == 0:  # if no such element
        addition = len(counts[np.nanargmax(counts):])

    max_overlaps = np.nanargmax(counts) + addition

    min_overlaps = min_overlaps - min_overlaps % grain
    max_overlaps = max_overlaps - max_overlaps % grain

    print("np.nanargmax(counts)", np.nanargmax(counts), "min_overlaps", min_overlaps, "max_overlaps", max_overlaps)
    if (max_overlaps <= min_overlaps):
        print(30 * 'XX', 'something wrong in borders')

    metrics['min_overlaps_right'] = int(min_overlaps)
    metrics['max_overlaps_right'] = int(max_overlaps)

    ds_min = 1 * 10 ** -6
    metrics['ds_min_right'] = ds_min

    flatness = 0.3
    metrics['flatness_right'] = flatness

    logging.info("running WL-right with min_overlaps=%i, max_overlaps=%i, ds_min=%f, alpha=%f, flattness=%3.2f"
                 % (int(min_overlaps), int(max_overlaps), ds_min, alpha, flatness))

    start_time = time.time()
    s_right, n_sweeps_right = WL(n, max_overlaps, min_overlaps, exclude=exclude, alpha=alpha,
                                 sweep_length=sweep_length, grain=grain, ds_min=ds_min, flatness=flatness)

    #     s_right,n_sweeps_right = None, None
    metrics['WL_right_run_time'] = time.time() - start_time
    metrics['n_sweeps_right'] = n_sweeps_right
    logging.info("Running WL-left took %4.0f seconds and %i sweeps" % (metrics['WL_right_run_time'], n_sweeps_right))

    return bins, counts, s_left, s_right, metrics



def run_all(n, n_steps):
    """
    global run
    """
    logging.info("running 'run'...")
    bins, counts, s_left, s_right, metrics = run(n, n_steps)

    logging.info("gluing left and right entropies with counts")
    total_s = glue_s(bins, counts, s_left[-1], s_right[-1])
    logging.info("saving results")
    experiment_folder = save_results(n, s_left, s_right, total_s, bins, counts, metrics)
    logging.info("loading saved data for plotting and calculating #SAW")
    s_left, s_right, s_total, bins, counts, metrics = load_data(experiment_folder)
    logging.info('plotting')
    plot(bins, counts, total_s, metrics, save_plot_to=experiment_folder)




if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(module)s.%(funcName)s:%(lineno)d %(message)s",
        handlers=[
            logging.FileHandler("debug.log"),
            logging.StreamHandler()
        ]
    )

    ns = [80, 80]
    n_stepss = [200000, 2000000]
    logging.info("looping number of beads")
    for n, n_steps in zip(ns, n_stepss):
        logging.info('number of beads %i, number of steps %i' %(n, n_steps))
        run_all(n, n_steps)