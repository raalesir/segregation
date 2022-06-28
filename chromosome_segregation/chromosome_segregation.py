"""Main module."""

try:
    from  aux import  plot, glue_s, rescale_s, save_results, load_data,\
        cache_n_conf, calculate_saw_fraction_all, prepare_entropy_plot, get_result_subfolders,plot_specific_entropy,\
        calculate_saw_fraction,  plot_specific_entropy_comparison

    from  simulation import WL, URW
    import consts

except:
    from chromosome_segregation.aux import   plot, glue_s, rescale_s, save_results, load_data,\
        cache_n_conf, calculate_saw_fraction_all, prepare_entropy_plot, get_result_subfolders, plot_specific_entropy,\
        calculate_saw_fraction, plot_specific_entropy_comparison

    from chromosome_segregation.simulation import WL, URW




import  logging
import  time
import  numpy as np
import shutil, os
import argparse
from argparse import RawTextHelpFormatter



LOG_FILE="log.log"

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
    logging.info("To run URW took %4.0f seconds" %metrics['URW_run_time'] )

    # bins_bias_100, counts_bias_100, out = URW_biased(n=n, n_steps=n_steps, alpha=70.)
    logging.info("normalizing overlap counts to 1 ...")
    counts = counts / np.nansum(counts)
    # print(counts)
    metrics['n_steps'] = n_steps

    # counts = np.array(counts)
    # print('counts', counts, )

    sweep_length = 1000
    metrics['sweep_length'] = sweep_length



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
        alpha = -0.8
        grain = 1

    exclude = ()
    if (n == 12): exclude = (23, 24, 26, 27, 28, 29)
    if (n == 10): exclude = (15, 17, 18, 19)
    if  (n==8):  exclude =(10,11)

    metrics['grain_right'] = grain
    metrics['alpha_right'] = alpha

    min_overlaps = np.nanargmax(counts) + np.nanargmax(counts[np.nanargmax(counts):] < 10 ** -2)

    addition = np.nanargmax(counts[np.nanargmax(counts):] <= 1 / n_steps)
    if addition == 0:  # if no such element
        addition = len(counts[np.nanargmax(counts):])

    max_overlaps = np.nanargmax(counts) + addition

    min_overlaps = min_overlaps - min_overlaps % grain
    max_overlaps = max_overlaps - max_overlaps % grain

    logging.info("np.nanargmax(counts)=%i, min_overlaps=%i; max_overlaps=%i"%(np.nanargmax(counts),min_overlaps, max_overlaps))
    if (max_overlaps <= min_overlaps):
        print(30 * 'XX', 'something wrong in borders')

    metrics['min_overlaps_right'] = int(min_overlaps)
    metrics['max_overlaps_right'] = int(max_overlaps)

    ds_min = 10 **-4
    metrics['ds_min_right'] = ds_min

    flatness = 0.25
    metrics['flatness_right'] = flatness

    logging.info("running WL-right with min_overlaps=%i, max_overlaps=%i, ds_min=%f, alpha=%f, flattness=%3.2f"
                 % (int(min_overlaps), int(max_overlaps), ds_min, alpha, flatness))

    start_time = time.time()
    s_right, n_sweeps_right = WL(n, max_overlaps, min_overlaps, exclude=exclude, alpha=alpha,
                                 sweep_length=sweep_length, grain=grain, ds_min=ds_min, flatness=flatness)

    metrics['WL_right_run_time'] = time.time() - start_time
    metrics['n_sweeps_right'] = n_sweeps_right
    logging.info("Running WL-right took %4.0f seconds and %i sweeps" % (metrics['WL_right_run_time'], n_sweeps_right))


    #########
    # WL  left
    #########
    if n > 30:
        max_overlaps = np.argmax(np.array(counts) > 10 ** (-2))
        alpha = 1.5  # 2.0#1.7
        ds_min = 10**-5#10 ** -7  # 0.0000001
        flatness = 0.25  # 0.1
    else:
        max_overlaps = np.argmax(np.array(counts)) + 5
        alpha = 0.0
        ds_min = 10 ** -4  # 1 * 10 ** -8  # 0.0000001
        flatness = 0.05

    logging.info("running WL-left with max_overlaps=%i, ds_min=%e, alpha=%f, flatness=%3.2f"
                 % (int(max_overlaps), ds_min, alpha, flatness))
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
    logging.info("Running WL-left took %4.0f seconds and %i sweeps" % (metrics['WL_left_run_time'], n_sweeps_left))

    metrics['n_sweeps_left'] = n_sweeps_left

    return bins, counts, s_left, s_right, metrics



def run_all(n, n_steps):
    """
    global run
    """
    logging.info("running 'run'...")
    bins, counts, s_left, s_right, metrics = run(n, n_steps)

    logging.info("gluing left and right entropy with counts")
    total_s = glue_s(bins, counts, s_left[-1], s_right[-1])

    logging.info("calculating SAW fraction")
    reverse_n, saw_fraction, fitted_s = calculate_saw_fraction(n, total_s)


    logging.info("saving results")
    experiment_folder = save_results(n, s_left, s_right, total_s, bins, counts, metrics, saw_fraction, fitted_s)
    logging.info("loading saved data for plotting and calculating #SAW")
    s_left, s_right, s_total, bins, counts, metrics = load_data(experiment_folder)
    logging.info('plotting')
    plot(bins, counts, total_s, metrics, save_plot_to=experiment_folder)
    # reverse_n_, specific_excess_entropy_, entropy = calculate_saw_fraction_all(experiment_folder)

    return  experiment_folder




if __name__ == "__main__":

    my_parser = argparse.ArgumentParser(description="""Calculating  the specific excess entropy for a ring grid
polymer as a  function of a reciprocal number  of beads""", formatter_class=RawTextHelpFormatter)

    # Add the arguments
    # my_parser.add_argument('Path',
    #                        metavar='path',
    #                        type=str,
    #                        help='the path to list')
    #
    # # Execute the parse_args() method
    # args = my_parser.parse_args()
    #
    # input_path = args.Path


    # logging.basicConfig(
    #     level=logging.INFO,
    #     format="%(asctime)s [%(levelname)s] %(module)s.%(funcName)s:%(lineno)d %(message)s",
    #     handlers=[
    #         logging.FileHandler(LOG_FILE),
    #         logging.StreamHandler()
    #     ]
    # )

    ns = [10,10]
    n_stepss = [50000, 50000]
    print("caching the number of confs... can take several minutes")

    consts.caches = cache_n_conf(N_=max(ns), dx=30, dy=30 , dz=30)
    print("done caching. The cache shape is %s" % str(consts.caches.shape))

    # print("looping number of beads")
    # for n, n_steps in zip(ns, n_stepss):
    #
    #     logging.basicConfig(
    #         level=logging.INFO,
    #         format="%(asctime)s [%(levelname)s] %(module)s.%(funcName)s:%(lineno)d %(message)s",
    #         handlers=[
    #             logging.FileHandler(LOG_FILE),
    #             logging.StreamHandler()
    #         ]
    #     )
    #
    #
    #     logging.info("logging to: %s"%(LOG_FILE))
    #     logging.info('running simulation for number of beads %i, number of steps %i' %(n, n_steps))
    #     experiment_folder = run_all(n, n_steps)
    #     shutil.move(LOG_FILE, os.path.join(experiment_folder, LOG_FILE))
    #     logging.shutdown()


    subfolders = get_result_subfolders(path=consts.RESULTS_FOLDER_FREE)
    print(subfolders)
    x, y, errs = prepare_entropy_plot(subfolders)
    # print(x,y,errs)
    plot_specific_entropy(x,y,errs=errs, save_to=consts.RESULTS_FOLDER_FREE)

    subfolders = get_result_subfolders(path=consts.RESULTS_FOLDER_HALFSPACE)
    print(subfolders)
    x_half, y_half, errs_half = prepare_entropy_plot(subfolders)
    # print(x,y,errs)
    plot_specific_entropy(x_half,y_half,errs=errs_half, save_to=consts.RESULTS_FOLDER_HALFSPACE)

    subfolders = get_result_subfolders(path=consts.RESULTS_FOLDER_BOX)
    print(subfolders)
    x_box, y_box, errs_box = prepare_entropy_plot(subfolders)
    # print(x,y,errs)
    plot_specific_entropy(x_box, y_box, errs=errs_box, save_to=consts.RESULTS_FOLDER_BOX)


    plot_specific_entropy_comparison(x, y, errs, x_half, y_half, errs_half, fname='specific_entropy_free_half_comparison.png')
    plot_specific_entropy_comparison(x, y, errs, x_box, y_box, errs_box,fname='specific_entropy_free_box_comparison.png')




