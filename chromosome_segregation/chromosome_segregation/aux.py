"""
Helpers
"""
import  os
import  json
import  numpy as np

import  math
import  logging
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

RESULTS_FOLDER = 'results'


def get_subfolder(path):
    """
    Given path returns a folder for the next experiment.
    Looks into ``path``, checks for the existing  subfolders and creates a new one with incremented by one
    number. Example: `run_6` if `run_0.... run_5` already present

    :param path: path
    :type path: str
    :return: folder name
    :rtype: str
    """


    experiment_folder = None
    for root, dirs, files in os.walk(path):
        dirs = sorted([el for el in dirs if el.startswith('run')], key=lambda x: int(x.split('_')[1]))
        if len(dirs) == 0:
            experiment_folder = 'run_0'
            os.makedirs(os.path.join(path, experiment_folder))
        else:
            number = int(dirs[-1].split('_')[-1]) + 1
            experiment_folder = 'run_' + str(number)

            os.makedirs(os.path.join(path, experiment_folder))

        return experiment_folder



def save_results(n, s_left, s_right, s_total, bins, counts, metrics):
    """
    Given the number ``n`` -- the number of monomers  -- saves the data to ``results`` folder

    :param n:
    :type n:
    :param s_left:
    :type s_left:
    :param s_right:
    :type s_right:
    :param s_total:
    :type s_total:
    :param bins:
    :type bins:
    :param counts:
    :type counts:
    :param metrics:
    :type metrics:
    :return:
    :rtype:
    """

    OUT_FOLDER = os.path.join(RESULTS_FOLDER, str(n))

    if not os.path.exists(OUT_FOLDER):
        os.makedirs(OUT_FOLDER)

    experiment_folder = get_subfolder(OUT_FOLDER)

    OUT_FOLDER = os.path.join(OUT_FOLDER, experiment_folder)

    names = ['s_left.txt', 's_right.txt', 's_total.txt', 'bins.txt', 'counts.txt']
    data = [s_left, s_right, s_total, bins, counts]
    save_dict = dict(zip(names, data))

    for item in save_dict:
        f_name = os.path.join(OUT_FOLDER, item)
        np.savetxt(f_name, save_dict[item])

    with open(os.path.join(OUT_FOLDER, 'metrics.txt'), 'w') as f:
        f.write(json.dumps(metrics))

    return OUT_FOLDER



def load_data(experiment_folder):
    names = ['s_left.txt', 's_right.txt', 's_total.txt', 'bins.txt', 'counts.txt', 'metrics.txt']

    FOLDER = experiment_folder

    bins = np.loadtxt(os.path.join(FOLDER, names[3]))
    counts = np.loadtxt(os.path.join(FOLDER, names[4]))
    s_left = np.loadtxt(os.path.join(FOLDER, names[0]))
    s_right = np.loadtxt(os.path.join(FOLDER, names[1]))
    s_total = np.loadtxt(os.path.join(FOLDER, names[2]))

    with open(os.path.join(FOLDER, 'metrics.txt'), 'r') as f:
        metrics = json.load(f)

    return s_left, s_right, s_total, bins, counts, metrics



def plot(bins, counts, collect_s, metrics, save_plot_to):
    nsteps = metrics['n_steps']
    nsweeps = metrics['n_sweeps_left'] + metrics['n_sweeps_right']
    sw_len = metrics['sweep_length']
    n = metrics['n']

    plt.rc('font', size=20)
    plt.figure(figsize=(16, 10))

    plt.scatter(np.arange(len(collect_s)), collect_s, facecolor='black', marker='.', s=50, label='WL')

    # plt.scatter(bins11_80, counts_11_80,facecolor='blue', alpha=0.5, label='URW-biased', marker='x', s=150)

    plt.scatter(bins, counts, facecolor='none', edgecolor='red', alpha=0.5, label='URW', marker='o', s=150)
    plt.yscale('log')
    plt.legend()
    plt.title('%i confs for n=%i beads; number of sweeps %i of %i length' % (nsteps, n, nsweeps, sw_len))
    plt.grid()
    plt.xlabel('number of overlaps')
    plt.ylabel('fraction of configurations')

    plt.savefig(os.path.join(save_plot_to, 'overlaps.png'), )



def rescale_s(s_, bins, counts, bin_index, s_index):

    # index in s equals to the number of overlaps. I.e. index=k means k-overlaps.
    #     first_non_zero = np.argmax(s>0) #np.where(s >0)[0][0]
    #     print('first_non_zero', first_non_zero)

    # finding the  index in "bins" of the element equal to first non zero in s.
    #     index = np.where(bins == first_non_zero )[0][0] #index in bins for the last_s [-1] number of overlaps
    #     index = np.argmax(bins == first_non_zero )
    #     print('index', index)
    #     print('len(counts)', len(counts))
    #     print('len(s)', len(s_))
    delta = s_[s_index] - np.log(counts[bin_index])  # len(last_s) is equal to the max overlap number
    #     print('delta', delta,'s[first_non_zero]', s_[first_non_zero], 'first_noon_zero', first_non_zero)
    s = np.where(s_ > 0, np.exp(s_ - delta), np.nan)

    return s


def glue_s(bins, counts, s_left, s_right):
    """
    gluing together 3 pieces of S:  s_left, counts=exp(s) and s_right.

    first we need to identify the glue points, after -- rescale the pieces to match the 'counts' and
    finally glue them together
    """

    left_index = np.argmax(bins == len(s_left) - 1)  # index in bins for the last_s [-1] number of overlaps

    rescaled_s_left = rescale_s(s_left, bins, counts, bin_index=left_index, s_index=len(s_left) - 1)

    first_non_zero = np.argmax(np.array(s_right) > 0)
    right_index = np.argmax(bins == first_non_zero)  # index in bins for the last_s [-1] number of overlaps

    rescaled_s_right = rescale_s(s_right, bins, counts, bin_index=right_index, s_index=first_non_zero)

    if (left_index > right_index):
        delta = left_index - right_index + 2
    else:
        delta = 0

    print('delta', delta)
    total_s = np.concatenate((rescaled_s_left[:len(rescaled_s_left) + 1 - delta],
                              counts[left_index + 1:right_index],
                              rescaled_s_right[first_non_zero:])
                             ).astype('float')
    total_s[total_s == 0.0] = np.nan
    # norming
    total_s = total_s / np.nansum(total_s)

    return total_s


def n_conf(N, dx, dy, dz):
    """
    calculates the number of conformations  of ideal grid polymer given
    number of bonds and displacements along the grid.
    """
    dx = abs(dx); dy = abs(dy); dz = abs(dz)


    if ((N - dx - dy + dz) % 2 != 0) | ((N - dx - dy - dz) % 2 != 0):
        return 0
    else:

        n_plus = int((N - dx - dy + dz) / 2)
        n_minus = int((N - dx - dy - dz) / 2)

        numerator = math.factorial(N)
        res = 0.0
        for x in range(n_minus + 1):
            for y in range(n_minus - x + 1):
                res += numerator / math.factorial(x) / math.factorial(x + dx) / math.factorial(y) / math.factorial(
                    y + dy) / \
                       math.factorial(n_plus - x - y) / math.factorial(n_minus - x - y)

        return res



def cache_n_conf(N_, dx, dy, dz):
    """
    caches the n_conf for each point on the grid given by dz, dy, dz
    """
    res = []
    for n in range(N_):
        for i in range(dx):
            for j in range(dy):
                for k in range(dz):
                    res.append(n_conf(n + 1, i, j, k))
                    # print(n+1, i, j, k, n_conf(n+1,i,j,k))
    return np.array(res).reshape(N_, dx, dy, dz)  # .astype(int)
    # return re



def calculate_saw_fraction(path):
    """
    produces fraction of SAW for a given ``path`` where results of a single experiment  are stored. e.g 'results/10/run_1'.
    The glued `s` can contain NANs, since `s_right` can be calculated not for  every overlap in order to speed up convergence.
    The right-hand side tail looks like  a straight line, so we fit the right tail with a line and filling the missed data
    (NANs) with the data after fitting a line.
    """

    def func(x, a, b):
        return a * x + b

    logging.info("extracting `n` from the path")
    n = int(path.split('/')[1])
    logging.info("n=%i" %n)
    logging.info('loading  saved results...')
    s_left, s_right, s_total, bins, counts, metrics = load_data(path)
    logging.info("done loading saved data")

    #     print('sum of s_total', np.nansum(s_total), path)
    first_index = np.nanargmax(np.isnan(s_total)) -1 # the  index before the  first  NAN in s_total
    if first_index > 0: # it there are NaNs!
        last_index = len(s_total)
        logging.info('first  index to fit %i, last is  %i'%(first_index, last_index))
        y_data = s_total[first_index:last_index]
        x_data = range(len(y_data))
        real_data_length = len(x_data)

        # filtering NaNs from y_data and keep sync with  x_data
        # it scipy can not fit curve with nans
        logging.info("removing NANs from data to fit...")
        data = [(p[0], p[1]) for p in zip(y_data, x_data) if not p[0] != p[0]]

        y_data, x_data = list(zip(*data))
        y_data = np.log(np.array(y_data))  # transform to ln
        x_data = np.array(x_data)

        logging.info("data is ready for fitting.")
        # fitting
        refilled_data = []
        try:
            popt, pcov = curve_fit(func, x_data, y_data)
            # adding points missed
            refilled_data = [func(x, *popt) for x in range(real_data_length)]
            refilled_data = np.exp(refilled_data)
        except:
            pass
        logging.info("joining back filled piece with added points")

        s_total_joined = np.concatenate((s_total[:first_index], refilled_data, s_total[last_index:]))

        # renorming
        s_total_joined_normed = s_total_joined / np.sum(s_total_joined)
        return round(1 / n, 4), np.log(s_total_joined_normed[0]) / n, s_total_joined_normed
    else:
        return round(1 / n, 4), np.log(s_total[0]) / n, s_total

