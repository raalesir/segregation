"""
Helpers
"""
import  os
import  json
import  numpy as np
import gzip

import  math
import  logging
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import  statistics
from collections import  Counter

try:
    import consts
except:
    from chromosome_segregation import consts

RESULTS_FOLDER = consts.RESULTS_FOLDER

ln_fact_cache = []

def list_to_arr(boxes):
    items = Counter(boxes).items()
#     keys = Counter(boxes).keys()
    arr = np.zeros(30)
    for item in Counter(boxes).items():
#     print(item)
        arr[item[0]] = item[1]
    return arr/np.sum(arr)


def is_inside_box(x,y,z):
    """
    checks if a  point inside  the box
    :param x:
    :type x:
    :param y:
    :type y:
    :param z:
    :type z:
    :return:
    :rtype:
    """

    if consts.RUN_FREE:
        return  1
    elif consts.RUN_HALFSPACE:
        if (abs(x)<=consts.size_x/2) & \
                (abs(y)  <= consts.size_y/2) & \
                (z <=consts.size_z)  &\
                (z>=0):
            return 1
        else:
            return 0
    else:
        if (abs(x) <= consts.size_x / 2) & \
                (abs(y) <= consts.size_y / 2) & \
                (abs(z) <= consts.size_z / 2):
            return 1
        else:
            return 0


def plot_specific_entropy_comparison(x, y, errs, x_half, y_half, errs_half, save_to='.', fname='tmp.png'):
    def func1(x, a, b, c):
        #     return a * np.exp(-b * x)
        #     return a*x +b
        return a * x * np.log(x) + b * x + c

    popt, pcov = curve_fit(func1, x, y)
    print(popt)
    xx = np.linspace(0.0001, 0.15, 100)
    plt.rc('font', size=12)

    plt.figure(figsize=(16, 10))

    plt.plot(xx, func1(xx, *popt),
             label=' fit to data: $f(x) = %1.4f [x\ln(x)] +%1.4f x %1.4f$' % (popt[0], popt[1], popt[2]))
    plt.scatter(x, y, facecolor='red', edgecolor='black', s=100, label='simulation results')
    plt.errorbar(x, y, yerr=errs, fmt=".", color='black', capsize=5)

    popt, pcov = curve_fit(func1, x_half, y_half)
    plt.plot(xx, func1(xx, *popt), linestyle='-.', label=' fit for half-space')
    plt.scatter(x_half, y_half, facecolor='blue', alpha=0.5, edgecolor='black', s=100,
                label='simulation results for half-space')
    plt.errorbar(x_half, y_half, yerr=errs_half, fmt=".", color='black', capsize=5)

    plt.scatter(0, -0.2476, s=200, label='limit for $n\\to\infty$', marker='X', color='black')
    plt.scatter(1 / 12., -0.31273, s=250, label='exact for n=12', marker='8', facecolor='none', edgecolor='black')
    plt.scatter(1 / 10., -0.31907, s=250, label='exact for  n=10', marker='D', facecolor='none', edgecolor='black')
    plt.scatter(1 / 8., -0.32539, s=250, label='exact for  n=8', marker='s', facecolor='none', edgecolor='black')

    plt.grid()
    plt.xlabel('1/N')
    plt.ylabel('$\Delta S/N$')
    plt.title("""Specific excess entropy $\Delta S/N$ for free polymer and a polymer confined in half-space.
    The limit value for scaling is -0.2476""")
    plt.legend()
    # plt.xticks(list(plt.xticks()[0])+ [.01])
    plt.xticks(np.arange(0, 0.15, .01))
    plt.xlim(-0.003, 0.15)

    plt.savefig(os.path.join(save_to, fname), )



def plot_specific_entropy(x,y,errs, save_to=RESULTS_FOLDER):
    def func1(x, a, b, c):
        #     return a * np.exp(-b * x)
        #     return a*x +b
        return a * x * np.log(x) + b * x + c

    popt, pcov = curve_fit(func1, x, y)
    print(popt)
    xx = np.linspace(0.0001, 0.15, 100)
    plt.rc('font', size=14)

    plt.figure(figsize=(12, 8))
    plt.plot(xx, func1(xx, *popt),
             label=' fit to data: $f(x) = %1.4f [x\ln(x)] +%1.4f x %1.4f$' % (popt[0], popt[1], popt[2]))
    plt.scatter(x, y, facecolor='red', edgecolor='black', s=100, label='simulation results')

    plt.errorbar(x, y, yerr=errs, fmt=".", color='black', capsize=5)

    plt.scatter(0, -0.2476, s=200, label='limit for $n\\to\infty$', marker='x', color='black')
    plt.scatter(1 / 12., -0.31273, s=250, label='exact for n=12', marker='8', facecolor='none', edgecolor='black')
    plt.scatter(1 / 10., -0.31907, s=250, label='exact for  n=10', marker='D', facecolor='none', edgecolor='black')
    plt.scatter(1 / 8., -0.32539, s=250, label='exact for  n=8', marker='s', facecolor='none', edgecolor='black')

    plt.grid()
    plt.xlabel('1/N')
    plt.ylabel('$\Delta S/N$')
    plt.title("""Specific excess entropy $\Delta S/N$.  The limit value for scaling is -0.2476""")
    plt.legend()
    # plt.xticks(list(plt.xticks()[0])+ [.01])
    plt.xticks(np.arange(0, 0.15, .01))
    plt.xlim(-0.003, 0.15)

    plt.savefig(os.path.join(save_to, 'specific_entropy.png'), )





def get_result_subfolders(path=RESULTS_FOLDER):
    lst = []
    for root, subdir, files in os.walk(path):
        if (len(subdir) >0):
            for sd  in subdir:
                if (sd.startswith('run')):
                    lst.append(os.path.join(root, sd))
    return lst



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



def save_results(n, s_left, s_right, s_total, bins, counts, metrics, saw_fraction, fitted_s):
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

    names = ['s_left.txt', 's_right.txt', 's_total.txt', 'bins.txt', 'counts.txt','fitted_s.txt', 'saw.txt']
    data = [s_left, s_right, s_total, bins, counts, fitted_s, np.array((round(1/n, 4), saw_fraction))]
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


def get_grow_caches(fname, params):
    if os.path.isfile(fname):
        with gzip.open(fname) as  f:
            header = f.readline().decode("utf-8")[2:-1]
#             print('raw header', header)
            header = tuple([int(el) for el in header.split(',')])
#             print('header', header)
        d = np.loadtxt(fname)
#         d = np.load(fname)
        try:
            d = d.reshape(header)
#             print(d.shape, 'fff')

        except: pass
        if any(el[0]<el[1] for el in zip(d.shape,params)):
            logging.info("need to recalculate")
            logging.info("shape of existing cache is %s, shape of query is: %s" % (str(d.shape), str(params)))

            d = cache_n_conf(*params)
            np.savetxt(fname, d.ravel(), header=','.join([str(el) for el in d.shape]))
        logging.info("shape of existing cache is %s, shape of query is: %s"%(str(d.shape), str(params)))
    else:
        d = cache_n_conf(*params)
#         print(d.shape)
        np.savetxt(fname, d.ravel(), header=','.join([str(el) for el in d.shape]))
#         np.save(fname, d)


    return d

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


def ln_fact(x):

        if x <= 1 :
            return  0
        else:
        # return sum([math.log(el) for el in range(1, x)])
            return .5*np.log(2*np.pi*x)+ x*(np.log(x)-1)

def cache_ln_factorial(n):
    """
    caching factorials
    :param n:
    :type n:
    :return:
    :rtype:
    """
    return  [ln_fact(el) for el in range(0,n+1)]

def n_conf_large(N, dx, dy, dz):
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

        numerator = ln_fact_cache[N]
        res = 0.0
        for x in range(n_minus + 1):
            for y in range(n_minus - x + 1):
                res += np.exp(
                    numerator -  ln_fact_cache[x]- ln_fact_cache[x + dx] - ln_fact_cache[y] - ln_fact_cache[y + dy] - \
                    ln_fact_cache[n_plus - x - y] - ln_fact_cache[n_minus - x - y]
                )

        return res


def cache_n_conf(N_, dx, dy, dz):
    """
    caches the n_conf for each point on the grid given by dz, dy, dz
    """
    res = []
    global ln_fact_cache

    ln_fact_cache = cache_ln_factorial(N_)

    for n in range(N_):
        for i in range(dx):
            for j in range(dy):
                for k in range(dz):
                    # res.append(n_conf(n + 1, i, j, k))
                    res.append(n_conf_large(n + 1, i, j, k))
        print('current %i out of %i'%(n, N_))
                    # print(n+1, i, j, k, n_conf(n+1,i,j,k))
    return np.array(res).reshape(N_, dx, dy, dz)  # .astype(int)
    # return re



def calculate_saw_fraction(n, s_total):
    """

    :param n:
    :type n:
    :param s_total:
    :type s_total:
    :return:
    :rtype:
    """

    def func(x, a, b):
        return a * x + b

    first_index = 0
    if n > 20:
        first_index = np.nanargmax(np.isnan(s_total)) - 1  # the  index before the  first  NAN in s_total
    if first_index > 0:  # it there are NaNs!
        last_index = len(s_total)
        logging.info('first  index to fit %i, last is  %i' % (first_index, last_index))
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



def calculate_saw_fraction_all(path):
    """
    produces fraction of SAW for a given ``path`` where results of a single experiment  are stored. e.g 'results/10/run_1'.
    The glued `s` can contain NANs, since `s_right` can be calculated not for  every overlap in order to speed up convergence.
    The right-hand side tail looks like  a straight line, so we fit the right tail with a line and filling the missed data
    (NANs) with the data after fitting a line.
    """


    logging.info("extracting `n` from the path")
    n = int(path.split('/')[1])
    logging.info("n=%i" %n)
    logging.info('loading  saved results...')
    s_left, s_right, s_total, bins, counts, metrics = load_data(path)
    logging.info("done loading saved data")

    #     print('sum of s_total', np.nansum(s_total), path)
    reverse_n, saw_fraction, s = calculate_saw_fraction(n, s_total=s_total)

    return  reverse_n, saw_fraction, s


def prepare_entropy_plot(paths):
    """
    reads results directory and for each N calculates statistics for a plot
    """
    res = []
    for folder in paths:
        n = int(folder.split('/')[1])
        if os.path.isfile(os.path.join(folder, 'saw.txt')):
            logging.info("open 'saw.txt' for reading reverse n and saw fraction")
            reverse_n, specific_entropy = np.loadtxt(os.path.join(folder, 'saw.txt'))
        else:
            logging.info("saw.txt is missing. recalculating")
            reverse_n, specific_entropy, s = calculate_saw_fraction_all(folder)
            logging.info("saving to saw.txt")
            np.savetxt(os.path.join(folder, 'saw.txt'), np.array((reverse_n, specific_entropy)))
            logging.info("saved...")

        logging.info("%s, %i, %4.3f, %f"%(folder, n, reverse_n, specific_entropy))
        res.append((reverse_n, specific_entropy))

    # for each N get list of values
    dict_ = {}
    for k, v in sorted(res):
        if k not in dict_:
            dict_[k] = []
        dict_[k].append(v)

    # get statistics on the values
    for k in dict_:
        if len(dict_[k]) > 1:
            dict_[k] = sum(dict_[k]) / len(dict_[k]), statistics.stdev(dict_[k])
        else:
            dict_[k] = sum(dict_[k]) / len(dict_[k]), 0

    x = np.array(list(dict_.keys()))
    y, errs = zip(*list(dict_.values()))

    y = np.array(y)
    errs = np.array(errs)

    return x, y, errs


# def get_box(sx, sy, sz, l):
#     variants = [(l[0], l[1], l[2]),
#                 (l[0], l[2], l[1]),
#                 (l[1], l[0], l[2]),
#                 (l[1], l[2], l[0]),
#                 (l[2], l[0], l[1]),
#                 (l[2], l[1], l[0])
#                 ]
#     min_box = 1000
#     for v in variants:
#         box_x = sx.index(v[0])
#         box_y = sy.index(v[1])
#         box_z = sz.index(v[2])
#         box = max(box_x, box_y, box_z)
#         if box < min_box:
#             min_box = box
#
#     return min_box


def get_box(sx, sy, sz, l):
    variants = [(l[0], l[1], l[2]),
                (l[0], l[2], l[1]),
                (l[1], l[0], l[2]),
                (l[1], l[2], l[0]),
                (l[2], l[0], l[1]),
                (l[2], l[1], l[0])
                ]
    min_box = 1000
    for v in variants:
        try:

            if v[0] < sx[0]:
                box_x = 0
            else:
                box_x = sx.index(v[0])

            if v[1] < sy[0]:
                box_y = 0
            else:
                box_y = sy.index(v[1])

            if v[2] < sz[0]:
                box_z = 0
            else:
                box_z = sz.index(v[2])

        # try:
        #     box_x = sx.index(v[0])
        #     box_y = sy.index(v[1])
        #     box_z = sz.index(v[2])
            box = max(box_x, box_y, box_z)

            if box < min_box:
                min_box = box
        except ValueError as e:
            # print(v)
            pass

    if min_box == 1000:
        return None
    else:
        return min_box



def get_indexes(box, extend_to_left=1, extend_to_right=5, length=30):
    min_ = min(box)
    xs = range(box[0] - min_, length + box[0] - min_)
    ys = range(box[1] - min_, length + box[1] - min_)
    zs = range(box[2] - min_, length + box[2] - min_)

    index_in_focus = xs.index(box[0])

    a, b = max(0, index_in_focus - extend_to_left), index_in_focus + extend_to_right
    return list(xs)[a:b + 1], list(ys)[a:b + 1], list(zs)[a:b + 1]



def get_en_n_to_infty(data, filters={'min_n': 20, 'area': 25.0}):
    """
    calculating specific free energy for  the limit for n to infty. Using scipy.curve_fit
    """
    min_energy = {}
    print(filters)
    area_filter = filters.pop('area', None)
    #     print(type(area_filter))

    for k, group in aggregated.groupby('area'):
        if area_filter:
            group = group[(group['n'] > filters['min_n']) & (group['area'] == area_filter)][['n', 'en_mean']].dropna()
        else:
            group = group[group['n'] > filters['min_n']][['n', 'en_mean']].dropna()

        # get fit parameters from scipy curve fit
        if len(group['en_mean']) > 0:
            try:
                parameters, covariance = curve_fit(fun, 1 / group['n'], group['en_mean'])
                min_energy[k] = parameters[1]
            except:
                pass

    return min_energy
