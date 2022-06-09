"""
Helpers
"""
import  os
import  json
import  numpy as np



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


