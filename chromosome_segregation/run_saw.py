import  numpy as np

from  chromosome_segregation.simulation import URW_saw
from chromosome_segregation.aux import cache_n_conf

from chromosome_segregation  import overlaps

from  chromosome_segregation import  consts

import  matplotlib.pyplot as plt

from statistics import Counter


def get_n(box, density):
    t = density*box[0]*box[1]*box[2]*8
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


def plot():
    plt.figure(figsize=(12, 8))
    # thickness=3

    for result, thickness in zip(total_results, (1, 2, 3)):
        x = []
        y = []
        boxes = [[i, thickness, thickness] for i in range(n_boxes, 0, -1)]
        for i in range(len(boxes)):
            n = get_n(boxes[i])
            x.append(1 - n / get_n(boxes[0]))
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

        max_n = get_n([n_boxes, thicknesses[-1], thicknesses[-1]], density)
        print('max_n=%i' % (max_n))

        consts.caches = cache_n_conf(N_=max_n + 1, dx=20, dy=20, dz=20)

        total_results = []
        for thickness in thicknesses:
            boxes = [[i, thickness, thickness] for i in range(n_boxes, 0, -1)]
            print('boxes: %s' % boxes)
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

        plot()


if __name__ == "__main__":

    density = 0.3
    n_boxes = 3
    thicknesses = list(range(1, 4))

    run(density, n_boxes, thicknesses)

