"""
Changing chain conformation
"""

import logging
import sys

import numpy as np

try:
    import consts, aux
except:
    from chromosome_segregation import consts, aux


def regrow(n, dx, dy, dz, res):
    """
    recursive uniform regrow.
    The generated coordinates are uniform random
    """

    if n == 1:
        return res + [[0, 0, 0]]
    else:
        neighbours = [[n - 1, dx - 1, dy, dz],
                      [n - 1, dx + 1, dy, dz],
                      [n - 1, dx, dy - 1, dz],
                      [n - 1, dx, dy + 1, dz],
                      [n - 1, dx, dy, dz - 1],
                      [n - 1, dx, dy, dz + 1]
                      ]
        counts = []
        for n_ in neighbours:

            # checking if outside the box

            is_inside = aux.is_inside_box(*n_[1:])

            if is_inside == 1:
                counts.append(consts.caches[n_[0] - 1, abs(n_[1]), abs(n_[2]), abs(n_[3])])
            else:
                counts.append(0)

        # normalising
        counts = [c / sum(counts) for c in counts]
        # uniform sampling
        #         nonzeros = [v for v in counts if v>0]
        #         counts = [1/len(nonzeros) if c !=0 else c for c in counts]
        # wrong sampling
        #         counts = [1/6. for i in range(6)]

        #         print(counts)

        # making cumulative sum
        counts_ = np.cumsum(counts)
        if sum(counts_) == 0:
            print('failed to grow')
            sys.exit()
        # selecting one of neigbours
        selected = np.argmax(counts_ > np.random.rand())
        #         print(neighbours[selected])
        res.append(neighbours[selected][1:])

        return regrow(*neighbours[selected], res)


def regrow_biased(n, dx, dy, dz, res, w, alpha, k):
    """
    recursive biased regrow.
    The generated coordinates are biased
    """

    if n == 1:
        return res + [[0, 0, 0]], w, k
    else:
        neighbours = [[n - 1, dx - 1, dy, dz],
                      [n - 1, dx + 1, dy, dz],
                      [n - 1, dx, dy - 1, dz],
                      [n - 1, dx, dy + 1, dz],
                      [n - 1, dx, dy, dz - 1],
                      [n - 1, dx, dy, dz + 1]
                      ]
        counts = []
        tmp = 0
        all_coincidence = []
        for neighbour in neighbours:

            # Lets do the following: a trial coordinate can bring 0 or 1 as an additional contribution
            # to the number of coordinate coinsidence.
            n_coincide = 0
            if (neighbour[1:] in res) or neighbour[1:] == [0, 0, 0]:
                n_coincide = 1
            all_coincidence.append(n_coincide)

            # checking if outside the box
            is_inside = aux.is_inside_box(*neighbour[1:])
            # if (is_inside ==0): print("OUTSiDE")


            #             print(coords, n_coincide)
            count = consts.caches[neighbour[0] - 1, abs(neighbour[1]), abs(neighbour[2]), abs(neighbour[3])]
            counts.append(np.exp(-alpha * n_coincide) * count * is_inside)
            tmp += count * is_inside  # accumulating the denominator

        # calculating W
        w = w * sum(counts) / tmp
        #         print(w)
        # normalising p_i
        counts = [c / sum(counts) for c in counts]

        # making cumulative sum
        counts_ = np.cumsum(counts)
        if sum(counts_) == 0:
            print('failed to grow')
            sys.exit()
        # selecting one of neigbours
        selected = np.argmax(counts_ > np.random.random())
        res.append(neighbours[selected][1:])
        if counts[selected] == 0:
            print('SELECTED NOT POSSIBLE')
        k += all_coincidence[selected]
        return regrow_biased(*neighbours[selected], res, w, alpha, k)


def regrow_saw(n, dx, dy, dz, res, w, alpha, k ,prob):
    """
    recursive  regrow for SAW only.
    """

    if n == 1:
        return res + [[0, 0, 0]], w, k, prob
    else:
        neighbours = [[n - 1, dx - 1, dy, dz],
                      [n - 1, dx + 1, dy, dz],
                      [n - 1, dx, dy - 1, dz],
                      [n - 1, dx, dy + 1, dz],
                      [n - 1, dx, dy, dz - 1],
                      [n - 1, dx, dy, dz + 1]
                      ]
        counts = []
        tmp = 0
        all_coincidence = []
        collect_closeness_to_axis=[]
        for neighbour in neighbours:

            # if there is already such a point, set n_coincide=0 to filter out that trial
            n_coincide = 1
            # if (neighbour[1:] in res) or neighbour[1:] == [0, 0, 0]:
            #     n_coincide = 0
            all_coincidence.append(n_coincide)

            # checking if outside the box
            # is_inside = 1  # aux.is_inside_box(*neighbour[1:])
            # sum of distances to OX axis

            closeness_to_axis = abs(neighbour[2]) + abs(neighbour[3])


            collect_closeness_to_axis.append(closeness_to_axis)

            #             print(coords, n_coincide)
            count = consts.caches[neighbour[0] - 1, abs(neighbour[1]), abs(neighbour[2]), abs(neighbour[3])]
            # if count >1: count=1
            counts.append(count * n_coincide * np.exp(-alpha*closeness_to_axis))

            tmp += count #* n_coincide  # accumulating the denominator
        if tmp == 0:
            # print('failed to grow... restarting')
            return res + [[0, 0, 0]], w, k, prob
        # calculating W
        w = w * sum(counts) / tmp
        # w = w + np.log(sum(counts) / tmp)
        #         print(w)

        if sum(counts) == 0:
            return res + [[0, 0, 0]], w, k, prob
            # print('failed to grow')
            # sys.exit()

        # normalising p_i
        counts = [c / sum(counts) for c in counts]

        # making cumulative sum
        counts_ = np.cumsum(counts)


        # selecting one of neigbours
        selected = np.argmax(counts_ > np.random.random())
        # if tmp ==0: print('all zeros, no saws', selected, res, neighbours[selected][1:])
        res.append(neighbours[selected][1:])

        prob = prob*counts[selected]
        # prob  = prob + np.log(counts[selected])

        if counts[selected] == 0:
            print('SELECTED NOT POSSIBLE')

        k += collect_closeness_to_axis[selected]
        return regrow_saw(*neighbours[selected], res, w, alpha, k, prob)


def regrow_saw_segregation_prove(n, dx, dy, dz, res, w, alpha, k, coords):
    """
    recursive  regrow for SAW only.
    """

    if n == 1:
        return res + [[0, 0, 0]], w, k
    else:
        neighbours = [[n - 1, dx - 1, dy, dz],
                      [n - 1, dx + 1, dy, dz],
                      [n - 1, dx, dy - 1, dz],
                      [n - 1, dx, dy + 1, dz],
                      [n - 1, dx, dy, dz - 1],
                      [n - 1, dx, dy, dz + 1]
                      ]
        counts = []
        tmp = 0
        all_coincidence = []
        for neighbour in neighbours:

            # if there is already such a point, set n_coincide=0 to filter out that trial
            n_coincide = 1
            if (neighbour[1:] in res) or (neighbour[1:] == [0, 0, 0]) or (neighbour[1:] in coords):
                n_coincide = 0
            all_coincidence.append(n_coincide)
            # checking if outside the box
            is_inside = aux.is_inside_box(*neighbour[1:])

            count = consts.caches[neighbour[0] - 1, abs(neighbour[1]), abs(neighbour[2]), abs(neighbour[3])]

            counts.append(count * n_coincide * is_inside)
            tmp += count * n_coincide * is_inside  # accumulating the denominator
        if sum(counts) == 0:
            logging.info("failed to grow... restarting")
            # print('failed to grow... restarting')
            return res + [[0, 0, 0]], w, k
        # calculating W
        # w = w * sum(counts) / tmp
        #         print(w)

        # normalising p_i
        # counts = [c / sum(counts) for c in counts]

        # making cumulative sum
        counts_ = np.cumsum([c / sum(counts) for c in counts])

        if sum(counts_) == 0:
            print('failed to grow')
            sys.exit()
        # selecting one of neighbours
        selected = np.argmax(counts_ > np.random.random())
        # if tmp ==0: print('all zeros, no saws', selected, res, neighbours[selected][1:])
        res.append(neighbours[selected][1:])

        # calculating the probability of the configuration
        w = w * counts[selected] / tmp

        if counts[selected] == 0:
            print('SELECTED NOT POSSIBLE')

        k += all_coincidence[selected]
        return regrow_saw_segregation_prove(*neighbours[selected], res, w, alpha, k, coords)



def build_one_monomer(neighbours, res1, res2):

    counts = []
    tmp = 0
    all_coincidence = []
    for neighbour in neighbours:

        # if there is already such a point, set n_coincide=0 to filter out that trial
        n_coincide = 1
        if (neighbour[1:] in res1) or (neighbour[1:] == [0, 0, 0]) or (neighbour[1:] in res2):
            n_coincide = 0
        all_coincidence.append(n_coincide)
        # checking if outside the box
        is_inside = aux.is_inside_box(*neighbour[1:])

        count = consts.caches[neighbour[0] - 1, abs(neighbour[1]), abs(neighbour[2]), abs(neighbour[3])]

        counts.append(count * n_coincide * is_inside)
        tmp += count * n_coincide * is_inside  # accumulating the denominator

    if sum(counts) == 0:
        logging.debug("failed to grow... restarting")
        return -1 #res1 + [[0, 0, 0]]

    # making cumulative sum
    counts_ = np.cumsum([c / sum(counts) for c in counts])

    # selecting one of neighbours
    selected = np.argmax(counts_ > np.random.random())

    return selected




def regrow_saw_segregation_prove_two_chains(n1, dx1, dy1, dz1, n2, dx2, dy2, dz2, res1, res2):
    """
    recursive SAW regrow  for two chains simultaneously
    """

    if n1 == 1:
        return res1 + [[0, 0, 0]], res2 + [[0, 0, 0]] #w, k
    else:
        neighbours1 = [[n1 - 1, dx1 - 1, dy1, dz1],
                      [n1 - 1, dx1 + 1, dy1, dz1],
                      [n1 - 1, dx1, dy1 - 1, dz1],
                      [n1 - 1, dx1, dy1 + 1, dz1],
                      [n1 - 1, dx1, dy1, dz1 - 1],
                      [n1 - 1, dx1, dy1, dz1 + 1]
                      ]
        neighbours2 = [[n2 - 1, dx2 - 1, dy2, dz2],
                       [n2 - 1, dx2 + 1, dy2, dz2],
                       [n2 - 1, dx2, dy2 - 1, dz2],
                       [n2 - 1, dx2, dy2 + 1, dz2],
                       [n2 - 1, dx2, dy2, dz2 - 1],
                       [n2 - 1, dx2, dy2, dz2 + 1]
                       ]
        # growing a monomer for  one  chain
        selected1 = build_one_monomer(neighbours1, res1, res2)
        if  selected1 >-1:
            res1.append(neighbours1[selected1][1:])
        else:
            logging.debug("selected1 is: %i" %selected1)
            return res1 + [[0, 0, 0]], res2 + [[0, 0, 0]]

        # growing a monomer for  the another  chain
        selected2 = build_one_monomer(neighbours2, res1, res2)
        if selected2 >-1:
            res2.append(neighbours2[selected2][1:])
        else:
            return res1 + [[0, 0, 0]], res2 + [[0, 0, 0]]


        return regrow_saw_segregation_prove_two_chains(*neighbours1[selected1], *neighbours2[selected2],  res1, res2)

# """
# Changing chain conformation
# """
#
# import math
# import  numpy as np
# import  sys
# try:
#     import consts, aux
# except:
#     from chromosome_segregation import consts, aux

# def regrow(n, dx, dy, dz, res):
#     """
#     recursive uniform regrow.
#     The generated coordinates are uniform random
#     """
#
#     if n == 1:
#         return res + [[0, 0, 0]]
#     else:
#         neighbours = [[n - 1, dx - 1, dy, dz],
#                       [n - 1, dx + 1, dy, dz],
#                       [n - 1, dx, dy - 1, dz],
#                       [n - 1, dx, dy + 1, dz],
#                       [n - 1, dx, dy, dz - 1],
#                       [n - 1, dx, dy, dz + 1]
#                       ]
#         counts = []
#         for n_ in neighbours:
#
#             # checking if outside the box
#             is_inside = aux.is_inside_box(*n_[1:])
#
#             if is_inside ==1:
#                 counts.append(consts.caches[n_[0] - 1, abs(n_[1]), abs(n_[2]), abs(n_[3])])
#             else:
#                 counts.append(0)
#
#         # normalising
#         counts = [c / sum(counts) for c in counts]
#         # uniform sampling
#         #         nonzeros = [v for v in counts if v>0]
#         #         counts = [1/len(nonzeros) if c !=0 else c for c in counts]
#         # wrong sampling
#         #         counts = [1/6. for i in range(6)]
#
#         #         print(counts)
#
#         # making cumulative sum
#         counts_ = np.cumsum(counts)
#         if sum(counts_) == 0:
#             print('failed to grow')
#             sys.exit()
#         # selecting one of neigbours
#         selected = np.argmax(counts_ > np.random.rand())
#         #         print(neighbours[selected])
#         res.append(neighbours[selected][1:])
#
#
#         return regrow(*neighbours[selected], res)



# def regrow_biased(n, dx, dy, dz, res, w, alpha, k):
#     """
#     recursive biased regrow.
#     The generated coordinates are biased
#     """
#
#     if n == 1:
#         return res + [[0, 0, 0]], w, k
#     else:
#         neighbours = [[n - 1, dx - 1, dy, dz],
#                       [n - 1, dx + 1, dy, dz],
#                       [n - 1, dx, dy - 1, dz],
#                       [n - 1, dx, dy + 1, dz],
#                       [n - 1, dx, dy, dz - 1],
#                       [n - 1, dx, dy, dz + 1]
#                       ]
#         counts = []
#         tmp = 0
#         all_coincidence = []
#         for neighbour in neighbours:
#
#             # Lets do the following: a trial coordinate can bring 0 or 1 as an additional contribution
#             # to the number of coordinate coinsidence.
#             n_coincide = 0
#             if neighbour[1:] in res:
#                 n_coincide = 1
#             all_coincidence.append(n_coincide)
#
#             # checking if outside the box
#             is_inside = aux.is_inside_box(*neighbour[1:])
#             # if (is_inside ==0): print("OUTSiDE")
#
#
#             #             print(coords, n_coincide)
#             count = consts.caches[neighbour[0] - 1, abs(neighbour[1]), abs(neighbour[2]), abs(neighbour[3])]
#             counts.append(np.exp(-alpha * n_coincide) * count * is_inside)
#             tmp += count* is_inside  # accumulating the denominator
#
#         # calculating W
#         w = w * sum(counts) / tmp
#         #         print(w)
#         # normalising p_i
#         counts = [c / sum(counts) for c in counts]
#
#         # making cumulative sum
#         counts_ = np.cumsum(counts)
#         if sum(counts_) == 0:
#             print('failed to grow')
#             sys.exit()
#         # selecting one of neigbours
#         selected = np.argmax(counts_ > np.random.random())
#         res.append(neighbours[selected][1:])
#         if counts[selected] == 0:
#             print('SELECTED NOT POSSIBLE')
#         k += all_coincidence[selected]
#         return regrow_biased(*neighbours[selected], res, w, alpha, k)
