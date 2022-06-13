"""
Changing chain conformation
"""

import math
import  numpy as np
import  sys



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
            counts.append(caches[n_[0] - 1, abs(n_[1]), abs(n_[2]), abs(n_[3])])

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
        #         print(c)
        #         print('res', res)
        #         sel.append(selected)

        #         nz.append(counts[selected])


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
            if neighbour[1:] in res:
                n_coincide = 1

                # calculating number of coordinate coincidence  for the  grown part of the coords
            #             coords = res+[tuple(neighbour[1:]) , (0,0,0)] # adding a trial position
            #             c = Counter(coords).values()

            #             coords = np.array(res+[neighbour[1:] , [0,0,0]]) # adding a trial position
            #             u, c = np.unique(coords, axis=0, return_counts=True)
            # number of coincidential coords (NOT the number of pair-wise intersections)
            #             n_coincide = sum([el for el in c if  el>1])/len(coords)
            #             n_coincide = sum([comb(v,2) for v in c])/len(coords)
            all_coincidence.append(n_coincide)
            #             print(coords, n_coincide)
            count = caches[neighbour[0] - 1, abs(neighbour[1]), abs(neighbour[2]), abs(neighbour[3])]
            counts.append(np.exp(-alpha * n_coincide) * count)
            tmp += count  # accumulating the denominator
        # print(count, np.exp(-alpha*n_coincide))

        #         print(10*'dd')
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
        selected = np.argmax(counts_ > np.random.rand())
        #         print(neighbours[selected])
        res.append(neighbours[selected][1:])
        #         print(c)
        #         print('res', res)
        #         sel.append(selected)

        #         nz.append(counts[selected])
        #         print(all_coincidence[selected])
        k += all_coincidence[selected]
        #         print(len(res), all_coincidence[selected])
        return regrow_biased(*neighbours[selected], res, w, alpha, k)

# sys.exit()
