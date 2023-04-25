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




def regrow_saw(n, dx, dy, dz, res, w, alpha, k, prob, prob_native, limits ,use_limits):
    """
    recursive  regrow for SAW only.
    """

    if n == 1:
        #print(indexes[0][box], indexes[1][box], indexes[2][box])
        #contact_list = [[ dx - 1, dy, dz],
        #              [ dx + 1, dy, dz],
        #              [ dx, dy - 1, dz],
        #              [ dx, dy + 1, dz],
        #              [ dx, dy, dz - 1],
        #              [ dx, dy, dz + 1]
        #            ]
        #print(contact_list)
        #sys.exit()

        contact_list = [
                 [1 ,0, 0],\
                 [-1 ,0, 0],\
                 [0 ,1, 0],\
                 [0 ,-1, 0],\
                 [0 ,0, 1],\
                 [0 ,0, -1]\
               ]
        number_of_contacts = len([el for el in contact_list if el in res ])-1
        return res + [[0, 0, 0]], w, k+number_of_contacts, prob, prob_native

    else:
        neighbours = [[n - 1, dx - 1, dy, dz],
                      [n - 1, dx + 1, dy, dz],
                      [n - 1, dx, dy - 1, dz],
                      [n - 1, dx, dy + 1, dz],
                      [n - 1, dx, dy, dz - 1],
                      [n - 1, dx, dy, dz + 1]
                    ]
        if len(res)>2: 
           # removing the last in the result
           neighbours = [el for el in neighbours if el[1:] != res[-2]]
        counts = []
        counts_native = []
        tmp = 0
        number_non_overlap  = 0
        collect_closeness_to_axis=[]
        
        for neighbour in neighbours:

            # if there is already such a point, set n_coincide=0 to filter out that trial
            n_coincide = 1
            is_outside = 1
            number_of_contacts = 0

            if (neighbour[1:] in res) or (neighbour[1:] == [0, 0, 0]):
                n_coincide = 0
            else:
               contact_list = [
                 [neighbour[1]+1 ,neighbour[2], neighbour[3]],\
                 [neighbour[1]-1 ,neighbour[2], neighbour[3]],\
                 [neighbour[1] ,neighbour[2]+1, neighbour[3]],\
                 [neighbour[1] ,neighbour[2]-1, neighbour[3]],\
                 [neighbour[1] ,neighbour[2], neighbour[3]+1],\
                 [neighbour[1] ,neighbour[2], neighbour[3]-1]\
               ]
               contact_array = [el for el in contact_list if (el in res) and (el != [0,0,0]) ]
               number_of_contacts = len(contact_array)-1
               #if len(res) >2:
               #  is_minus2_in_contact = 1 if res[-3] in contact_array else 0
               #else:
               #  is_minus2_in_contact = 0
            if use_limits:
              if  (aux.is_in_box(np.array(res +  [neighbour[1:]]), limits=limits)):
                 is_inside = 1
              else:
                 is_inside = 0
            else:
              is_inside = 1

            collect_closeness_to_axis.append(number_of_contacts)

            #mask.append(n_coincide)

            # checking if outside the box
            # is_inside = 1  # aux.is_inside_box(*neighbour[1:])
            # sum of distances to OX axis

            #out_of_box = sum(np.where(t1 > t2 , t1-t2, 0))
            #out_of_box = (abs(neighbour[1]) + abs(neighbour[2]) + abs(neighbour[3]))/3



            #             print(coords, n_coincide)
            count = consts.caches[neighbour[0] - 1, abs(neighbour[1]), abs(neighbour[2]), abs(neighbour[3])]
            
            counts_native.append(n_coincide)
            
            if (count > 0) and (n_coincide == 1):
                number_non_overlap += 1
                #count =1
            
            tmp += count #* n_coincide  # accumulating the denominator

            counts.append(count * n_coincide * is_inside*  np.exp(alpha*number_of_contacts))
            #counts.append(count  *  is_inside* np.exp(alpha*number_of_contacts))

           # counts.append(count * n_coincide)
        #print(counts, counts_native)
        #sys.exit()
        
        
        if tmp == 0:
            # print('failed to grow... restarting')
            return res + [[0, 0, 0]], w, k, prob, prob_native
        # calculating W
        #w = w * sum(counts) / tmp
        #         print(w)
        #counts_native = np.multiply(mask, counts_native) # leaving just normilized entries for SAWs.
        #counts = np.multiply(mask, counts) # leaving just normilized entries for SAWs.
        #counts = np.multiply(counts, [np.exp(-el*alpha) for el in collect_closeness_to_axis]) # leaving just normilized entries for SAWs.
        
        if sum(counts) == 0:
            #print('failed to grow')
            #sys.exit()
            return res + [[0, 0, 0]], w, k, prob, prob_native
        
        # normalising p_i
        #counts_native = [c / sum(counts_native) for c in counts_native]
        counts_ = [c / sum(counts) for c in counts]


        # making cumulative sum
        counts__ = np.cumsum(counts_)

        # selecting one of neigbours
        selected = np.argmax(counts__ > np.random.random())#*sum(counts_))

        # if overlap selected -- return
        if counts_native[selected] == 0:
           return res + [[0, 0, 0]], w, k, prob, prob_native

        # if tmp ==0: print('all zeros, no saws', selected, res, neighbours[selected][1:])
        res.append(neighbours[selected][1:])

        # if the grown part escapes the box --> exit
        if (not (aux.is_in_box(np.array(res), limits=limits))) and  use_limits:
            return res + [[0, 0, 0]], w, k, prob, prob_native

        #print('counts native  %s' %counts_native)

        #print('selected ', selected)
        #if (number_non_overlap == 1) and (counts_[selected] != 1):
        #      logging.error("grow errori.Number non everlap is %i"%number_non_overlap)
        #      sys.exit()
        # configuration probability
        prob = prob  *  counts_[selected] #* number_non_overlap 
        #if (prob > 1.01) and (alpha == 0.0):
        #   print('prob: %f, counts   %s, %s; number non overlaps: %i, alpha: %f' %(prob, counts_,  counts, number_non_overlap, alpha))
        #   sys.exit()
        prob_native = prob_native + number_non_overlap #* counts_native[selected]
        if counts[selected] == 0:
            print(counts, selected)
            print('SELECTED NOT POSSIBLE')

        k += collect_closeness_to_axis[selected]
        return regrow_saw(*neighbours[selected], res, w, alpha, k, prob, prob_native, limits, use_limits)


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
