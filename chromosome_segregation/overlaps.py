"""
    calculates number of configurations for a ring grid phantom polymer and
    overlap distribution for the chain.
"""

import math
#from math import comb
from scipy.special import comb
from collections import Counter
import json
import os, sys
import matplotlib.pyplot as plt


class Overlap:
    """
    calculating number of overlaps  for a ring grid polymer of given N.
    N should be  even to make a ring polymer.
    """

    def __init__(self, N):
        """
        N -- number of monomers
        """
        self.n = N
        self.n_conformations = self.n_conform()
        self.overlaps_hist = None
        self.indexes = self.calculate_steps()
        # print('self.indexes created')
        self.dict = self.make_steps()
        print('self.dict is calculated. the length is %i'%len(self.indexes))
        self.encoded_conformations = None

        self.numbers_of_contact = None
        self.trajectories  = None

    #         self.keep_result = []




    def __str__(self):
        return "overlaps for %i beads has %i conformations" % (self.n, self.n_conformations)

    def n_conform(self):
        """
        """
        r = 0
        for i in range(self.n // 2 + 1):
            for j in range(self.n // 2 - i + 1):
                r = r + math.factorial(self.n) / math.factorial(i)**2 / math.factorial(j)**2 / math.factorial(
                    self.n // 2 - i - j)**2
        return r


    def fun(self, d, res):
        """
        returns a letter representation for a number combination.
        For example for  n=4 it  is (2,0,0) -->> ['i+i+i-i-', 'i+i-i+i-', 'i+i-i-i+', '-+-+', '-++-'  ,'--++']

        :param d: numeric representation of the closed trajectory e.g. (2,0,1)
        :type d: dict
        :param res: string representation for all possible configurations  for the given numeric representation, see above in the function doc
        :type res: str
        :return: string representation for all possible configurations
        :rtype:
        """
        if sum(d.values()) == 0: # if (0,0,0) that is we are at the initial point
            # Overlap.keep_result(res)
            yield res
        else:
            for k in [item for item in d.keys() if d[item] > 0]:
                r = res
                r += k

                tmp = d.copy()
                tmp[k] -= 1

                # https: // stackoverflow.com / questions / 38254304 / can - generators - be - recursive
                yield  from  self.fun(tmp, r)


    def calculate_steps(self):
        """
        given number of monomers, n, produce the indexes (i,j,k)
        as the number of steps to make in positive and negative direction
        """
        res = []
        for i in range(self.n // 2 + 1):
            for j in range(self.n // 2 - i + 1):
                res.append((i, j, self.n // 2 - i - j))
        return res


    def make_step(self, tup):
        """
        encodes  single index
        :return:
        :rtype:
        """
        res = []
        d = {}
        d['i+'] = tup[0]
        d['i-'] = tup[0]
        d['j+'] = tup[1]
        d['j-'] = tup[1]
        d['k+'] = tup[2]
        d['k-'] = tup[2]

        res.append(d)

        return res


    def make_steps(self):
        """
        encodes indexes in dict
        """
        res = []
        for tup in self.indexes:
            d = {}
            d['i+'] = tup[0]
            d['i-'] = tup[0]
            d['j+'] = tup[1]
            d['j-'] = tup[1]
            d['k+'] = tup[2]
            d['k-'] = tup[2]

            res.append(d)
        return res

    #     @static
    def keep_result(data):
        """
        """
        Overlap.keep_result.all.append(data)


    def calculate_all_conformations(self):
        Overlap.keep_result.all = []

        # tmp = []
        #
        # i = 0
        # for entry in self.dict:
        #     # generating trajectory representation  for each combination
        #     print('working with %i entry %s out of %i' %(i, entry, len(self.dict)))
        #     res = self.fun(entry, '')
        #     tmp.append(res)
        #     # self.fun(entry, '')
        #     print(next(res))
        #     if not res:
        #         print('weird conf')
        #     # Overlap.keep_result(res)
        #     print('number of trajectories generated: %i, %i ' %(len(Overlap.keep_result.all), len(tmp)))
        #     i +=1
        self.trajectories = [self.fun(entry, '') for entry in self.dict]
        print(self.trajectories)
        # return self.trajectories

        # yield [self.fun(entry, '') for entry in self.dict]


    def encode_single_conformation(self, conformation):
        """

        :param conformation:
        :type conformation:
        :return:
        :rtype:
        """

        conf_encoded = []
        start = [0, 0, 0]
        for symbol in [conformation[i:i + 2] for i in range(0, len(conformation), 2)]:
            if symbol == 'k+':
                start[2] += 1
            elif symbol == 'k-':
                start[2] -= 1
            elif symbol == 'i+':
                start[0] += 1
            elif symbol == 'i-':
                start[0] -= 1
            elif symbol == 'j+':
                start[1] += 1
            elif symbol == 'j-':
                start[1] -= 1
            conf_encoded.append(tuple(start.copy()))

        return conf_encoded



    def tmpfun(self, conformation):
        """

        :return:
        :rtype:
        """

        conf_encoded = []
        start = [0, 0, 0]
        for symbol in [conformation[i:i + 2] for i in range(0, len(conformation), 2)]:
            if symbol == 'k+':
                start[2] += 1
            elif symbol == 'k-':
                start[2] -= 1
            elif symbol == 'i+':
                start[0] += 1
            elif symbol == 'i-':
                start[0] -= 1
            elif symbol == 'j+':
                start[1] += 1
            elif symbol == 'j-':
                start[1] -= 1
            conf_encoded.append(tuple(start.copy()))

        return conf_encoded


    def encode_to_coords(self):
        """
        encodes string conformation representation to the numeric representation
        e.g. 'i+i+i-i-' --> [ (1,0,0), (2,0,0), (1,0,0), (0,0,0) ]
        """

        res = (self.tmpfun(conformation) for trajectory in self.trajectories  for conformation in trajectory )# Overlap.keep_result.all)
        # res = [self.tmpfun(conformation) for conformation in self.trajectories]
        # i = 0
        # res = []
        # print('Overlap.keep_result.all', len(Overlap.keep_result.all))

        # for conformation in Overlap.keep_result.all:
        # for trajectory in self.trajectories:
        #     for conformation in trajectory:
        #     # print('conformation', conformation)
        #         if i%10000 == 0 :
        #             print('%5.2f' %(i/self.n_conformations * 100), end='\r')
        #         i+=1
        #         res.append(self.tmpfun(conformation) )

        return res




    def get_overlaps(self):
        """
        """
       # overlaps = []
       # ncs = []

        self.overlaps_hist = {}
        self.numbers_of_contact = {}

        length = self.n_conformations # len(Overlap.keep_result.all)
        i=0
        for conf in self.encoded_conformations:
             # conf is like [ (1,0,0), (2,0,0), (1,0,0), (0,0,0) ]
            if i%10000 == 0:
                print('passed', '%6.4f' %(i/length * 100) ,end='\r')
            i +=1
            number_overlaps = sum([comb(lst, 2) for lst in Counter(conf).values()])

            if number_overlaps in self.overlaps_hist:
                self.overlaps_hist[number_overlaps] +=1
            else:
                self.overlaps_hist[number_overlaps] = 1

            #overlaps.append(
            #    number_overlaps
            #)
            if number_overlaps == 0:
              nc = self.get_number_of_contacts(conf)
              if nc  in self.numbers_of_contact:
                    self.numbers_of_contact[nc] +=1
              else:
                    self.numbers_of_contact[nc] = 1
         #   ncs.append(self.get_number_of_contacts(conf))


        #counts = Counter(overlaps)
        #self.overlaps_hist = dict(counts)

        #self.numbers_of_contact = dict(Counter(ncs))


    def get_overlaps_histogram(self):

        fname = "counts_%i.json" % (self.n)
        if not os.path.isfile(fname):

            # self.calculate_all_conformations() # all confs are encoded as a list of strings like 'i+i+i-i-'

            self.trajectories =  [self.fun(entry, '') for entry in self.dict]
            print(self.trajectories)

            # print('all confs generated', len(self.trajectories))
            # print("first conf",  next(self.trajectories))
            # print('there are: %i confs'%len(Overlap.keep_result.all))

            print('encoding to coords')
            self.encoded_conformations = self.encode_to_coords()
            print('done encoding to coords')
            self.get_overlaps()
            # self.get_number_of_contacts()


        else:
            dct = open(fname, 'r').read()
            dct = json.loads(dct)
            self.overlaps_hist = dict(zip([int(el) for el in dct.keys()], dct.values()))

        return self.overlaps_hist



    def get_number_of_contacts(self, c):
        """
        """

        nc = 0
        for i in range(len(c)):
                contacts = [(c[i][0] + 1, c[i][1], c[i][2]), \
                            (c[i][0] - 1, c[i][1], c[i][2]), \
                            (c[i][0], c[i][1] + 1, c[i][2]), \
                            (c[i][0], c[i][1] - 1, c[i][2]), \
                            (c[i][0], c[i][1], c[i][2] + 1), \
                            (c[i][0], c[i][1], c[i][2] - 1) \
                            ]

                nc += len([cn for cn in contacts if cn in c]) - 2
        if nc%2 != 0:
            print(nc, c)
            sys.exit()

        return  nc/2




    def save_overlaps_histogram(self):
        fname = "counts_%i.json" % (self.n)
        if not os.path.isfile(fname):
            json.dump(self.overlaps_hist, open(fname, 'w'))

        fname = "contacts_%i.json" % (self.n)
        if not os.path.isfile(fname):
            json.dump(self.numbers_of_contact, open(fname, 'w'))

    def plot_overlaps_histogram(self):

        self.overlaps_hist = dict(zip(self.overlaps_hist.keys(), [v/sum(self.overlaps_hist.values()) for v in self.overlaps_hist.values()]))
        plt.bar(self.overlaps_hist.keys(), self.overlaps_hist.values())
        plt.yscale('log')
        plt.xlabel('number of overlaps')
        plt.ylabel('number of conformations')



if __name__=="__main__":

    overlaps = Overlap(16)
    print(overlaps)
    print('calculating trajectories')

    print(overlaps.get_overlaps_histogram())
    print(overlaps.numbers_of_contact)

    overlaps.save_overlaps_histogram()

    overlaps.plot_overlaps_histogram()
