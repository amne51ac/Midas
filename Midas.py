# -*- coding: utf-8 -*-
"""
Created on Sat Jul 12 09:47:07 2014

@author: amne51ac

All content licensed under GPL unless otherwise noted or required.

GPL: http://www.gnu.org/copyleft/gpl.html
"""
'''
from tkfiledialog import askopenfilename

class Midas:
        
    def __init__(self):
        #i don't really know what needs to be initialized but okay
    
    def import(self):
        #use params to understand the format of the import file
        filename = askopenfilename()
        with open(filename) as myfiel
'''

from math import log

class Midas():

    values = []    

    #def __init__(self):    

    def __init__(self):
        newmidas = []
        count = 0
        midas = [line.strip().split(',') for line in open("Midas Raw\
 Data.csv", 'r')]
        for i in midas[1:]:
            newmidas.append({})
            for index, j in enumerate(i):
                if midas[0][index]:
                    if 'ID' in midas[0][index]:
                        newmidas[count]['ID'] = int(j)
                    else:
                        newmidas[count][midas[0][index]] = float(j)
                else:
                    newmidas[count]['unk' + str(index)] = float(j)
                    
            count += 1
        self.values = newmidas

    def __absolute_mag(self, distance_pc=470):
        for i in range(len(self.values)):
            self.values[i]['mv'] = (self.values[i]['V'] -
                                    5*log((distance_pc/10), 10))

    def __b_minus_v(self):
        for i in range(len(self.values)):
            self.values[i]['bv'] = self.values[i]['B'] - self.values[i]['V']
            '''except:
                print 'The necessary values for the b-v operation are missing'
                break'''

    def __expected_b_minus_v(self):
        for i in range(len(self.values)):
            self.values[i]['xbv'] = ((0.0000176008*(self.values[i]['mv']**6))-
                                     (0.0005951859*(self.values[i]['mv']**5))+
                                     (0.0076215956*(self.values[i]['mv']**4))-
                                     (0.0467622819*(self.values[i]['mv']**3))+
                                     (0.1423842899*(self.values[i]['mv']**2))-
                                     (0.0185234255*self.values[i]['mv'])-
                                     0.1413034474)

    def __b_minus_v_deviation(self):
        for i in range(len(self.values)):
            self.values[i]['bvdev'] = (self.values[i]['bv'] -
                                     self.values[i]['xbv'])

    def __binary_expected_b_minus_v(self, offset=0.753):
        for i in range(len(self.values)):
            self.values[i]['bxbv'] = ((0.0000176008*((self.values[i]['mv']+offset)**6))-
                                     (0.0005951859*((self.values[i]['mv']+offset)**5))+
                                     (0.0076215956*((self.values[i]['mv']+offset)**4))-
                                     (0.0467622819*((self.values[i]['mv']+offset)**3))+
                                     (0.1423842899*((self.values[i]['mv']+offset)**2))-
                                     (0.0185234255*(self.values[i]['mv']+offset))-
                                     0.1413034474)
                                     
    def __binary_b_minus_v_deviation(self):
        for i in range(len(self.values)):
            self.values[i]['binbvdev'] = (self.values[i]['bv'] -
                                          self.values[i]['bxbv'])

    def __q_value(self, offset=0.753):
        for i in range(len(self.values)):
            self.values[i]['Q'] = (-self.values[i]['mv']+
                                   (-7.5071*(self.values[i]['bv']**6)+
                                    33.231*(self.values[i]['bv']**5)-
                                    52.543*(self.values[i]['bv']**4)+
                                    35.909*(self.values[i]['bv']**3)-
                                    9.8787*(self.values[i]['bv']**2)+
                                    6.1966*(self.values[i]['bv'])+1.351
                                    ))/offset
                                                
    def analyze(self, distancepc=470, offset=0.753):
        self.__absolute_mag(distancepc)
        self.__b_minus_v()
        self.__expected_b_minus_v()
        self.__b_minus_v_deviation()
        self.__binary_expected_b_minus_v(offset)
        self.__binary_b_minus_v_deviation()
        self.__q_value(offset)
        
    def headers(self):
        return self.values[0].keys()