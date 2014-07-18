# -*- coding: utf-8 -*-
"""
Created on Sat Jul 12 09:47:07 2014

@author: amne51ac
"""

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
                    newmidas[count][index] = float(j)
                    
            count += 1
        self.values = newmidas

    def bv(self):
        for i in range(len(self.values)):
            self.values[i]['BV'] = self.values[i]['B'] - self.values[i]['V']
            '''except:
                print 'The necessary values for the b-v operation are missing'
                break'''
    
    def mv(self, distance_pc):
        for i in range(len(self.values)):
            self.values[i]['MV'] = (self.values[i]['V'] -
                                    5*log((distance_pc/10), 10))
    
    def xmv(self):
        for i in range(len(self.values)):
            self.values[i]['xMV'] = ((0.0000176008*(self.values[i]['BV']**6))-
                                    (0.0005951859*(self.values[i]['BV']**5))+
                                    (0.0076215956*(self.values[i]['BV']**4))-
                                    (0.0467622819*(self.values[i]['BV']**3))+
                                    (0.1423842899*(self.values[i]['BV']**2))-
                                    (0.0185234255*self.values[i]['BV'])-
                                    0.1413034474)
    
    
