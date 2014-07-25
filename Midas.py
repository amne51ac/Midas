# -*- coding: utf-8 -*-
"""
Created on Sat Jul 12 09:47:07 2014

@author: amne51ac

All content licensed under GPL unless otherwise noted or required.

GPL: http://www.gnu.org/copyleft/gpl.html
"""

from math import log
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

class Midas():
    """
Import and append analysis to .csv table of data with headers.

Midas(filename='Raw Midas Data.csv')

Available methods are:

        headers()
            Displays column headers for the currently loaded file.
        
        get_values()
            Returns all data contained within self.__values.
            
        x_y_map()
            Returns a map of all of the stars in the base, color formatted by
            absolute magnitude.
            
        
Appropriate .csv format is (not order restricted, additional columns permitted):

        ID Number, X Position, Y Position, B,     V,     R,     I,     RA,    Declination...
        int,       float,      float...
        
There must be the same number of entries in each record, including header.
Blank cells are permitted in the header, not in data (use 0.0).
0 value not permitted in ID column.
    """
    
    __values = []    

    #def __init__(self):    

    def __init__(self, filename='Midas Raw Data.csv'):

        newmidas = []
        count = 0
        try:
            midas = [line.strip().split(',') for line in open(filename, 'r')]
        except TypeError:
            raise TypeError('Invalid File Format, must be .csv')
        except IOError:
            raise
        try:
            for i in midas[1:]:
                newmidas.append({})
                for index, j in enumerate(i):
                    if midas[0][index]:
                        if 'ID' in midas[0][index]:
                            if int(j) > 0:
                                newmidas[count]['ID'] = int(j)
                            else:
                                newmidas[count]['ID'] = count + 1000000
                        else:
                            newmidas[count][midas[0][index]] = float(j)
                    else:
                        newmidas[count]['unk' + str(index)] = float(j)
                count += 1
        except:
            raise TypeError('Invalid File Format, check headers and data types')
                    
        self.__values = newmidas
        self.__analyze()
        #self.x_y_map()
        

    def __absolute_mag(self, distance_pc=470):
        for i in range(len(self.__values)):
            self.__values[i]['mv'] = (self.__values[i]['V'] -
                                    5*log((distance_pc/10), 10))

    def __b_minus_v(self):
        for i in range(len(self.__values)):
            self.__values[i]['bv'] = self.__values[i]['B'] - self.__values[i]['V']
            '''except:
                print 'The necessary values for the b-v operation are missing'
                break'''

    def __expected_b_minus_v(self):
        for i in range(len(self.__values)):
            self.__values[i]['xbv'] = ((0.0000176008*(self.__values[i]['mv']**6))-
                                     (0.0005951859*(self.__values[i]['mv']**5))+
                                     (0.0076215956*(self.__values[i]['mv']**4))-
                                     (0.0467622819*(self.__values[i]['mv']**3))+
                                     (0.1423842899*(self.__values[i]['mv']**2))-
                                     (0.0185234255*self.__values[i]['mv'])-
                                     0.1413034474)

    def __b_minus_v_deviation(self):
        for i in range(len(self.__values)):
            self.__values[i]['bvdev'] = (self.__values[i]['bv'] -
                                     self.__values[i]['xbv'])

    def __binary_expected_b_minus_v(self, offset=0.753):
        for i in range(len(self.__values)):
            self.__values[i]['bxbv'] = ((0.0000176008*((self.__values[i]['mv']+offset)**6))-
                                     (0.0005951859*((self.__values[i]['mv']+offset)**5))+
                                     (0.0076215956*((self.__values[i]['mv']+offset)**4))-
                                     (0.0467622819*((self.__values[i]['mv']+offset)**3))+
                                     (0.1423842899*((self.__values[i]['mv']+offset)**2))-
                                     (0.0185234255*(self.__values[i]['mv']+offset))-
                                     0.1413034474)
                                     
    def __binary_b_minus_v_deviation(self):
        for i in range(len(self.__values)):
            self.__values[i]['binbvdev'] = (self.__values[i]['bv'] -
                                          self.__values[i]['bxbv'])

    def __q_value(self, offset=0.753):
        for i in range(len(self.__values)):
            self.__values[i]['Q'] = (-self.__values[i]['mv']+
                                   (-7.5071*(self.__values[i]['bv']**6)+
                                    33.231*(self.__values[i]['bv']**5)-
                                    52.543*(self.__values[i]['bv']**4)+
                                    35.909*(self.__values[i]['bv']**3)-
                                    9.8787*(self.__values[i]['bv']**2)+
                                    6.1966*(self.__values[i]['bv'])+1.351
                                    ))/offset
                                    
    def __verify_input_data(self):
        len_check = len(self.__values[1])
        for i in self.__values:
            if type(i) is not dict:
                return False
            if len(i) != len_check:
                return False
            for key,j in i.iteritems():
                if type(j) not in (float, int):
                    return False
                elif type(key) is not str:
                    return False
        return True
                                                
    def __analyze(self, distancepc=470, offset=0.753):
        if self.__verify_input_data():
            self.__absolute_mag(distancepc)
            self.__b_minus_v()
            self.__expected_b_minus_v()
            self.__b_minus_v_deviation()
            self.__binary_expected_b_minus_v(offset)
            self.__binary_b_minus_v_deviation()
            self.__q_value(offset)
        else:
            raise TypeError('Invalid File Format, please ensure each column '+
                            'is properly headed, entry lengths are equal')
    
    def get_values(self):
        return self.__values
    
    def headers(self):
        return self.__values[0].keys()

    def x_y_map(self):

        x = []
        y = []
        c = []
        s = []
        m = []
        
        for i in self.__values:
            x.append(i['X Position'])
            y.append(i['Y Position'])
            m.append(i['mv'])
            if 1 >= i['Q'] >= 0:
                c.append(i['Q'])
            else:
                c.append(0)
            s.append(i['bv'])
            
        for i in c:
            i = i/max(c)
        for i in m:
            i = i/max(m)
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, c=m, edgecolors='None', cmap='OrRd', alpha = .5)
        ax.grid(True)
        fig.tight_layout()
        plt.ion()
        plt.show()

    def import_iso(age=.2):
        iso = []
        with open("ISO.csv") as myfile: 
            iso_presets = myfile.readline().split(',')
            iso_headings = myfile.readline().split(',')
            while True:
                temp = myfile.readline().split(',')
                try:
                    temp = float(temp[1].split()[0])
                except:
                    continue
                if age == temp:
                    while True:
                        temp2 = myfile.readline().split(',')
                        if temp2[1]:
                            iso.append(temp2)
                        else:
                            break
                    break
            
        newiso = []
        
        for i in iso:
            newi = {}
            for j in range(len(i)):
                newi[iso_headings[j]] = float(i[j])
            newiso.append(newi)
        
        isomap = []
        
        for i in newiso:
            isomap.append([i['Mv'], i['B-V']])
        
        return isomap

    def fit_iso(self):
        x = []
        y = []
        isomap = import_iso()
        for i,j in isomap:
            if (i < 12) and (i > 1):
                y.append(i)
                x.append(j)
        
        # fit the data with a 4th degree polynomial
        z4 = polyfit(x, y, 6) 
        p4 = poly1d(z4) # construct the polynomial 
        
        z5 = polyfit(x, y, 11)
        p5 = poly1d(z5)
        
        xx = linspace(-0.1, 1.65)
        pylab.plot(x, y, 'o', xx, p4(xx),'-g', xx, p5(xx),'-b')
        pylab.legend(['Isochrone', '6th degree poly', '11th degree poly'])
        pylab.gca().invert_yaxis()
        pylab.show()