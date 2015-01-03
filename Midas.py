# -*- coding: utf-8 -*-
"""
Created on Sat Jul 12 09:47:07 2014

@author: amne51ac

All content licensed under GPL unless otherwise noted or required.

GPL: http://www.gnu.org/copyleft/gpl.html
"""

from tabulate import tabulate
from math import log, acos, sin, cos, radians, asin, atan
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
#import pylab

class Midas():
    
    """
Import and append analysis to .csv table of data with headers.

Midas(filename='Raw Midas Data.csv')

This program is not intended to be saving previous sessions, which is why input
is accepted as CSV while output is in a readable tabular format for checking
for errors and consistency only.

Available methods are:

        headers()
            Displays column headers for the currently loaded file.
        
        get_values()
            Returns all data contained within self.__values.
            
        x_y_map()
            Returns a map of all of the stars in the base, color formatted by
            absolute magnitude.
            
        save_it()
            Saves a tabular txt copy of the cluster data for reviewing only, 
            can not be re-imported.
            
        
Appropriate .csv format is (not order restricted, additional columns permitted):

        ID Number, X Position, Y Position, B,     V,     R,     I,     RA,    Declination...
        int,       float,      float...
        
There must be the same number of entries in each record, including header.
Blank cells are permitted in the header, not in data (use 0.0).
0 value not permitted in ID column.
    """
    
    __values = []
    __iso_values = []
    __membership = []

    def __init__(self, filename='Midas Raw Data.csv', distancepc=470, offset=0.753):
        self.__import_data(filename)
        if self.__verify_input_data():
            self.__absolute_mag(distancepc)
            self.__b_minus_v()
            self.__expected_b_minus_v()
            self.__b_minus_v_deviation()
            self.__binary_expected_b_minus_v(offset)
            self.__binary_b_minus_v_deviation()
            self.__q_value(offset)
            self.__add_member_mate()
            self.__import_members()
            self.__b1950_j2000()
        else:
            raise TypeError('Invalid File Format, please ensure each column '+
                            'is properly headed, entry lengths are equal')
            

    def __import_data(self, filename='Midas Raw Data.csv'):

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

    def __absolute_mag(self, distance_pc=470):
        for i in range(len(self.__values)):
            self.__values[i]['mv'] = (self.__values[i]['V'] -
                                    5*log((distance_pc/10), 10))

    def __b_minus_v(self):
        for i in range(len(self.__values)):
            self.__values[i]['bv'] = (self.__values[i]['B'] -
                                      self.__values[i]['V'])
            '''except:
                print 'The necessary values for the b-v operation are missing'
                break'''

    def __expected_b_minus_v(self):
        fit = self.__fit_iso_xbv()
        for i in range(len(self.__values)):
            self.__values[i]['xbv'] = ((fit[0]*(self.__values[i]['mv']**11))+
                                       (fit[1]*(self.__values[i]['mv']**10))+
                                       (fit[2]*(self.__values[i]['mv']**9))+
                                       (fit[3]*(self.__values[i]['mv']**8))+
                                       (fit[4]*(self.__values[i]['mv']**7))+
                                       (fit[5]*(self.__values[i]['mv']**6))+
                                       (fit[6]*(self.__values[i]['mv']**5))+
                                       (fit[7]*(self.__values[i]['mv']**4))+
                                       (fit[8]*(self.__values[i]['mv']**3))+
                                       (fit[9]*(self.__values[i]['mv']**2))+
                                       (fit[10]*self.__values[i]['mv'])+
                                        fit[11])

    def __b_minus_v_deviation(self):
        for i in range(len(self.__values)):
            self.__values[i]['bvdev'] = (self.__values[i]['bv'] -
                                     self.__values[i]['xbv'])

    def __binary_expected_b_minus_v(self, offset=0.753):
        fit = self.__fit_iso_xbv()
        for i in range(len(self.__values)):
            self.__values[i]['bxbv'] = ((fit[0]*((self.__values[i]['mv']+offset)**11))+
                                        (fit[1]*((self.__values[i]['mv']+offset)**10))+
                                        (fit[2]*((self.__values[i]['mv']+offset)**9))+
                                        (fit[3]*((self.__values[i]['mv']+offset)**8))+
                                        (fit[4]*((self.__values[i]['mv']+offset)**7))+
                                        (fit[5]*((self.__values[i]['mv']+offset)**6))+
                                        (fit[6]*((self.__values[i]['mv']+offset)**5))+
                                        (fit[7]*((self.__values[i]['mv']+offset)**4))+
                                        (fit[8]*((self.__values[i]['mv']+offset)**3))+
                                        (fit[9]*((self.__values[i]['mv']+offset)**2))+
                                        (fit[10]*(self.__values[i]['mv']+offset))+
                                         fit[11])
                                     
    def __binary_b_minus_v_deviation(self):
        for i in range(len(self.__values)):
            self.__values[i]['binbvdev'] = (self.__values[i]['bv'] -
                                          self.__values[i]['bxbv'])

    def __q_value(self, offset=0.753):
        fit = self.__fit_iso_xmv()
        for i in range(len(self.__values)):
            self.__values[i]['Q'] = (-self.__values[i]['mv']+
                                     ((fit[0]*(self.__values[i]['bv']**11))+
                                      (fit[1]*(self.__values[i]['bv']**10))+
                                      (fit[2]*(self.__values[i]['bv']**9))+
                                      (fit[3]*(self.__values[i]['bv']**8))+
                                      (fit[4]*(self.__values[i]['bv']**7))+
                                      (fit[5]*(self.__values[i]['bv']**6))+
                                      (fit[6]*(self.__values[i]['bv']**5))+
                                      (fit[7]*(self.__values[i]['bv']**4))+
                                      (fit[8]*(self.__values[i]['bv']**3))+
                                      (fit[9]*(self.__values[i]['bv']**2))+
                                      (fit[10]*(self.__values[i]['bv']))+
                                       fit[11]))/offset
                                    
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
        
        
    def __add_member_mate(self):
        for i in range(len(self.__values)):
            self.__values[i]['member mate'] = ''
                                                
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

    def hr_diagram(self):

        x = []
        y = []
        
        for i in self.__values:
            if (i['bvdev'] < 1):
                x.append(i['mv'])
                y.append(i['bv'])
            
        fig, ax = plt.subplots()
        ax.scatter(y, x)
        #ax.grid(True) (i['bv'] < 20) and
        #fig.tight_layout()
        #plt.ion()
        plt.gca().invert_yaxis()
        plt.show()

    def __import_iso(self, age=.2):
        iso = []
        with open("ISO.csv") as myfile: 
            iso_headings = myfile.readline().split(',')
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
        x = []
        y = []
        for i,j in isomap:
            if (i < 12) and (i > 1):
                y.append(i)
                x.append(j)        
        
        return x, y

    '''def fit_iso(self):
        x = []
        y = []
        isomap = self.__import_iso()
        for i,j in isomap:
            if (i < 12) and (i > 1):
                y.append(i)
                x.append(j)
        # fit the data with a 4th degree polynomial
        z4 = np.polyfit(x, y, 6) 
        p4 = np.poly1d(z4) # construct the polynomial
        z5 = np.polyfit(x, y, 11)
        p5 = np.poly1d(z5)
        print z5
        print p5(12)
        
        xx = np.linspace(-0.1, 1.65)
        pylab.plot(x, y, 'o', xx, p4(xx),'-g', xx, p5(xx),'-b')
        pylab.legend(['Isochrone', '6th degree poly', '11th degree poly'])
        pylab.gca().invert_yaxis()
        pylab.show()'''
        
    def __fit_iso_xbv(self, age = .2):
        x, y = self.__import_iso(age)
        return np.polyfit(y, x, 11)
        
    def __fit_iso_xmv(self, age = .2):
        x, y = self.__import_iso(age)
        return np.polyfit(x, y, 11)
        
    def save_it(self, filename = 'Midas_Output.txt'):
        with open(filename, 'w') as myfile:
            myfile.write(tabulate([i.values() for i in self.get_values()], self.headers()))
            
    def __import_members(self, memfilename='Members.csv'):
        temp = []
        membership = []
        with open (memfilename, 'r') as myfile:
            for i, l in enumerate(myfile):
                pass
            member_file_length = i
        with open (memfilename, 'r') as myfile:
            membership_headings = myfile.readline().strip().split(',')
            for i in range(member_file_length):
                temp0 = myfile.readline().split(',')
                temp.append(temp0)
        
        for i in range(len(temp)):
            newmember = {}
            for j in range(len(temp[i])):
                if j == 0:
                    newmember[membership_headings[j]] = int(temp[i][j])
                elif j == 1:#RA1950
                    temp1 = temp[i][j].split()
                    newmember[membership_headings[j]] = 15*(float(temp1[0])+
                                              float(temp1[1])/60+
                                              float(temp1[2])/3600)
                elif j == 2:#DEC1950
                    temp1 = temp[i][j].split()
                    newmember[membership_headings[j]] = (float(temp1[0])+
                                            float(temp1[1])/60+
                                            float(temp1[2])/3600)
                elif j == 11:#RA2
                    temp1 = temp[i][j].split()
                    newmember[membership_headings[j]] = 15*(float(temp1[0])+
                                           float(temp1[1])/60+
                                           float(temp1[2])/3600)
                elif j == 12:#DEC2
                    temp1 = temp[i][j].split()
                    newmember[membership_headings[j]] = (float(temp1[0])+
                                         float(temp1[1])/60+
                                         float(temp1[2])/3600)
                elif j in [8, 9]:
                    newmember[membership_headings[j]] = temp[i][j]
                else:
                    newmember[membership_headings[j]] = float(temp[i][j])
            membership.append(newmember)
        
        self.__membership = membership
    
    def __separation(self, ra1, dec1, ra2, dec2):
        return acos((sin(radians(dec1))*
                     sin(radians(dec2)))+
                    (cos(radians(dec1))*
                     cos(radians(dec2))*
                     cos(radians(ra1)-
                         radians(ra2))))
    
    def __b1950_j2000(self):
        for i in range(len(self.__membership)):
            x = (cos(self.__membership[i]['RA1950'])*
                 cos(self.__membership[i]['DE1950']))
            y = (sin(self.__membership[i]['RA1950'])*
                 cos(self.__membership[i]['DE1950']))
            z = sin(self.__membership[i]['DE1950'])
            X = 0.999926*x-0.011179*y-0.004859*z
            Y = 0.011179*x+0.999938*y-0.000027*z
            Z = 0.004859*x-0.000027*y+0.999988*z
            r = atan(Y/X)
            if X < 0:
                self.__membership[i]['RA2000'] = r + 180
            elif (Y < 0) and (X > 0):
                self.__membership[i]['RA2000'] = r + 360
            else:
                self.__membership[i]['RA2000'] = r
                print "RA2000 ERROR 1"
            self.__membership[i]['DE2000'] = asin(Z)
    
    def mating(self):
        for i,m in enumerate(self.__membership):
            mate = {'best' : 0, 'score' : 99999}
            for j,v in enumerate(self.__values):
                if 'mate' in v:
                    continue
                score = self.__separation(v['RA'],v['Declination '],m['RA2000'],m['DE2000'])
                if score < mate['score']:
                    mate['score'] = score
                    mate['best'] = j
            self.__values[mate['best']]['member mate'] += '%s-%s' %(i, mate['score'])