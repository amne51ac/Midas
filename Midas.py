# -*- coding: utf-8 -*-
"""
Created on Sat Jul 12 09:47:07 2014

@author: amne51ac

All content licensed under GPL unless otherwise noted or required.

GPL: http://www.gnu.org/copyleft/gpl.html
"""

from tabulate import tabulate
from math import log, acos, sin, cos, radians
import numpy as np
from numpy import deg2rad, transpose, dot, arcsin, arctan2, zeros, ndarray, array, rad2deg, pi
from premat import premat
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
            
        display_all_jp()
            Displays both a Cartesian and Herzsprung-Russel depiction of all of
            the stars in the Jones-Prosser survey, no filter.  Returns a count
            of all of the stars fitting the above criteria.
        
        display_all_mated_members(mem=0)
            Displays both a Cartesian and Herzsprung-Russel depiction of all of
            the stars that were successfully mated to Midas data points.  Only
            filterable by Jones-Prosser membership probability. Returns a count
            of all of the stars fitting the above criteria.
        
        display_all_mates_membership(bvdev=0.1, lowq=0, highq=1, mem=0)
            Displays both a Cartesian and Herzsprung-Russel depiction of all of
            the stars that were successfully mated to Midas data points.  Can
            be filtered by membership probability, Q value, and B-V deviation.
            Returns a count of all of the stars fitting the above criteria.
        
        hr_diagram(bvdev=0.1)
            Displays a Hertzsprung-Russel diagram of the stars in the Midas 
            (values) dataset.
        
        show_unmated()
            Returns each of the Jones-Prosser stars that are unmated to a Midas
            data point.
            
        
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
        print "Welcome to Midas"
        self.__import_data(filename)
        if self.__verify_input_data():
            self.__absolute_mag(distancepc)
            self.__b_minus_v()
            self.__expected_b_minus_v()
            self.__b_minus_v_deviation()
            self.__binary_expected_b_minus_v(offset)
            self.__binary_b_minus_v_deviation()
            self.__q_value(offset)
            self.__import_members()
            self.__b1950_j2000()
            self.__add_member_mate()
            self.__add_member_match_count()
            self.__distance_and_visual_mating()
        else:
            raise TypeError('Invalid File Format, please ensure each column '+
                            'is properly headed, entry lengths are equal')
            

    def __import_data(self, filename='Midas Raw Data.csv'):
        print "Importing Data"

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
        print "Calculating Absolute Magnitudes"
        for i in range(len(self.__values)):
            self.__values[i]['mv'] = (self.__values[i]['V'] -
                                    5*log((distance_pc/10), 10))

    def __b_minus_v(self):
        print "Calculating B-V Values"
        for i in reversed(range(len(self.__values))):
            if self.__values[i]['B'] < 30:
                self.__values[i]['bv'] = (self.__values[i]['B'] -
                                          self.__values[i]['V'])
            else:
                self.__values.pop(i)
            '''except:
                print 'The necessary values for the b-v operation are missing'
                break'''

    def __expected_b_minus_v(self):
        print "Calculating Expected B-V Values"
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
        print "Finding B-V Deviation"
        for i in range(len(self.__values)):
            self.__values[i]['bvdev'] = (self.__values[i]['bv'] -
                                     self.__values[i]['xbv'])

    def __binary_expected_b_minus_v(self, offset=0.753):
        print "Projecting Expected B-V for Binary"
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
        print "Finding Binary B-V Deviation"
        for i in range(len(self.__values)):
            self.__values[i]['binbvdev'] = (self.__values[i]['bv'] -
                                          self.__values[i]['bxbv'])

    def __q_value(self, offset=0.753):
        print "Calculating Q Value"
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
        print "Verifying Input Data"
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
        print "Adding Member Mate Field"
        for i, k in enumerate(self.__values):
            self.__values[i]['mate_candidates'] = []
            
    def __add_member_match_count(self):
        print "Adding Member Match Count Field"
        for i, k in enumerate(self.__membership):
            self.__membership[i]['match_count'] = 0
                                                
    def get_values(self):
        return self.__values
    
    def headers(self):
        return self.__values[0].keys()

    def x_y_map(self):
        print "Displaying Cartesian Map"
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

    def hr_diagram(self, bvdev=0.1):
        print "Genertaing Hertzsprung-Russel Diagram"
        x = []
        y = []
        
        for i in self.__values:
            if (abs(i['bvdev']) < bvdev):
                x.append(i['mv'])
                y.append(i['bv'])
                
        ix, iy = self.__import_iso()
            
        fig, ax = plt.subplots()
        ax.scatter(y, x)
        plt.plot(ix, iy, 'ro')
        #ax.scatter(ix, iy, edgecolors='None', cmap='ro')
        #ax.grid(True) (i['bv'] < 20) and
        #fig.tight_layout()
        #plt.ion()
        plt.gca().invert_yaxis()
        plt.show()

    def __import_iso(self, age=.2):
        print "Importing Yonsei-Yale Isochrone Data"
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
        print "Generating 11th Degree Polynomial Regression for X-BV"
        x, y = self.__import_iso(age)
        return np.polyfit(y, x, 11)
        
    def __fit_iso_xmv(self, age = .2):
        print "Generating 11th Degree Polynomial Regression for X-MV"
        x, y = self.__import_iso(age)
        return np.polyfit(x, y, 11)
        
    def save_it(self, filename = 'Midas_Output.txt'):
        print "Saving Output"
        with open(filename, 'w') as myfile:
            myfile.write(tabulate([i.values() for i in self.get_values()], self.headers()))

    def __import_members(self, memfilename='Members.csv'):
        print "Importing Jones-Prosser Membership Data"
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
                         
    def __precess(self, ra0, dec0, equinox1, equinox2, doprint=False, fk4=False, radian=False):
        scal = True
        if isinstance(ra0, ndarray):
            ra = ra0.copy()
            dec = dec0.copy()
            scal = False
        else:
            ra = array([ra0])
            dec = array([dec0])
        npts = ra.size
    
        if not radian:
            ra_rad = deg2rad(ra)     # Convert to double precision if not already
            dec_rad = deg2rad(dec)
        else:
            ra_rad = ra
            dec_rad = dec
    
        a = cos(dec_rad)
    
        x = zeros((npts, 3))
        x[:, 0] = a * np.cos(ra_rad)
        x[:, 1] = a * np.sin(ra_rad)
        x[:, 2] = np.sin(dec_rad)
    
        # Use PREMAT function to get precession matrix from Equinox1 to Equinox2
    
        r = premat(equinox1, equinox2, fk4=fk4)
    
        x2 = transpose(dot(transpose(r), transpose(x)))      # rotate to get
        # output direction cosines
    
        ra_rad = zeros(npts) + arctan2(x2[:, 1], x2[:, 0])
        dec_rad = zeros(npts) + arcsin(x2[:, 2])
    
        if not radian:
            ra = rad2deg(ra_rad)
            ra = ra + (ra < 0.) * 360.e0            # RA between 0 and 360 degrees
            dec = rad2deg(dec_rad)
        else:
            ra = ra_rad
            dec = dec_rad
            ra = ra + (ra < 0.) * 2.0e0 * pi
    
        if doprint:
            print 'Equinox (%.2f): %f,%f' % (equinox2, ra, dec)
        if scal:
            ra, dec = ra[0], dec[0]
        return ra, dec
    
    def __b1950_j2000(self):
        print "Precessing"
        for i, k in enumerate(self.__membership):
            a = self.__precess(k['RA1950'], k['DE1950'], 1950, 2000)
            self.__membership[i]['RA'], self.__membership[i]['Declination'] = a
            
    def __distance_mating(self, dist=0.000075):
        print "Mating By Distance"
        for c, d in enumerate(self.__values):
            for i, k in enumerate(self.__membership):
                b = self.__separation(d['RA'], d['Declination '], k['RA'], k['Declination'])
                if b < dist:
                    self.__values[c]['mate_candidates'].append([b, k])
                    self.__membership[i]['match_count'] += 1
        
    def __mate_check(self):
        print "Performing Mate Check"
        count1 = 0
        count2 = 0
        for i in self.__values:
            if i['mate_candidates']:
                count1 += len(i['mate_candidates'])

        for i in self.__membership:
            if i['match_count'] == 0:
                count2 += 1
        
        print "Successfully Mated %i stars out of 650! \n%i unmated data points remaining from Jones-Prosser." %(count1, count2)
        return count1, count2
            
    def __distance_and_visual_mating(self, dist=0.000025, vdev=0.457): #empirical values
        print "Performing Mating Ritual by Distance and Apparent Magnitude"
        for c, d in enumerate(self.__values):
            for i, k in enumerate(self.__membership):
                xdist = self.__separation(d['RA'], d['Declination '], k['RA'], k['Declination'])
                xvdev = abs(d['V'] - k['Vmag'])
                yvdev = d['V'] - k['Vmag']
                if (xdist < dist) and (xvdev < vdev):
                    self.__values[c]['mate_candidates'].append([xdist, xvdev, yvdev, k])
                    self.__membership[i]['match_count'] += 1
        self.__mate_check()

    def display_mates_membership(self, bvdev=0.1, lowq=0, highq=1, mem=0):
        print "Displaying Jones-Prosser Data Mated to Midas Data That Falls Within Membership and Isochrone/Q Value Threshold"
        x = []
        y = []
        c = []
        s = []
        m = []
        x1 = []
        y1 = []
        count = 0
        
        for i in self.__values:
            if i['mate_candidates'] and ((abs(i['bvdev']) < bvdev) or
                                         (i['Q'] <= highq and
                                          i['Q'] > lowq)):
                if i['mate_candidates'][0][3]['Mem'] > mem:
                    x.append(i['X Position'])
                    y.append(i['Y Position'])
                    m.append(i['mv'])
                    x1.append(i['mv'])
                    y1.append(i['bv'])
                    if 1 >= i['Q'] >= 0:
                        c.append(i['Q'])
                    else:
                        c.append(0)
                    s.append(i['bv'])
                    count += 1
            
        print '%i records were accepted' %count
            
        for i in c:
            i = i/max(c)
        for i in m:
            i = i/max(m)
            
        fig, ax = plt.subplots(nrows=2)
        ax1, ax2 = ax.ravel()
        
        ax1.scatter(x, y, c=m, edgecolors='None', cmap='OrRd', alpha = .5)
        ax1.grid(True)
        ax1.set_title('Spatial arrangement')
        ax1.set_xlabel('X Position')
        ax1.set_ylabel('Y Position')
        
        ix, iy = self.__import_iso()
        ax2.scatter(y1, x1)
        ax2.plot(ix, iy, 'y--')
        ax2.invert_yaxis()
        ax2.set_title('Hertzsprung-Russell Diagram')
        ax2.set_ylabel('Mv')
        ax2.set_xlabel('B-V')
        
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle('Midas and Jones-Prosser Mated Data Points Within 0 < Q <= 1 and .1 B-V Deviation and J-P Membership > 0')
        
        #fig.tight_layout()
        #plt.ion()
        #plt.gca().invert_yaxis()
        plt.show()
        
        return count
        
    def display_all_mated_members(self, mem=0):
        print "Displaying All Jones-Prosser Points Mated to Midas Points Lying Within J-P Membership Threshold"
        x = []
        y = []
        c = []
        s = []
        m = []
        x1 = []
        y1 = []
        count = 0
        
        for i in self.__values:
            if i['mate_candidates']:
                if i['mate_candidates'][0][3]['Mem'] > mem:
                    x.append(i['X Position'])
                    y.append(i['Y Position'])
                    m.append(i['mv'])
                    x1.append(i['mv'])
                    y1.append(i['bv'])
                    if 1 >= i['Q'] >= 0:
                        c.append(i['Q'])
                    else:
                        c.append(0)
                    s.append(i['bv'])
                    count += 1
            
        print '%i records were accepted' %count
            
        for i in c:
            i = i/max(c)
        for i in m:
            i = i/max(m)
            
        fig, ax = plt.subplots(nrows=2)
        ax1, ax2 = ax.ravel()
        
        ax1.scatter(x, y, c=m, edgecolors='None', cmap='OrRd', alpha = .5)
        ax1.grid(True)
        ax1.set_title('Spatial arrangement')
        ax1.set_xlabel('X Position')
        ax1.set_ylabel('Y Position')
        
        ix, iy = self.__import_iso()
        ax2.scatter(y1, x1)
        ax2.plot(ix, iy, 'y--')
        ax2.invert_yaxis()
        ax2.set_title('Hertzsprung-Russell Diagram')
        ax2.set_ylabel('Mv')
        ax2.set_xlabel('B-V')
        
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle('All Successfully Mated Jones-Prosser Data Points With J-P Membership > 0')
        
        #fig.tight_layout()
        #plt.ion()
        #plt.gca().invert_yaxis()
        plt.show()
        
        return count
        
    def display_all_jp(self):
        print "Displaying All Jones-Prosser Data Points Mated to Midas Points"
        x = []
        y = []
        c = []
        s = []
        m = []
        x1 = []
        y1 = []
        count = 0
        
        for i in self.__values:
            if i['mate_candidates']:
                x.append(i['X Position'])
                y.append(i['Y Position'])
                m.append(i['mv'])
                x1.append(i['mv'])
                y1.append(i['bv'])
                if 1 >= i['Q'] >= 0:
                    c.append(i['Q'])
                else:
                    c.append(0)
                s.append(i['bv'])
                count += 1
            
        print '%i records were accepted' %count
            
        for i in c:
            i = i/max(c)
        for i in m:
            i = i/max(m)
            
        fig, ax = plt.subplots(nrows=2)
        ax1, ax2 = ax.ravel()
        
        ax1.scatter(x, y, c=m, edgecolors='None', cmap='OrRd', alpha = .5)
        ax1.grid(True)
        ax1.set_title('Spatial arrangement')
        ax1.set_xlabel('X Position')
        ax1.set_ylabel('Y Position')
        
        ix, iy = self.__import_iso()
        ax2.scatter(y1, x1)
        ax2.plot(ix, iy, 'y--')
        ax2.invert_yaxis()
        ax2.set_title('Hertzsprung-Russell Diagram')
        ax2.set_ylabel('Mv')
        ax2.set_xlabel('B-V')
        
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle('All Successfully Mated Jones-Prosser Data Points')
        
        #fig.tight_layout()
        #plt.ion()
        #plt.gca().invert_yaxis()
        plt.show()
        
        return count
        
    def show_unmated(self):
        for i in self.__membership:
            if not i['match_count']:
                print i
    
            
if __name__ == '__main__':
    m = Midas()