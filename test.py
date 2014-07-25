# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 18:49:19 2014

@author: amne51ac

All content licensed under GPL unless otherwise noted or required.
"""
from tabulate import tabulate
from math import *
temp = []
membership = []
with open ('Members.csv', 'r') as myfile:
    membership_headings = myfile.readline().strip().split(',')
    for i in range(len(myfile)-1):
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
            print i, j
            newmember[membership_headings[j]] = float(temp[i][j])
    membership.append(newmember)

with open('Memberoutput.txt', 'w') as myfile:
    myfile.write(tabulate([reversed(i.values()) for i in membership],
                          reversed(membership[1].keys())))

'''print degrees(acos((sin(radians(dec1))*
                    sin(radians(dec2)))+
                   (cos(radians(dec1))*
                    cos(radians(dec2))*
                    cos(radians(ra1)-
                        radians(ra2)))))'''