#-------------------------------------------------------------------------------
# Name:        Ipop, Moran's I adjusted by population
# Purpose:
#
# Author:      Liang-Huan Chin
#
# Created:     10/09/2013
# Copyright:   (c) Leo Chin 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
'''
This script is trying to implement Global Moran's I adjusted by population
N. Oden. 1995. Adjusting Moran's I for population density.
            Statistics in Medicine, 14:17-26.
Python package is required in this script:
    numpy: https://pypi.python.org/pypi/numpy
'''
import numpy as NUM



if __name__ == '__main__':
    # read infected polygon ID and population info
    pop_csv = open("Polygon_Population.csv")
    pop = NUM.array([row.split(',') for row in pop_csv], int)
    # pop[blockgroups_id, population]
    total_pop = NUM.sum(pop[:,1])
    pop_csv.close()
##    print pop[0,1]
    # read infected_polygon ID and its neighbors
    neigh_csv = open("polygon_relation.csv")
    neighbor = NUM.array([row.split(',') for row in neigh_csv], int)
    # neighbor[blockgroups_id, neighbor_id]
    neigh_csv.close()
##    print neighbor[0]
    # read every_date and total_infections in each every_date
    inf_csv = open("date_infection.csv")
    infection = NUM.array([row.split(',') for row in inf_csv], int)
    # infection[date, infection numbers]
    inf_csv.close()
##    print infection
    # read every_date and total_infections in each every_date
    dateInf_csv = open("date_inf_in_blockgp.csv")
    dateInf = NUM.array([row.split(',') for row in dateInf_csv], int)
    # dateInf[date, blockgroups_id, infection numbers]
    dateInf_csv.close()

##    # write a r, p recording csv file
##    rpfile = open("rp.csv", 'w')
##
##    for i in xrange(len(infection)):
##        if infection[i,1] > 0:
##            for j in xrange(len(dateInf)):
##                if infection[i,0] == dateInf[j,0]:
##                    # ri is the observed proportions of all cases falling in region i
##                    ri = float(dateInf[j,2]) / infection[i,1]
##                    for n in xrange(len(pop)):
##                        if dateInf[j,1] == pop[n,0]:
##                            # pi is the expected proportions of all cases falling in region i
##                            pi = float(pop[n,1]) / total_pop
##                            for k in xrange(len(neighbor)):
##                                if neighbor[k,0] == dateInf[j,1]:
##                                    itemindex = NUM.where(pop[:,0] == neighbor[k,1])
##                                    # pj is the expected proportions of all cases falling in neighbor regions j
##                                    pj = float(pop[itemindex,1]) / total_pop
####                                    print pop[itemindex,1]
##                                    itemindex = NUM.where((dateInf[:,1] == neighbor[k,1]) & (dateInf[:,0] == infection[i,0]))
##                                    # rj is the observed proportions of all cases falling in neighbor regions j
##                                    rj = float(dateInf[itemindex,2]) / infection[i,1]
##                                    rpfile.write(str(infection[i,0])+",")
##                                    rpfile.write(str(ri)+",")
##                                    rpfile.write(str(rj)+",")
##                                    rpfile.write(str(pi)+",")
##                                    rpfile.write(str(pj)+",")
##                                    rpfile.write("\n")
##    rp.close()

    inputfile = open("rp.csv", 'r')
    rp = NUM.array([row.split(',') for row in inputfile])
    # rp[blockgroups_id, population]
    date = rp[:,0].astype(int)
    ri = rp[:,1].astype(float)
    rj = rp[:,2].astype(float)
    pi = rp[:,3].astype(float)
    pj = rp[:,4].astype(float)
    inputfile.close()
    [A,B,C,E,E2, sigma, sigmari, sigmapi, S0, n, temp, temp2] = NUM.repeat(0.00, 12)
    Ipop=[]
    V=[]
    for i in xrange(len(infection)):
        n = infection[i,1]
        exp = float(-1) / total_pop
        if infection[i,1] >0:
            for j in xrange(len(date)):
                if j > 0:
                    if date[j] == infection[i,0]:
                        A=A+(pi[j]*pj[j])
                        C=C+(pi[j]*pj[j]*4)
                        E2=E2+((pj[j]*2)*(pj[j]*2))
                        sigma=sigma+((ri[j]-pi[j])*(rj[j]-pj[j]))
                        if pi[j]!=pi[j-1]:
                            B=B+pi[j]
                            E=E+(pi[j]*E2)
                            sigmari=sigmari+ri[j]
                            sigmapi=sigmapi+pi[j]
                        else:
                            B=pi[j]
                            E=pi[j]*E2
                            sigmari=ri[j]
                            sigmapi=pi[j]
                else:
                    A=pi[0]*pi[0]
                    B=pi[0]
                    C=pi[0]*pj[0]*4
                x = total_pop
                S0=(x*x*A)-(x*B)
                #print str(b_bar)
            temp=((n*n*sigma)-(n*(1-(2*(n/x)))*sigmari)-(n*(n/x)*sigmapi))/(S0*n/x*(1-(n/x)))
            Ipop.append(float(temp))
            temp2=((2*A*A)+(C/2)-E)/(A*A*x*x)
            V.append(float(temp2))
            [A,B,C,E,E2, sigma, sigmari, sigmapi, S0, n, temp, temp2] = NUM.repeat(0.0, 12)
        else:
            Ipop.append(-9999)
            V.append(-9999)

    outputfile=open("Ipop.csv",'w')
    outputfile.write("date"+",")
    outputfile.write("Ipop"+",")
    outputfile.write("E[Ipop]"+",")
    outputfile.write("V[Ipop]"+",")
    outputfile.write("\n")

    for i in range(len(Ipop)):
        outputfile.write(str(infection[i,0])+",")
        outputfile.write(str(Ipop[i])+",")
        outputfile.write(str(exp)+",")
        outputfile.write(str(V[i])+",")
        outputfile.write("\n")

    outputfile.close()



