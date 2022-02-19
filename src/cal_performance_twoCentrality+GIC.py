import sys
import csv
import scipy.stats as stats
import math
from scipy.stats import norm

# obtain the degree of each lncRNA
def getS(filename):
    fLPI = open(filename,'r')
    lnc_protein = {}
    s = {}
    
    for line in fLPI:
        line = line.strip()
        data = line.split(',')
        lncName = data[0]
        proteinID = data[1]
        if lncName == 'lncRNA':
            continue
        if lncName not in lnc_protein:
            lnc_protein[lncName] = []
            lnc_protein[lncName].append(proteinID)
        else:
            lnc_protein[lncName].append(proteinID)

    for k,v in lnc_protein.items():
        s[k] = len(v)
    return s


# obtain the ranking of each lncRNA
def getRank(filename):
    lnc = []
    f = open(filename, 'r')
    i = 0
    rank = {}
    for line in f:
        line = line.strip()
        data = line.split(',')
        if data[0] == 'lncRNA':
            continue
        else:
            i += 1
            rank[data[0]] = i
    return rank


# Obtain the essential lncRNAs predicted under a parameter combination
def getEssentialLncRNAs(K, T, z):
    lnc_bigger_z = []
    lnc_smaller_z = []
    for k, v in s.items():
        if v >= z:
            lnc_bigger_z.append(k)
        else:
            lnc_smaller_z.append(k)

    essenLnc_centrality = []
    for i in lnc_bigger_z:
        if i in r_cen_A and r_cen_A[i] < len(r_cen_A)*K*0.01 and r_cen_B[i] < len(r_cen_B)*K*0.01:
            essenLnc_centrality.append(i)

    essenLnc_GIC = []
    for i in lnc_smaller_z:
        if i in r_g and r_g[i] < len(r_g)*T*0.01:
            essenLnc_GIC.append(i)

    essenLnc = []
    for i in essenLnc_centrality:
        essenLnc.append(i)
    for i in essenLnc_GIC:
        essenLnc.append(i)

    return essenLnc


# obtain known essential lncRNAs
def readess(filename):
    res = []
    f = open(filename,'r')
    for line in f:
        res.append(line.strip())
    return res



if __name__ == '__main__':
    
    if sys.argv[1] == 'mouseHomologousOfHuman':
        dataPath1 = '../data/mouse/'
        dataPath2 = '../result/mouse/'
    else:
        dataPath1 = '../data/'+ sys.argv[1] + '/'
        dataPath2 = '../result/'+ sys.argv[1] + '/'
    savePath = '../result/performance/'
    
    s = getS(dataPath1+'LPI.csv')
    r_cen_A = getRank(dataPath2+ sys.argv[2] + '_score.csv')
    r_cen_B = getRank(dataPath2+ sys.argv[3] + '_score.csv')
    r_g = getRank(dataPath2+'GIC_score.csv')


    if sys.argv[1] == 'human': # parameter setting of human
        ess = readess(dataPath1+'esslnc_homo.txt')
        K_arr = [5,10,15,20,25]
        T_arr = [5,10,15,20,25]
        flag_str = ',5,5,5,5,5,10,10,10,10,10,15,15,15,15,15,20,20,20,20,20,25,25,25,25,25\n'
    if sys.argv[1] == 'mouse': # parameter setting of mouse 
        ess = ['Xist','Gas5','Tsix','Meg3','Fendrr','Dnm3os','Gt(ROSA)26Sor','Braveheart']
        K_arr = [1,3,5,7,9]
        T_arr = [1,3,5,7,9]
        flag_str = ',1,1,1,1,1,3,3,3,3,3,5,5,5,5,5,7,7,7,7,7,9,9,9,9,9\n'
    if sys.argv[1] == 'mouseHomologousOfHuman': # parameter setting of homologous of human
        ess = ['Xist','Gas5','Airn','Meg3','Tug1','HOTAIR','Fendrr','Stard13','Neat1','Nron','TINCR']
        K_arr = [1,3,5,7,9]
        T_arr = [1,3,5,7,9]
        flag_str = ',1,1,1,1,1,3,3,3,3,3,5,5,5,5,5,7,7,7,7,7,9,9,9,9,9\n'
    
    N = len(r_cen_A)
    M = len(ess)
    print(N,M)

    # create a file to save performance indicators
    f = open(savePath+sys.argv[1]+'_'+sys.argv[2]+'+'+sys.argv[3]+'+GIC_performance.csv','w')

    # calculate performance indicators and write them to a file
    for z in  [5,10,15,20]:
        f.write('度数：'+str(z)+'\n')
        f.write(flag_str)
        f.write(',X,Y,Sen,FET_p,FPR,X,Y,Sen,FET_p,FPR,X,Y,Sen,FET_p,FPR,X,Y,Sen,FET_p,FPR,X,Y,Sen,FET_p,FPR\n')
        
        for K in K_arr:
            strline = str(K)
            for T in T_arr:
                resultSet = getEssentialLncRNAs(K, T, z)
                n = 0
                for i in ess:
                    if i in resultSet:
                        n += 1
                x = len(resultSet)
                y = n
                n00 = y
                n01 = M-y
                n10 = x-y
                n11 = N-y-(M-y)-(x-y)
                sen = n00/M #sensitivity

                oddsRatio, FET_p = stats.fisher_exact([[n00, n01], [n10, n11]]) # calculate fisher's exact test index
                FET_p = -1 * math.log10(FET_p)

                FPR = round(n10/(N-M),6)

                everystr= str(x)+','+str(y)+','+str(sen)+',' +str(FET_p)+ ','+ str(FPR)
                strline+= ',' + everystr
                
            strline += '\n'
            f.write(strline)
    f.close()
