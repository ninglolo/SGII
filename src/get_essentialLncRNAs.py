import sys
import csv

# obtain the degree of each lncRNA
def getS(filename):
    fLPI = open(filename,'r')
    lnc_protein = {} #preserve the interaction between lncRNA and protein
    s = {} #save the degree of lncRNAs
    
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
    rank = {} #Preserve the rank of lncRNAs
    for line in f:
        line = line.strip()
        data = line.split(',')
        if data[0] == 'lncRNA':
            continue
        else:
            i += 1
            rank[data[0]] = i
    return rank


# obtain SGII prediction result under different parameter combinations
def getEssentialLncRNAs(K, T, z):
    lnc_bigger_z = []  #Preserve the lncRNAs that the degree is equal to and bigger than z.
    lnc_smaller_z = [] #Preserve the lncRNAs that the degree is smaller than z.
    for k, v in s.items():
        if v >= z:
            lnc_bigger_z.append(k)
        else:
            lnc_smaller_z.append(k)

    essenLnc_centrality = []  #essential lncRNAs predicted by centralities
    for i in lnc_bigger_z:
        if i in r_b and r_b[i] < len(r_b)*K*0.01 and r_c[i] < len(r_c)*K*0.01 and r_d[i] < len(r_d)*K*0.01 and r_e[i] < len(r_e)*K*0.01:
            essenLnc_centrality.append(i)

    essenLnc_GIC = [] #essential lncRNAs predicted by GIC
    for i in lnc_smaller_z:
        if i in r_g and r_g[i] < len(r_b)*T*0.01:
            essenLnc_GIC.append(i)

    essenLnc = [] #essential lncRNAs predicted by SGII
    for i in essenLnc_centrality:
        essenLnc.append(i)
    for i in essenLnc_GIC:
        essenLnc.append(i)

    return essenLnc


if __name__ == '__main__':
    
    dataPath1 = '../data/'+ sys.argv[1] + '/' #path 1 of source files
    dataPath2 = '../result/'+ sys.argv[1] + '/' #path 2 of source files
    savePath = '../result/'+ sys.argv[1] + '/' #save path of result files

    s = getS(dataPath1+'LPI.csv')
    r_b =  getRank(dataPath2+'BC_score.csv')
    r_c =  getRank(dataPath2+'CC_score.csv')
    r_d =  getRank(dataPath2+'DC_score.csv')
    r_e =  getRank(dataPath2+'EC_score.csv')
    r_g = getRank(dataPath2+'GIC_score.csv')

    resultSet = getEssentialLncRNAs(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))
                                                         
    fnew = open(savePath+'essentialLncRNAs.csv', 'w')
    fnew.write('lncRNA\n')
    for i in resultSet:
        fnew.write(i+'\n')
    fnew.close()

