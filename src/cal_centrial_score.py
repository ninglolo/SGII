import sys
import csv
import numpy as np
import networkx as nx

# Read lncrNA-protein interactions
def readLPI(filename):
    fLPI = open(filename,'r')
    lnc_protein = {}
    
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
    return lnc_protein

    

if __name__ == '__main__':

    dataPath = '../data/'+ sys.argv[1] + '/'
    savePath = '../result/'+ sys.argv[1] + '/'
    
    lnc_protein = readLPI(dataPath+'LPI.csv')
    
    # construct lncRNA-protein-protein interaction network
    g=nx.Graph()
    fnet = open(dataPath+'LPPI.csv', 'r')
    for line in fnet:
        line = line.strip()
        data = line.split(',')
        g.add_edge(data[0], data[1])

    # calculate DC scores
    fDC = open(savePath+'DC_score_allLncRNAs.csv','w')
    fDC.write('lncRNA,score'+'\n')
    den_cn = nx.degree_centrality(g)
    for k in sorted(den_cn,key=den_cn.__getitem__,reverse=True):
        if k in lnc_protein:
            fDC.write(str(k)+','+str(den_cn[k])+'\n')
    fDC.close()
    print('DC finished')

    # calculate BC scores
    fBC = open(savePath+'BC_score_allLncRNAs.csv','w')
    fBC.write('lncRNA,score'+'\n')
    bet_cn = nx.betweenness_centrality (g)
    for k in sorted(bet_cn,key=bet_cn.__getitem__,reverse=True):
        if k in lnc_protein:
            fBC.write(str(k)+','+str(bet_cn[k])+'\n')
    fBC.close()
    print('BC finished')

    # calculate CC scores
    fCC = open(savePath+'CC_score_allLncRNAs.csv','w')
    fCC.write('lncRNA,score'+'\n')
    cen_cn = nx.closeness_centrality(g)
    for k in sorted(cen_cn,key=cen_cn.__getitem__,reverse=True):
        if k in lnc_protein:
            fCC.write(str(k)+','+str(cen_cn[k])+'\n')
    fCC.close()
    print('CC finished')

    # calculate EC scores
    fEC = open(savePath+'EC_score_allLncRNAs.csv','w')
    fEC.write('lncRNA,score'+'\n')
    eig_cen = nx.eigenvector_centrality(g, max_iter=3000)
    for k in sorted(eig_cen,key=eig_cen.__getitem__,reverse=True):
        if k in lnc_protein:
            fEC.write(str(k)+','+str(eig_cen[k])+'\n')
    fEC.close()
    print('EC finished')

