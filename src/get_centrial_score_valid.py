import sys
import csv


def get_GIC_lncRNAs(GIC_score):
    fGIC = open(GIC_score,'r')
    GIC_lnc = []
    for line in fGIC:
        data = line.split(',')
        if data[0] == 'lncRNA':
            continue
        if data[0] not in GIC_lnc:
            GIC_lnc.append(data[0])
    print(len(GIC_lnc))
    return GIC_lnc

def get_centrality_score_valid(filename_all, filename_valid):
    f_all = open(filename_all,'r')
    f_valid = open(filename_valid,'w')
    lnc = []
    for line in f_all:
        data = line.split(',')
        if data[0] == 'lncRNA':
            f_valid.write(line)
            continue
        if data[0] in GIC_lnc:
            if data[0] not in lnc:
                lnc.append(data[0])
                f_valid.write(line)
    f_valid.close()


if __name__ == '__main__':
    
    dataPath = '../result/'+ sys.argv[1] + '/'
    savePath = '../result/'+ sys.argv[1] + '/'

    GIC_lnc = get_GIC_lncRNAs(dataPath+'GIC_score.csv')
    DC_score = get_centrality_score_valid(dataPath+'DC_score_allLncRNAs.csv',dataPath+'DC_score.csv')
    BC_score = get_centrality_score_valid(dataPath+'BC_score_allLncRNAs.csv',dataPath+'BC_score.csv')
    CC_score = get_centrality_score_valid(dataPath+'CC_score_allLncRNAs.csv',dataPath+'CC_score.csv')
    EC_score = get_centrality_score_valid(dataPath+'EC_score_allLncRNAs.csv',dataPath+'EC_score.csv')
