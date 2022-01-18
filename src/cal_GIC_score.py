import os
import sys
import re
import csv
import time
import math
import string
import xlwt

'''
LRMODEL_FEATURES_7 = ['intercept', 'length', 'eng/L', 'cga', 'gcg', 'tcg', 'acg', 'tca']
LRMODEL_COEFS_7 = [0.1625, 2.638e-04, 2.194, 19.88, 37.59, 50.37, 35.44, -64.66] # +/- sample ratio 1:1
LRMODEL_7 = dict(zip(LRMODEL_FEATURES_7, LRMODEL_COEFS_7))
'''

BASES = ['a', 't', 'c', 'g']
TRIPLETS = []
for B1 in BASES:
    for B2 in BASES:
        for B3 in BASES:
            TRIPLETS.append(B1+B2+B3)

mers = {} #先构造一个列表，存放三联体的个数都初始化为0

def readlnc_transID(filename):
    f = open(filename,'r')
    lnc_transID = {}
    ncName = ''
    transID = ''
    for line in f:
        line = line.strip()
        data = line.split(',')
        
        ncName = data[0]
        transID = data[2]
        if ncName == 'ncName':
            continue
        else:
            if ncName not in lnc_transID:
                lnc_transID[ncName] = []
                lnc_transID[ncName].append(transID)
            else:
                lnc_transID[ncName].append(transID)
    return lnc_transID

def readlnc_eng(filename):
    f = open(filename,'r')
    lnc_eng = {}
    for line in f:
        line = line.strip()
        data = line.split(',')
        
        transID = data[0]
        eng = data[1]
        if transID == 'transID':
            continue
        else:
            lnc_eng[transID] = float(eng)
    return lnc_eng

def readfasta(filename):
    f = open(filename,'r')
    res = {}
    for line in f:
        if line.startswith('>'):
            line = line.strip()
            ID = line.split('>',2)[1]
            res[ID] = ''
        else:
            res[ID] += line
    return res

def slidingWindow(seq, l, win, step=1): #滑动窗口
    length = l
    mod = divmod((length-win), step)[1]
    if (win >= length):
        return seq
    else:
        start = 0
        end = win
        fragments = []
        while (len(seq[start:end]) == win):
            fragments.append(seq[start:end])
            start += step
            end += step
        if (mod > 0):
            fragments.append(seq[(length-win):])
        return fragments
 
def stat3mer(seq, l):   #计算频率
    freq = {}
    for item in TRIPLETS:
        mers[item] = 0
    num3mer = float(l-2)
    all3mer = slidingWindow(seq, l, win=3, step=1)
    for i in set(TRIPLETS):
        mers[i] = all3mer.count(i)
    for triplet in TRIPLETS:
        freq[triplet] = mers[triplet]/num3mer
    return freq


if __name__ == '__main__':

    dataPath = '../data/'+ sys.argv[1] + '/'
    savePath = '../result/'+ sys.argv[1] + '/'
    
    lnc_transID = readlnc_transID(dataPath+'ncName_ncID_transID.csv')
    lnc_eng = readlnc_eng(dataPath+'eng.csv')
    trans_seq = readfasta(dataPath+'transcripts_seq.fasta')


    LRMODEL_FEATURES_7 = ['intercept', 'length', 'eng/L', 'cga', 'gcg', 'tcg', 'acg', 'tca']
    if sys.argv[1] == 'mouse':
        LRMODEL_COEFS_7 = [0.1625, 2.638e-04, 2.194, 19.88, 37.59, 50.37, 35.44, -64.66] # +/- sample ratio 1:1
    if sys.argv[1] == 'human':
        LRMODEL_COEFS_7 = [0.7417, 2.612e-04, 4.295, 48.66, 15.64, 76.23, -1.113, -60.29] # +/- sample ratio 1:1
    LRMODEL_7 = dict(zip(LRMODEL_FEATURES_7, LRMODEL_COEFS_7))

           
    features_trans = {}  #计算出每个转录本的特征值
    for k,v in trans_seq.items():
        feature = {}
        seq = v.replace('\n','')
        length = len(seq)
        feature['trans'] = k
        feature['length'] = length
        feature['eng/L'] = lnc_eng[k]/length
        freq = stat3mer(seq, length)
        for item in TRIPLETS:
            feature[item] = freq[item]
        tmp = LRMODEL_7['intercept']+sum([feature[_]*LRMODEL_7[_] for _ in LRMODEL_FEATURES_7[1:]])
        gic = math.exp(tmp)/(math.exp(tmp)+1)
        feature['GIC_7'] = gic
        features_trans[k] = feature

    lncRNA_GIC_score = {}  #计算出每个lnc基因的特征值（一个基因可能包括多个转录本）
    for k,v in lnc_transID.items():
        feature = {}
        feature['length'] = 0
        feature['eng/L'] = 0
        for item in TRIPLETS:
            feature[item] = 0
    
        for trans in v:
            if trans not in features_trans:
                feature = {}
                break
            else:
                tranfeature = features_trans[trans]
                feature['length'] += tranfeature['length']
                feature['eng/L']  += tranfeature['eng/L']
                for item in TRIPLETS:
                    feature[item] += tranfeature[item]
        if len(feature)>0:
            feature['ncName'] = k
            feature['length'] = feature['length']/len(v)
            feature['eng/L'] = feature['eng/L']/len(v)
            for item in TRIPLETS:
                feature[item] =  feature[item]/len(v)
            tmp = LRMODEL_7['intercept']+sum([feature[_]*LRMODEL_7[_] for _ in LRMODEL_FEATURES_7[1:]])
            gic = math.exp(tmp)/(math.exp(tmp)+1)
            feature['GIC_score'] = gic
            lncRNA_GIC_score[k] = feature['GIC_score']

    lncRNA_GIC_score = sorted(lncRNA_GIC_score.items(),key=lambda d:d[1],reverse=True)
    f = open(savePath+'GIC_score.csv', 'w')
    f.write('lncRNA'+','+'score'+'\n')
    for k in lncRNA_GIC_score:
        f.write(k[0] + ',' + str(k[1])+'\n')
    f.close()
