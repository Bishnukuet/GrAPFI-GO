from pickle import dump, load
from operator import itemgetter
from matplotlib import  pyplot as plt
import pandas as pd
from  scipy.spatial.distance import euclidean
import random
import numpy as np


RootPath=""
#========================================================
# Read the trainng graph saved in two  files
# domains_per_protein.p contains  domains per each proteins
# proteins_per_domain contains proteins  per each domains
#Loading the kidera vectors for the Training proteins
#========================================================
domains_per_protein=load(open("Data/domains_per_protein.p", 'rb'))
proteins_per_domain=load(open("Data/proteins_per_domain.p",'rb'))
#kidera_per_protein=load(open("graph data/kidera_per_protein.p"))
#=======================================================
# Read the ground truth EC labels associated with each
# protein in the pre-build graph.
#======================================================
training_GO=load(open("Data/Training_GO.p", 'rb'))
#================GO Accuracy=====================================
print(len(domains_per_protein))
print(len(proteins_per_domain))
print(len(training_GO))

#==================================================
#
#Comparing GO terms
#==================================================

def precision_go(A,P):
    M=len(set.intersection(set(A), set(P)))
    N=len(set(P))
    return float(M)/float(N)

def recall_go(A,P):
    M=len(set.intersection(set(A), set(P)))
    N=len(set(A))
    return float(M)/float(N)

def score_prediction_go(A,P):
    return precision_go(A,P), recall_go(A,P)
def f_measure_go(pre, rec):
    if pre+rec>0:
        return 2*(pre*rec)/(pre+rec)
    else:
        return 0.0
def score_GO(A,P):
    r,c =score_prediction_go(A,P)
    return r, c, f_measure_go(r,c) 
     


    

#===========================Similarity Functions==================
def dot(A,B): 
    return (sum(a*b for a,b in zip(A,B)))
def cosine_similarity(a,b):
    return dot(a,b) / ( (dot(a,a) **.5) * (dot(b,b) ** .5) )
def eu_distance(A,B):
    return euclidean(A,B)
def link_strength(D1, D2):
    D1=set(D1)
    D2=set(D2)
    match=len(set.intersection( D1, D2))
    union=len(set.union(D1,D2))
    return float(match)/float(union)


#=======================================================
# Label Propagation
#============================================================

#========================================================
# Given a set fo domains, the following function
# finds the  neighbors in the graph who shares the dommain
#========================================================
def find_neighbors(domains ):

    neigbors=set()
    for dom in domains:

        if dom in proteins_per_domain:

            #P= [x[0] for x in proteins_per_domain[dom]]
            P=proteins_per_domain[dom]

            neigbors=set.union(neigbors,set(P))
    #print neigbors
    return neigbors


#===============================================================
# Computation of various neighborhood-based similarity measure
#===============================================================

# Expansion based only direct neighbors
def CN(setA, setB):
    # Common neighbor
    score=0.0
    anb=set.intersection(setA, setB)
    score=len(anb)
    return score 
def JA(setA, setB):
    # Jaccard Similarity
    score=0.0
    aub=set.union(setA, setB)
    anb=set.intersection(setA, setB)
    score=len(anb)/len(aub)
    return score

def PA(setA, setB):
    # Preferential Attachment
    score=0.0
    score=len(setA)*len(setB)
    return score 

def SA(setA, setB):
    # Salton Index
    score=0.0
    anb=set.intersection(setA, setB)
    pa=np.sqrt(len(setA)*len(setB))
    score=len(anb)/pa
    return score 
def SO(setA, setB):
    # Soronsen Index
    score=0.0
    anb=set.intersection(setA, setB)

    score=2.0*len(anb)/(len(setA)+len(setB))
    return score 

def  HPI(setA, setB):
    # HUb promoted Index
    score=0.0
    anb=set.intersection(setA, setB)
    score=len(anb)/min([len(setA), len(setB)])
    return score 
def  HDI(setA, setB):
    # Hub depressed Index
    score=0.0
    anb=set.intersection(setA, setB)
    score=len(anb)/max([len(setA), len(setB)])
    return score

def LLHN(setA, setB):
    #ocal Leicht-Holme-Newman index (LLHN)
    score=0.0
    anb=set.intersection(setA, setB)
    score=len(anb)/(len(setA)*len(setB))

    return score 

# expansion based on common neighbors

def AA(setA, setB):
    # Adamic Adar
    score=0.0
    anb=set.intersection(setA, setB)
    print(len(anb))
    if len(anb)==0:
        return 0.0
    for z in anb:
        D3=domains_per_protein[z]
        N3=find_neighbors(D3)

        if not len(N3)==0:
            score=score+1.0/np.log(len(N3))
    return score 

def RA(setA, setB):
    # Resource Allocation
    score=0.0
    anb=set.intersection(setA, setB)
    #print(len(anb))
    if len(anb)==0:
        return 0.0
    for z in anb:
        D3=domains_per_protein[z]
        N3=find_neighbors(D3)

        if not len(N3)==0:
            score=score+1.0/len(N3)
    return score 




#============================================================================
# Core annotation task is accomplished by this function
#============================================================================
def annotate(node,domains, lth=0.1, hth=1.0):
    #miss_count=0
    
    #domains = Benchmark_Domains[node]#domains_per_protein[node]
    #print domains
    N=find_neighbors(domains)

    #print N
    if len(N)==0:
        print( "No Neigbors found\n")
    #print node, N
    #print 'Neighbors found', len(N)
    label = {}
    for j in N:

        if j==node:
            continue
        if j not in training_GO:
            continue
        # Retrieve the corresponding EC of the 
        ec= training_GO[j]
        if len(ec)==0:
            #print 'No EC'
            continue

        if j in domains_per_protein:
            D2=domains_per_protein[j]
        else:
            continue

        N2=find_neighbors(D2)
        #sim_score=JA(N, N2)
        #sim_score=CN(N, N2)
        
        #sim_score=PA(N, N2)
        #sim_score=SA(N, N2)

        #sim_score=SO(N, N2)
        #sim_score=HPI(N, N2)
        #sim_score=HDI(N, N2)
        #sim_score=LLHN(N, N2)
        #sim_score=RA(N, N2)
        sim_score=0.0


        link_weight =link_strength(domains, D2)  

        #=========================
        # Combining score

        #Filter the neighbor based on their respective link weight
        if link_weight>=hth:
         #   print domains, D2
            continue
        if link_weight <lth:
            continue
        link_weight=link_weight+sim_score
        for e in ec:
            if e in label:
                label[e]=label[e]+link_weight
            else:
                label.update({e:link_weight})
    if len(label)>0:
        label=rankLabels(label)

    else:
        return []

    return label


#============================================================================
# To rank a dictionary by it's value and return a list of keys
#============================================================================

def softmax(x):
        e_x = np.exp(x - np.max(x))  # for computation stability
        return e_x / e_x.sum()
def normalize_weight(T):
    W=max(T.values())*1.0
    #X=list(T.values())
    #W=softmax(X)
    for t in T.keys():
        T[t]=round(T[t]/W,3)
    return T

def rankLabels(G):
    #return list(set(sorted(G, key=G.get, reverse=False)))
    G=normalize_weight(G)
    if G=={}:
        return []
    A=sorted(G.items(), key=itemgetter(1), reverse=True)
    #Max = A[0][1]
    #Min = Max * .50
    #A = [x for x in A if x[1]]
    return A

#===============================================
# function Annotations
#===============================================


#===============================================
# function Annotations
#===============================================

def function_annotations(nodes, domains, outfile="GrAPFI.pred", lth=0.30,hth=1.0, top_k=1):
    
    W=open(outfile,'w')
    for protein in nodes:
        dom=domains[protein]
        labels = annotate(protein,dom, lth=lth, hth=hth)
        predicted_labels = labels
        if len(predicted_labels)>0:
            for label in predicted_labels:
                GO=label[0].split(':')[1]
                GO, asp=GO.split('_')
                W.write(protein+'\t'+GO+'\t'+asp+'\t'+str(label[1])+'\n')
                #print 'Prediction:',protein, ec[0], ec[1], '\n' 
        else:
            W.write(protein+'\n')
            print (protein, " No Prediction","\n")
    W.close()


if __name__ == '__main__':

    Benchmark_Domains=load(open("Data/MetaGO_Benchmark.p", 'rb'))   
    Nodes=Benchmark_Domains.keys()
    print (len(Nodes))
    function_annotations(Nodes, Benchmark_Domains, lth=0.13, hth=2.0, top_k=60)