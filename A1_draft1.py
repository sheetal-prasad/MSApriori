#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import itertools
import re,os
transact_item=[]
sdc = 0.0 #deault SDC , If SDC not explicitly provided

sc_count = {}
tail_count = {}
item_pos = {}
rest =0

def main(*args):
    global sc_count, tail_count
    print("enter absolute path if file is present in differnt folder")
    input_filepath = input("Enter the data files's path: ")
    transaction_file = os.path.abspath(input_filepath)
    
    param_filepath = input("Enter the patramter files's path: ")
    param_file = os.path.abspath(param_filepath)
    
    output_filepath = input("Enter the output files's path: ")
    output_file_name = os.path.abspath(output_filepath)
    
    
    T = [] #Transaction items list
    F = [] #FREQUENT ITEM SET
    MIS = {} #MIS DICTIONARY
    
    file_read(transaction_file, param_file, T, MIS)
    F, sc_count, tail_count =  MsApriori(T, MIS)
    
    write_output(output_file_name, F)

def file_read(transaction_file, param_file, transaction_list, mis_dict):
    global sdc

    with open(transaction_file) as t_file:
        for line in t_file:
            transaction = re.sub('[\{\}\s]', '', line).split(',')
            """
            forming a list of unique items from the transaction file,
            used to map MIS in next step
            """
            for ele in transaction:
                if ele not in transact_item:
                    transact_item.append(ele)
            transaction_list.append(transaction)
    
    with open(param_file) as p_file:
        line = p_file.readline()
        while line is not '':
            if 'MIS' in line:
                temp = line.split('=')
                temp[0] = re.sub('MIS|\(|\)|\s', '', temp[0])
                if temp[0]!="rest":
                    mis_dict[temp[0]] = float(temp[1].strip())
                else:
                    rest = float(temp[1].strip())
            if 'SDC' in line:
                sdc = float(line.split('=')[1].strip())
            line = p_file.readline()

    """
    If MIS is not individually present for item, map it to MIS(rest) """  
    for ele in transact_item:
        if ele not in mis_dict:
            mis_dict[ele]= rest
            
def init_pass(M, T, MIS):
    """
    step 2 : Produce the seeds L for generating candidate itemsets of length 2
    """
    global sc_count, item_pos
    n = len(T) 
    L = []
    default_sup_count = 0
    num_of_transactions = n
    num_of_items = len(M)
    
    for item in M:
        for transaction in T:
            if item in transaction:
                sc_count[item] = sc_count.get(item, default_sup_count) + 1

    i = 0
    while i < num_of_items:
        if M[i] in sc_count.keys() and sc_count[M[i]] / num_of_transactions >= MIS[M[i]]:
            L.append(M[i])
            break
        i = i + 1

    for j in range(i + 1, num_of_items):
        if M[j] in sc_count.keys() and sc_count[M[j]] / num_of_transactions >= MIS[M[i]]:
            L.append(M[j])

    for i in range(len(L)):
        item_pos[L[i]] = i
    print("item_pos",item_pos)

    return L

def level2_candidate_gen(L, T, MIS, sup_count):
    """
        Function to the get the level 2 Candidates
    """
    global sdc
    n = len(T) 
    C2 = []
    
    for l in range(len(L) - 1):
        if sup_count[L[l]] / n >= MIS[L[l]]:
            for h in range(l + 1, len(L)):
                if sup_count[L[h]] / n >= MIS[L[l]] and abs((sup_count[L[h]] / n) - (sup_count[L[l]] / n)) <= sdc:
                    can2 = []
                    can2.append(L[l])
                    can2.append(L[h])
                    C2.append(tuple(can2))
    
    return C2

def MIScandidate_gen(Fk_1, T, MIS, sup_count):
    global sdc, item_pos
    n = len(T) 
    
    Ck = []

    for l in range(len(Fk_1) - 1):
        for h in range(l + 1, len(Fk_1)):
            i = 0
            while Fk_1[l][i] == Fk_1[h][i] and i < len(Fk_1[l]) - 1:
                i = i + 1
            if i == len(Fk_1[l]) - 1 and item_pos[Fk_1[l][i]] < item_pos[Fk_1[h][i]]:
                if abs((sup_count[Fk_1[l][i]] / n) - (sup_count[Fk_1[h][i]] / n)) <= sdc:
                    c = Fk_1[l] + tuple([Fk_1[h][-1]])
                    subsets = set(itertools.combinations(c, len(c) - 1))
                    flag = True
                    for s in subsets:
                        if c[0] in s or MIS[c[0]] == MIS[c[1]]:
                            if tuple(s) not in Fk_1:
                                flag = False
                    if flag == True:
                        Ck.append(c)

    return Ck
    
    
def MsApriori(T, MIS):
    """
    Implementation of MS APriori Algorithm
    """
    global sdc, sc_count, tail_count
    
    L = []
    F = []
    FI = []
    n = len(T) 

    #step 1 : sorting on I according to the MIS value of each item
    M = sorted(MIS.keys(), key = lambda i:(MIS[i], i))
    print("M----",M)
    
    
    L = init_pass(M, T, MIS)
    
    #step 3 : Frequent 1-itemsets (F1)
    if len(L) != 0:
        for i in range(len(L)):
            if sc_count[L[i]] / n >= MIS[L[i]]:
                F.append(L[i])
                FI.append(L[i])

    Fk_1 = F

    k = 2
    while len(F) != 0:
        F = []
        if k == 2:
            Ck = level2_candidate_gen(L, T, MIS, sc_count)
        else:   
            Ck = MIScandidate_gen(Fk_1, T, MIS, sc_count)

        for t in T:
            for c in Ck:
                if set(c).issubset(set(t)):
                    sc_count[tuple(c)] = sc_count.get(tuple(c), 0) + 1
                if set(c[1:]).issubset(set(t)):
                    sc_count[tuple(c[1:])] = sc_count.get(tuple(c[1:]), 0) + 1
                    tail_count[tuple(c)] = tail_count.get(tuple(c), 0) + 1

        
        for c in Ck:
            if c in sc_count.keys() and sc_count[tuple(c)] / n >= MIS[c[0]]:
                F.append(c)
                FI.append(c)

        Fk_1 = F
        k = k + 1

    return FI, sc_count, tail_count

def write_output(output_file_name, F):
    """
    formatting and writing to the specified output file/location.
    """
    global sc_count, tail_count

    if len(F) > 0:
        FI = {}
        for i in F:
            if type(i) == type(''):
                if 1 not in FI.keys():
                    t = []
                    t.append(i)
                    FI[1] = list(t)
                else:
                    FI[1].append(i)
            else: 
                if len(i) in FI.keys():
                    FI[len(i)].append(i)
                else:
                    t = []
                    t.append(i)
                    FI[len(i)] = list(t)
            
        with open(output_file_name, 'w') as o_file:
            print(FI)
            for n in range(1, len(FI) + 1):
                o_file.write("(Length {}-{}\n".format(n,len(FI[n])))
                for i in range(len(FI[n])):
                    if n == 1:
                        o_file.write("\t({})\n".format(FI[n][i]))
                    
                    else:
                        o_file.write("\t(")
                        for k in range(len(FI[n][i]) - 1):
                            o_file.write("{} ".format(FI[n][i][k]))
                        o_file.write(FI[n][i][-1])
                        o_file.write(")\n")
                o_file.write(")\n")

    elif len(F) == 0:
        print("There are no frequent itemsets\n")
        with open(output_file_name, 'w') as o_file:
            o_file.write("There are no frequent itemsets\n")



if __name__ == "__main__": main()
