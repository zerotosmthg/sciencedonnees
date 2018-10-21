# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 00:56:01 2018

@author: charl
"""
import numpy as np
import math
class dist_norm():
    #Standardized distance matrix (for 2 elements)
    def pairw(x,y):
        #diff = np.asarray(np.abs(x)+1) - np.asarray(np.abs(y)+1)
        diff = np.asarray(x) - np.asarray(y)
        sum = np.asarray(x) + np.asarray(y)
        diff = np.abs(diff)
        #H = np.tile(range(1,9),(8,1))
        '''
        print("diff = ", diff)
        print("sum =" ,sum)
        print("1/sum =" ,1/sum)
        '''
        H = np.tile(diff,(len(diff),1))
        #H = np.abs(H) #+1
        np.place(diff, diff==0, [1])
        F = np.repeat(diff,len(diff)).reshape(len(diff),len(diff))
        #F = np.abs(F+1)
        #G = H/np.abs(diff)
        G = (H/F) /(1/sum)
        '''
        print("F = " ,F)
        print("H = " ,H)
        print("G = " ,G)
        print("Norm G = ", np.linalg.norm(G))
        '''
        return math.sqrt(np.linalg.norm(G))
    
#    rescale_dist(X2[0],X2[1])
#    dist.pairwise((X2[0],X2[1]))
    
    
    #Standardized distance matrix (for a matrix)
    def norm_distance(M):
        D = np.zeros((len(M),len(M)))
        #print("D is: ", D.shape)
        print("Processing...") 
        for j in range(M.shape[0]):
            for i in range(j):
                #print(M[i,:],M[j,:])
                D[i,j] = dist_norm.pairw(M[i,:],M[j,:])
                D[j,i] = D[i,j]
                print(("Processed %d percent of %d."
                       %(round(100*0.5*(j**2-j+2*i+2)/((len(M)**2 - len(M))/2))
                       , (len(M)**2 - len(M))/2)),   end='\r')
#                print("Processed %d of %d." %(0.5*(j**2-j+2*i+2), (len(M)**2 - len(M))/2))
#                print("Processed %d of %d." %(0.5*(j**2-j+2*i+2), (len(M)**2 - len(M))/2))
        #print("Processing has completed!")        
        return D


