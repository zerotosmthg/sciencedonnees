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

        #sum = 1/np.asarray(x) + 1/np.asarray(y)
        u,v = np.asarray(x), np.asarray(y)
        np.place(u, u==0, [1])
        np.place(v, v==0, [1])
        
        sum = 1/u + 1/v
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
        
        #G = (H/F) #/(1/sum)        
        G = H/sum
        '''
        print("Sample F = " ,F[:20,:20])
        print("Sample H = " ,H[:20,:20])
        print("Sample G = " ,G[:20,:20])
        print("Norm G = ", np.linalg.norm(G))
        '''
        return math.sqrt(np.linalg.norm(G))
    
#    rescale_dist(X2[0],X2[1])
#    dist.pairwise((X2[0],X2[1]))
    
    
    #Standardized distance matrix (for a matrix)
    def norm_distance(M):
        D = np.zeros((len(M),len(M)))
        #print("D is: ", D.shape)
        #print("Processing...") 
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
        print('')        
        print("Processing has completed!")        
        return D

'''
S = dist_norm.norm_distance(X2[:200,:])
E = skl.neighbors.DistanceMetric.get_metric('euclidean').pairwise(X2[:200,:])

print("Norm S = %f" %math.sqrt(np.linalg.norm(S)))
print("Norm E = %f" %math.sqrt(np.linalg.norm(E)))

print("Norm Diff = %f" %math.sqrt(np.linalg.norm(E-S)))

#math.sqrt(np.linalg.norm(S-E))

import random
 
m = (0,0,0,0,0,0,0,0)
p = (1,1,1,1,1,1,1,1)
q = (2,2,2,2,2,2,2,2)
s = (3,3,3,3,3,3,3,3)

a = random.randrange(-10,10,step=8)

a,b = np.random.random_integers(-10,10,8)    , np.random.random_integers(-10,10,8)
print(a,b)
dist_norm.pairw(a,b)


'''