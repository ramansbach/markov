# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 13:13:02 2016

@author: Rachael Mansbach
Set of helper functions for analysis and other in python
"""
import numpy as np
#assume row vectors
def avg_by_var(var,val,throwout=5):
    uvar = np.unique(var)
    ureturn = np.inf*np.ones([3,len(uvar)])
    
    for v in range(0,len(uvar)):
        cval = val[var==uvar[v]]
        if len(cval) >= throwout:
            ureturn[0,v] = uvar[v]
            ureturn[1,v] = np.mean(cval)
            ureturn[2,v] = np.std(cval)
    ur = ureturn[ureturn!=np.inf]
    return np.reshape(ur,[3,len(ur)/3])
    

if __name__ == "__main__":
    v1 = np.array([1,1,1,1,1,2,2,2,2,2,3,3])
    v2 = np.array([1,1,1,1,1,4,4,4,4,4,9,9])
    print avg_by_var(v1,v2)
        
        