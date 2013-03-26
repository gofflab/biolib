'''
Created on Sep 2, 2010

@author: lgoff
'''

#Misc Tools and Utilities
def uniqify(seq): 
    # Not order preserving 
    keys = {} 
    for e in seq: 
        keys[e] = 1 
    return keys.keys()