#! /usr/bin/env python

def myExample():
    print
    print "------------------------------------"
    print "         PYTHON EXAMPLE             "
    print "------------------------------------"
    print

    A = ['apple', 'asparagus', 'cherry']
    B = {'fav fruit':'mango', 'fav veggie':'brussel sprout', 'DOB':'120492'}
    
    return A, B

if __name__=="__main__":

    A, B = myExample()

    print A
    print 
    print B
