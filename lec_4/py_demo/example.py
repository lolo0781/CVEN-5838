#! /usr/bin/env python

def myExample():
    print
    print "------------------------------------"
    print "         PYTHON EXAMPLE             "
    print "------------------------------------"
    print

    A = {'fruit':'apple', 'veggie':'asparagus', 'DOB':'120492'}
    print "A = ", A

    for a in A:
        print "A[", a, "] = ", A[a]

    

if __name__=="__main__":

    myExample()

