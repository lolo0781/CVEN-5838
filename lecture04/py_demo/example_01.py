#! /usr/bin/env python

def myExample():
    print
    print "------------------------------------"
    print "         PYTHON EXAMPLE             "
    print "------------------------------------"
    print

    A = ['apple', 'asparagus', 'DOB']
    
    print "A = ", A
    print 'Number of elements in A : ', len(A)

    for i in range(len(A)):
        print A[i]



    print

    

if __name__=="__main__":

    myExample()
