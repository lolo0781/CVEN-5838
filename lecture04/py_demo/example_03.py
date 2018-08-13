#! /usr/bin/env python

def myExample(desiredEntry):
    print
    print "------------------------------------"
    print "         PYTHON EXAMPLE             "
    print "------------------------------------"
    print

    B = {'fav fruit':'mango', 'fav veggie':'brussel sprout', 'DOB':'120492'}
    
    try:
        print B[desiredEntry]
        return B[desiredEntry]
    except:
        print 'Desired entry not in dictionary'
        exit(0)

if __name__=="__main__":

    B1 = myExample('DOB')
    B2 = myExample('DoB')

    print B1
    print 
    print B2
