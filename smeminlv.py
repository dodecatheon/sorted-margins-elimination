#!/usr/bin/env python
# % matplotlib notebook
# import matplotlib
# % config InlineBackend.figure_format = 'svg'
# % matplotlib inline
# import matplotlib.pyplot as plt
import numpy as np
from math import *
import csv
from pprint import pprint
from csvtoballots import csvtoballots

def sme_minlv(A,cands,cnames=[]):
    "Sorted Margins Elimination, MinLV (erw)"
    # erw = equal rated whole:
    # tied approval votes are counted as equal votes for each candidate
    # Enter a cnames list of the same length as cands to turn on verbosity

    # Returns candidates, sorted by method
    ncands = len(cands)
    verbose = (len(cnames) > 0)

    # Terminate recursion at 1 or 2 candidates:
    if ncands == 1:
        return(cands)

    if ncands == 2:
        a, b = cands
        # Ties preserve input order
        if A[b][a] > A[a][b]:
            return(cands[::-1])
        else:
            return(cands)

    # AA == pairwise sub-array, permuted by candidate input order
    AA = A[np.ix_(cands,cands)]

    # track the losing votes locations
    losing_votes = np.where(AA.T > AA, AA, 0)
    # winning_votes = np.where(AA.T < AA, AA, 0)
    tied_votes = np.where(AA.T == AA, AA, 0)
    wintie_votes = np.where(AA.T <= AA, AA, 0)
    np.fill_diagonal(tied_votes,0)
    np.fill_diagonal(wintie_votes,0)

    if verbose:
        print("----------")
        print("Permuted pairwise scores:\n    ", [cnames[q] for q in cands])
        for i, rowAA in enumerate(AA):
            print(cnames[cands[i]], rowAA)
        print("----------")

    # Set up for a descending sort by losing votes:
    scorelist = []
    for (c,
         rowLV, # rowWV,
         rowTV,
         rowWT,
         rowAA) in zip(cands,
                       losing_votes, # winning_votes,
                       tied_votes,
                       wintie_votes,
                       AA):

        cc = np.array([c])
        cands_minus_c = np.compress(cands!=c,cands)

        if np.count_nonzero(rowWT) == 0:
            # shortcut return if we find a defeated-by-all loser:
            if verbose:
                print("Candidate {} is defeated by all other candidates".format(cnames[c]),
                      [cnames[q] for q in cands_minus_c])
            return(np.concatenate((sme_minlv(A,cands_minus_c,cnames), cc)))

        if np.count_nonzero(rowLV) > 0:
            scorelist.append(np.compress(rowLV>0,rowLV).min())
        else:
            if np.count_nonzero(rowTV) == 0:
                # Terminate recursion if we find a beats-all winner
                if verbose:
                    print("Candidate {} defeats all candidates".format(cnames[c]),
                          [cnames[q] for q in cands_minus_c])
                # Shortcut return
                return(np.concatenate((cc, sme_minlv(A,cands_minus_c,cnames))))
                # Try to continue with winning votes
                # scorelist.append(np.compress(rowWV>0,rowWV).min())
            else:
                if verbose:
                    print("Candidate {} is tied with candidates".format(cnames[c]),
                          [cnames[q] for q in np.compress(rowTV>0,cands)])
                scorelist.append(np.compress(rowWT>0,rowWT).min())
                # return(False)

    # convert scorelist to np.array:
    minlv = np.array(scorelist)

    # position in cands when sorted by minimum losing score in descending order
    arglvsort = minlv.argsort()[::-1]

    if verbose:
        print("Sorted MinLV candidates:")
        pprint([z for z in zip([cnames[cands[q]] for q in arglvsort],
                               minlv[arglvsort])])

    # Loop over the minlv array, looking for out-of-(pairwise-defeat-)order pairs
    while True:
        lvsort = minlv[arglvsort]
        lvdiff = []
        outoforder = []
        for i in range(1,ncands):
            im1 = i - 1
            c_i = arglvsort[i]
            c_im1 = arglvsort[i-1]
            if AA[c_i][c_im1] > AA[c_im1][c_i]:
                outoforder.append(i-1)
                lvdiff.append(lvsort[im1] - lvsort[i])

        # terminate loop when no more pairs are out of order pairwise.
        if len(outoforder) == 0:
            break

        # find the minimum pairwise out of order pair, sorted by minimum losing votes
        mindiff = outoforder[sorted(range(len(outoforder)),
                                    key=lvdiff.__getitem__)[0]]

        if verbose:
            print("Minimum MinLV pairwise out-of-order difference at position",mindiff+1)
            print("Swapping pairwise out-of-order candidates {} and {}".format(cnames[cands[arglvsort[mindiff]]],
                                                                      cnames[cands[arglvsort[mindiff+1]]]))

        # ... and swap their order
        arglvsort[range(mindiff,mindiff+2)] = arglvsort[range(mindiff+1,mindiff-1,-1)]
        if verbose:
            print("New candidate order:", [cnames[q] for q in cands[arglvsort]])

    # We can terminate recursion here if ncands == 3
    if ncands == 3:
        return(cands[arglvsort])

    # Once the array is sorted, eliminate the last candidate and run again
    # on the remaining candidates recursively
    if verbose:
        print("Eliminating lowest-ordered candidate {}".format(cnames[cands[arglvsort[-1]]]))
        print("Running sme_minlv on candidates", [cnames[q] for q in cands[arglvsort[:-1]]])

    return(np.concatenate((sme_minlv(A,
                                     cands[arglvsort[:-1]],
                                     cnames),
                           cands[[arglvsort[-1]]])))

def test_sme(ballots,weight,cnames):
    [m,n] = np.shape(ballots)
    maxscore = ballots.max()
    A = np.zeros((n,n))

    # Important:  Tied votes above zero are counted as whole votes for and against
    # (Equal-Rated Whole = ERW)
    for ballot, w in zip(ballots,weight):
        for k in range(maxscore):
            # v ranges from maxscore down to 1
            v = maxscore - k
            WW = w * np.multiply.outer(np.where(ballot==v,1,0),
                                       np.where(ballot<=v,1,0))
            A += WW

    # Note that A's diagonal is effectively the total approval score for each candi
    cands = np.arange(n)
    winsort = sme_minlv(A,cands,cnames)
    print("SME_MinLV ranking = ",
          ",".join([cnames[winsort[k]] for k in range(n)]))

if __name__ == "__main__":
    import sys

    if (len(sys.argv) == 2):
        fname = sys.argv[1]
    else:
        fname = input("Enter csv filename: ")

    ballots, weight, cnames = csvtoballots(fname)

    # Figure out the width of the weight field, use it to create the format
    ff = '{{:{}d}}:'.format(int(log10(weight.max())) + 1)

    print("{}, {}".format('weight',','.join(cnames)))
    for ballot, w in zip(ballots,weight):
        print(ff.format(w),ballot)

    test_sme(ballots,weight,cnames)
