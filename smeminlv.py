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

__doc__ = """\
smeminlv --
Sorted Margins Elimination, Minimum Losing Votes (equal-rated whole) --
Returns ordered ranking with Condorcet (beats-all) winner as 0-th entry if one exists.
"""

def pairwise_erw(ballots, weight):
    # Important: Tied votes above the minimum ballot score are counted
    # as whole votes for and against (Equal-Rated Whole = ERW). Note
    # that the pairwise array's resulting diagonal is effectively the
    # total approval score for each candidate

    numballots, numcands = np.shape(ballots)
    maxscore = ballots.max()
    totalweight = weight.sum()
    A = np.zeros((numcands,numcands))
    maxscorep1 = maxscore + 1
    for ballot, w in zip(ballots,weight):
        vmin = int(ballot.min() + 1.5)
        for v in range(vmin,maxscorep1):
            A += np.multiply.outer(np.where(ballot==v,w,0),
                                   np.where(ballot<=v,1,0))
    return(totalweight,numcands,A)

def mycompress(arr):
    return(np.compress(arr>=0,arr))

def sme_minlv(ballots, weight, cands,cnames=[],scalar=False):
    "Sorted Margins Elimination, MinLV (erw)"
    # erw = equal rated whole:
    # tied approval votes are counted as equal votes for each candidate
    # Enter a cnames list of the same length as cands to turn on verbosity

    # Returns candidates, sorted by method
    ncands = len(cands)
    verbose = (len(cnames) > 0)

    # Terminate recursion at 1 or 2 candidates:
    if ncands == 1:
        if scalar:
            return(cands[0])
        else:
            return(cands)

    totalweight, numcands, AA = pairwise_erw(ballots[:,cands],weight)

    if ncands == 2:
        # Ties preserve input order
        if AA[1,0] > AA[0,1]:
            if scalar:
                return(cands[-1])
            else:
                return(cands[::-1])
        else:
            if scalar:
                return(cands[0])
            else:
                return(cands)

    # track the losing votes locations
    losing_votes = np.where(AA.T > AA, AA, -1)
    np.fill_diagonal(losing_votes,-1)

    tied_votes = np.where(AA.T == AA, AA, -1)
    np.fill_diagonal(tied_votes,-1)

    wintie_votes = np.where(AA.T <= AA, AA, -1)
    np.fill_diagonal(wintie_votes,-1)

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
                       losing_votes,
                       tied_votes,
                       wintie_votes,
                       AA):

        cc = np.array([c])
        cands_minus_c = np.compress(cands!=c,cands)

        WTcompress = mycompress(rowWT)

        if len(WTcompress) == 0:
            # shortcut return if we find a defeated-by-all loser:
            if verbose:
                print("Candidate {} is defeated by all other candidates".format(cnames[c]),
                      [cnames[q] for q in cands_minus_c])
            if scalar:
                return(sme_minlv(ballots, weight, cands_minus_c, cnames, scalar))
            else:
                return(np.concatenate((sme_minlv(ballots, weight, cands_minus_c, cnames), cc)))

        LVcompress = mycompress(rowLV)

        if len(LVcompress) > 0:
            LVmin = LVcompress.min()
            scorelist.append(LVmin)
            AAmin = rowAA.min()
            if (verbose and ( LVmin != AAmin )):
                print(("*** Candidate {} has different "
                       "Minimum Pairwise Support ({}) than MinLV ({})").format(cnames[c],
                                                                               AAmin,
                                                                               LVmin))
        else:
            if len(mycompress(rowTV)) == 0:
                # Terminate recursion if we find a beats-all winner
                if verbose:
                    print("Candidate {} defeats all candidates".format(cnames[c]),
                          [cnames[q] for q in cands_minus_c])
                # Shortcut return
                if scalar:
                    return(c)
                else:
                    return(np.concatenate((cc, sme_minlv(ballots, weight, cands_minus_c, cnames))))
            else:
                if verbose:
                    print("Candidate {} is tied with candidates".format(cnames[c]),
                          [cnames[q] for q in np.compress(rowTV>=0,cands)])
                scorelist.append(WTcompress.min())
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
        mindiffval = lvsort[0] - lvsort[ncands-1]
        mindiff = ncands
        # Loop backwards through the sorted array.
        # The mindiff calculation then ensures that tied differences
        # are resolved with preference to lower ranked pairs
        for i in range(ncands-1,0,-1):
            im1 = i - 1
            c_i = arglvsort[i]
            c_im1 = arglvsort[i-1]
            if AA[c_i,c_im1] > AA[c_im1,c_i]:
                outoforder.append(im1)
                lvdiff.append(lvsort[im1] - lvsort[i])
                if lvdiff[-1] < mindiffval:
                    mindiff = im1
                    mindiffval = lvdiff[-1]

        # terminate loop when no more pairs are out of order pairwise.
        if (len(outoforder) == 0) or (mindiff == ncands):
            break

        # find the minimum pairwise out of order pair, sorted by minimum losing votes
        # mindiff = outoforder[sorted(range(len(outoforder)),key=lvdiff.__getitem__)[0]]

        if verbose:
            print(("Swapping candidates {} and {} "
                   "at minimum pairwise out-of-order "
                   "MinLV difference").format(cnames[cands[arglvsort[mindiff]]],
                                              cnames[cands[arglvsort[mindiff+1]]]))

        # ... and swap their order
        arglvsort[range(mindiff,mindiff+2)] = arglvsort[range(mindiff+1,mindiff-1,-1)]
        if verbose:
            print("New candidate order:", [cnames[q] for q in cands[arglvsort]])

    # We can terminate recursion here if ncands == 3
    if ncands == 3:
        if scalar:
            return(cands[arglvsort[0]])
        else:
            return(cands[arglvsort])

    # Once the array is sorted, eliminate the last candidate and run again
    # on the remaining candidates recursively
    if verbose:
        print("Eliminating lowest-ordered candidate {}".format(cnames[cands[arglvsort[-1]]]))
        print("Running sme_minlv on candidates", [cnames[q] for q in cands[arglvsort[:-1]]])

    if scalar:
        return(sme_minlv(ballots,weight,cands[arglvsort[:-1]],cnames,scalar))
    else:
        return(np.concatenate((sme_minlv(ballots,
                                         weight,
                                         cands[arglvsort[:-1]],
                                         cnames),
                               cands[[arglvsort[-1]]])))

def test_sme(ballots,weight,cnames,scalar=False):
    numcands = np.shape(ballots)[1]
    cands = np.arange(numcands)
    if scalar:
        print("SME_MinLV winner =", cnames[sme_minlv(ballots,
                                                     weight,
                                                     cands,
                                                     cnames,
                                                     scalar=scalar)])
    else:
        winsort = sme_minlv(ballots,weight,cands,cnames)
        print("SME_MinLV ranking = ",
              " > ".join([cnames[winsort[k]] for k in range(numcands)]))

        # Test whether LIIA is violated (i.e., if SME_MinLV with winner excluded
        # is the same as the runner up from the full contest)

        if numcands > 2:
            winner, runnerup = winsort[0:2]

            cands_minus_winner = np.compress(cands!=winner,cands)

            newsort = sme_minlv(ballots,weight,cands_minus_winner)
            newwinner = newsort[0]

            if (newwinner != runnerup):
                print(("*** Potential LIIA violation found ***\n"
                       "Winner with winner {} excluded: {}; "
                       "original runner up: {}").format(cnames[winner],
                                                        cnames[newwinner],
                                                        cnames[runnerup]))

def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-i", "--inputfile",
                        type=str,
                        required=True,
                        help="REQUIRED: CSV Input file [default: none]")
    parser.add_argument("-w", "--winneronly",
                        action='store_true',
                        default=False,
                        help="Calculate winner-only, no ranking [default: False]")

    args = parser.parse_args()

    ballots, weight, cnames = csvtoballots(args.inputfile)

    # Figure out the width of the weight field, use it to create the format
    ff = '{{:{}d}}:'.format(int(log10(weight.max())) + 1)

    print("{}, {}".format('weight',','.join(cnames)))
    for ballot, w in zip(ballots,weight):
        print(ff.format(w),ballot)

    test_sme(ballots,weight,cnames,scalar=args.winneronly)

if __name__ == "__main__":
    main()
