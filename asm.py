#!/usr/bin/env python3
import numpy as np
from csvtoballots import *
from collections import deque

def smith_from_losses(losses,cands):
    sum_losses = losses.sum(axis=1)
    min_losses = sum_losses.min()

    # Initialize Smith set and queue
    smith = set(np.compress(sum_losses==min_losses,cands))
    queue = deque(smith)

    # Loop until queue is empty
    while (len(queue) > 0):
        # pop first item on queue
        c = queue.popleft()
        beats_c = np.compress(losses[c],cands)
        # for each candidate who defeats current candidate in Smith,
        # add that candidate to Smith and stick them on the end of the
        # queue
        for d in beats_c:
            if d not in smith:
                smith.add(d)
                queue.append(d)
    return(smith)

def asm(ballots, weight, cutoff=None):
    """T Sorted Margins, using Chris Benham's IBIFA to determine
    the natural approval threshold.
    """

    nballots, ncands = np.shape(ballots)
    cands = np.arange(ncands)
    tw = weight.sum()

    maxscore = int(ballots.max())
    maxscorep1 = maxscore + 1

    if cutoff == None:
        cutoff = maxscore // 2

    # ----------------------------------------------------------------------
    # Tabulation:
    # ----------------------------------------------------------------------
    # A: pairwise array, equal-rated-none
    # B: Same as A, but only including ratings above cutoff
    # T:  Total approval for candidate X
    A = np.zeros((ncands,ncands))
    B = np.zeros((ncands,ncands))
    T = np.zeros((ncands))
    for ballot, w in zip(ballots,weight):
        for r in range(maxscore,0,-1):
            A += np.multiply.outer(np.where(ballot==r,w,0),
                                   np.where(ballot<r ,1,0))
        for r in range(maxscore,cutoff,-1):
            B += np.multiply.outer(np.where(ballot==r,w,0),
                                   np.where(ballot<r ,1,0))
            T += np.where(ballot==r,w,0)

    # ----------------------------------------------------------------------
    # Rank and rating calculations:
    # ----------------------------------------------------------------------
    # Find the Smith set (the set of candidates that defeats all candidates
    # outside the set)
    smith = smith_from_losses(np.where(A.T > A, 1, 0),cands)
    nsmith = len(smith)
    bsmith = smith_from_losses(np.where(B.T > B, 1, 0),cands)
    nbsmith = len(bsmith)

    # Determine Approval Sorted Margins winner:
    # Seed the list by descending approval:
    ranking = np.array([c for c in T.argsort()[::-1] if c in smith])
    branking = np.array([c for c in T.argsort()[::-1] if c in bsmith])

    # Loop until no pairs are out of order pairwise
    maxdiff = T.max() - T.min() + 1
    while True:
        apprsort = T[ranking]
        apprdiff = []
        outoforder = []
        mindiffval = maxdiff
        mindiff = ncands
        for i in range(1,nsmith):
            im1 = i - 1
            c_i = ranking[i]
            c_im1 = ranking[im1]
            if (A[c_i,c_im1] > A[c_im1,c_i]):
                outoforder.append(im1)
                apprdiff.append(apprsort[im1] - apprsort[i])
                if apprdiff[-1] < mindiffval:
                    mindiff = im1
                    mindiffval = apprdiff[-1]
        # terminate when no more pairs are out of order pairwise:
        if (len(outoforder) == 0) or (mindiff == ncands):
            break

        # Do the swap
        ranking[range(mindiff,mindiff+2)] = ranking[range(mindiff+1,mindiff-1,-1)]

    winner = ranking[0]

    return(winner,tw,ncands,maxscore,cutoff,ranking,branking,T,A,B)
 
def test_asm(ballots,weight,cnames,cutoff=None):
    
    winner,tw,ncands,maxscore,cutoff,ranking,branking,T,A,B = asm(ballots,weight,cutoff)

    cands = np.arange(ncands)
    nsmith = len(ranking)
    nbsmith = len(branking)
    approval_ranking = T.argsort()[::-1]

    print("\nFull Pairwise Array:")
    for row in A:
        print(row)

    if (nsmith > 1) and (nbsmith == 1):
        print("\nApproval-only Pairwise Array:")
        for row in B:
            print(row)

    print("\nApproval, using cutoff @ {}:\n".format(cutoff),T)

    print("\nApproval rankings:")
    print("\t{}".format(' > '.join([cnames[c] for c in approval_ranking])))

    print("\nSmith set, ranked by Approval Sorted Margins:")
    print("\t{}".format(' > '.join([cnames[c] for c in ranking])))

    print("\nSmith set for approved ratings only, ranked by approval:")
    print("\t{}".format(' > '.join([cnames[c] for c in branking])))

    if nsmith > 1:
        print("\nApproval Sorted Margins pairwise results:")
        for i in range(1,nsmith):
            im1 = i - 1
            c_i = ranking[i]
            c_im1 = ranking[im1]
            cname_i = cnames[c_i]
            cname_im1 = cnames[c_im1]
            print("\t{}>{}: {} > {}".format(cname_im1,cname_i,A[c_im1,c_i],A[c_i,c_im1]))

    print("\nASM Winner: ", cnames[winner])

    # Condorcet-full//Condorcet-approved//Approval winner
    if nsmith == 1:
        bwinner = ranking[0]
    elif nbsmith == 1:
        bwinner = branking[0]
    else:
        bwinner = approval_ranking[0]

    print("\nC-full//C-approved//Approval Winner: ", cnames[bwinner])

    # Smith//Approval Winner
    if nsmith == 1:
        sawinner = ranking[0]
    else:
        sawinner = [c for c in approval_ranking if c in ranking][0]

    print("\nSmith//Approval Winner: ", cnames[sawinner])

    print("-----\n")
    return

def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-i", "--inputfile",
                        type=str,
                        required=True,
                        help="REQUIRED: CSV Input file")
    parser.add_argument("-c", "--cutoff",
                        type=int,
                        required=False,
                        default=None,
                        help="Optional: Approval cutoff [default: None]")
    args = parser.parse_args()

    ballots, weight, cnames = csvtoballots(args.inputfile)

    # Figure out the width of the weight field, use it to create the format
    ff = '\t{{:{}d}}:'.format(int(log10(weight.max())) + 1)

    print("Ballots:\n\t{}, {}".format('weight',','.join(cnames)))
    for ballot, w in zip(ballots,weight):
        print(ff.format(w),ballot)

    test_asm(ballots, weight, cnames, cutoff=args.cutoff)

if __name__ == "__main__":
    main()
