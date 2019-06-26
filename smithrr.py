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


# Check for qualifying, then if any pass, add them to the ranking
# As qualifying candidates are found, add their (r,s,t) ratings
# to the ratings dict
def checkQ(r,s,TG,MXvec,remaining,ratings,ranking):
    Q = []
    for c in remaining:
        t = TG[c]
        mx = MXvec[c]
        if (mx == 0) or (t > mx):
            Q.append(c)
            ratings[c] = (r,s,t)

    if len(Q) > 0:
        for q in Q:
            remaining.discard(q)
        ranking += sorted(Q,key=(lambda c:TG[c]),reverse=True)

def smithrr(ballots, weight):
    """Relevant Rating is an analog of Majority Judgment based on Chris
    Benham's IBIFA method, but does not pass LaterNoHelp because rating
    thresholds depend on scores given to other candidates.

    The method satisfies FBC like Majority Judgment, but is
    independent of irrelevant ballots; the addition of ballots with
    scores for only non-relevant candidates should not affect the winner's
    relevant rating.
    
    Like MJ, RR can fail participation when
    new ballots rate any relevant candidate higher than the previous
    winner's Relevant Rating.

    """

    nballots, ncands = np.shape(ballots)
    cands = np.arange(ncands)
    tw = weight.sum()

    maxscore = ballots.max()
    maxscorep1 = maxscore + 1

    # A: pairwise array, equal-rated-none
    # TEQ:  where X is equal-rated top with Y
    # BEQ:  where X is equal-rated bottom with Y
    # X:  for score r = 0..maxscore, save approval for X when score for Y is below r
    # T:  Total approval at score r = 0..maxscore for candidate X, except T[0] is zero
    A = np.zeros((ncands,ncands))
    TEQ = np.zeros((ncands,ncands))
    BEQ = np.zeros((ncands,ncands))
    X = np.zeros((maxscorep1,ncands,ncands))
    T = np.zeros((maxscorep1,ncands))

    # ----------------------------------------------------------------------
    # Tabulation:
    # ----------------------------------------------------------------------
    for ballot, w in zip(ballots,weight):
        rmin = int(ballot.min() + 0.5)
        rmax = int(ballot.max() + 0.5)
 
        TEQ += np.multiply.outer(np.where(ballot==rmax,w,0),
                                 np.where(ballot==rmax,1,0))

        BEQ += np.multiply.outer(np.where(ballot==rmin,w,0),
                                 np.where(ballot==rmin,1,0))

        for r in range(rmax,rmin,-1):
            A += np.multiply.outer(np.where(ballot==r,w,0),
                                   np.where(ballot<r ,1,0))

            X[r] += w * np.multiply.outer(np.where(ballot<r,1,0),
                                          np.where(ballot>0,1,0))

            T[r] += w * np.where(ballot==r,1,0)

    for r in range(1,maxscorep1):
        np.fill_diagonal(X[r],0)

    np.fill_diagonal(TEQ,0)
    np.fill_diagonal(BEQ,0)
    # ----------------------------------------------------------------------
    # Rank and rating calculations:
    # ----------------------------------------------------------------------
    # Condorcet and Improved Condorcet calculations:
    # Find the Smith set (the set of candidates that defeats all candidates
    # outside the set), the set of candidates that are undefeated when
    # Equal-Top scores are used, and the set of candidates that are undefeated
    # when both equal-top and equal-bottom scores are used symmetrically.

    # Calculate Smith, TEQ undefeated, symmetric undefeated:
    smith = smith_from_losses(np.where(A.T > A, 1, 0),cands)
    
    lossTEQ = np.where(A.T > (A + TEQ), 1, 0)
    losssym = np.where((A.T + BEQ) > (A + TEQ), 1, 0)

    TEQsum = lossTEQ.sum(axis=1)
    minTEQ = lossTEQ.min()
    symsum = losssym.sum(axis=1)
    minsym = symsum.min()

    undefeated_TEQ = set(np.compress(TEQsum==minTEQ,cands))
    undefeated_sym = set(np.compress(symsum==minsym,cands))

    if (len(undefeated_TEQ) == 0):
        undefeated_TEQ = set(cands)

    if (len(undefeated_sym) == 0):
        undefeated_sym = set(cands)

    # ----------------------------------------------------------------------
    # Relevant Ratings:
    # ----------------------------------------------------------------------
    # MX[r,i] is the maximum approval for any candidate on ballots that
    # rate candidate below rating r
    MX = X.max(axis=2)

    # ... and MXC is the candidate with that maximum excluded approval
    MXC = X.argmax(axis=2)

    # Run a series of rounds from MAXSCORE down to zero.
    # At each round, see which candidates have met their relative rating,
    # then add them to the ordered list, "ranking".

    r = maxscore
    TG = np.zeros((ncands))
    ratings = dict()
    ranking = []
    remaining = set(cands)

    # Ratings significance:  
    # 
    # ratings[c] = (r,s,t)
    #   
    # In each triple (r,s,t) for candidate c:
    # 
    # When r is greater than 0,
    # r is the rating level at which 
    # "t", the total number of ballots rating candidate c at level s and above
    # exceeds the maximum approval for any candidate rating c below r.

    # If r equals zero, then candidate c was not able to pass any qualifying threshold 
    # and was sorted at the end of the ranking by total approval.

    while (len(remaining) > 0):
        TG += T[r]            # Note that the T[0] list is all zeros, by construction above.
        s = r                 # "s" is used to save the ratings level used in TG
        checkQ(r,s,TG,MX[r],remaining,ratings,ranking)

        if len(remaining) > 0:
            r -= 1
            checkQ(r,s,TG,MX[r],remaining,ratings,ranking)
            # Note that this terminates when r == 0, since MX[0] is a list of zeros
            # by construction

    # For Benham's IBIFA method, just comment out the second checkQ in the while loop
    # ----------------------------------------------------------------------
    # Use Relevant Ratings to sort the Condorcet sets:
    # ----------------------------------------------------------------------
    # invert the rankings list:  rrindex[i] is the placement in the rankings for candidate i
    rrindex = [r for r in ranking]
    for i, r in enumerate(ranking):
        rrindex[r] = i

    rr_smith = sorted(list(smith),key=(lambda c:rrindex[c]))
    rr_TEQ = sorted(list(undefeated_TEQ),key=(lambda c:rrindex[c]))
    rr_sym = sorted(list(undefeated_sym),key=(lambda c:rrindex[c]))

    return(tw,ncands,maxscore,ratings,ranking,rr_smith,rr_TEQ,rr_sym,A,MX,MXC)

def test_smithrr(ballots,weight,cnames):
    
    tw,ncands,maxscore,ratings,ranking,rr_smith,rr_TEQ,rr_sym,A,MX,MXC = smithrr(ballots,weight)

    print("\nPairwise Array:\n", A)

    print("\nRelevant Ratings:\n")
    for c in ranking:
        r, s, t = ratings[c]
        mx = MX[r,c]
        mxc = MXC[r,c]
        print("\t{}: ({}, {}, {}), {}: {}".format(cnames[c],r,s,t,cnames[mxc],mx))

    print("\nRelevant Ratings ranking:\n")
    print("\t{}".format(' > '.join([cnames[c] for c in ranking])))
        
    print("\nSmith//RR:\n")
    print("\t{}".format(' > '.join([cnames[c] for c in rr_smith])))

    print("\nTEQ//RR:\n")
    print("\t{}".format(' > '.join([cnames[c] for c in rr_TEQ])))

    print("\nSym//RR:\n")
    print("\t{}".format(' > '.join([cnames[c] for c in rr_sym])))

    print("\n-----\n")
    return

def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-i", "--inputfile",
                        type=str,
                        required=True,
                        help="REQUIRED: CSV Input file [default: none]")
    args = parser.parse_args()

    ballots, weight, cnames = csvtoballots(args.inputfile)

    # Figure out the width of the weight field, use it to create the format
    ff = '\t{{:{}d}}:'.format(int(log10(weight.max())) + 1)

    print("Ballots:\n\t{}, {}".format('weight',','.join(cnames)))
    for ballot, w in zip(ballots,weight):
        print(ff.format(w),ballot)

    test_smithrr(ballots, weight, cnames)

if __name__ == "__main__":
    main()
