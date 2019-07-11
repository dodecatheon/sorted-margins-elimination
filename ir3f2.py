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
def checkQ(r,TG,MXvec,remaining,ratings,ranking):
    Q = []
    for c in remaining:
        t = TG[c]
        mx = MXvec[c]
        if (mx == 0) or (t > mx):
            Q.append(c)
            ratings[c] = (r,t)

    if len(Q) > 0:
        for q in Q:
            remaining.discard(q)
        ranking += sorted(Q,key=(lambda c:TG[c]),reverse=True)

def ir3f2(ballots, weight):
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

    # ----------------------------------------------------------------------
    # Tabulation:
    # ----------------------------------------------------------------------
    # X:  for score r = 0..maxscore, save approval for X when score for Y is below r
    # T:  Total approval at score r = 0..maxscore for candidate X, except T[0] is zero
    B = np.zeros((maxscorep1,ncands,ncands))
    X = np.zeros((maxscorep1,ncands,ncands))
    T = np.zeros((maxscorep1,ncands))
    for ballot, w in zip(ballots,weight):
        rmin = int(ballot.min() + 0.5)
        rmax = int(ballot.max() + 0.5)
        for r in range(rmax,rmin,-1):
            B[r] += np.multiply.outer(np.where(ballot==r,w,0),
                                      np.where(ballot<r ,1,0))
            X[r] += np.multiply.outer(np.where(ballot<r,w,0),
                                      np.where(ballot>0,1,0))
            T[r] += np.where(ballot==r,w,0)

    for r in range(1,maxscorep1):
        np.fill_diagonal(X[r],0)


    # A: pairwise array, equal-rated-none
    # B: pairwise array, equal-rated-none, for candidates with rate r
    A = B.sum(axis=0)

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
    # ratings[c] = (r,t)
    #   
    # In each pair (r,t) for candidate c:
    # 
    # When r is greater than 0,
    # r is the rating level at which "t", the total number of ballots
    # rating candidate c at level r and above, exceeds the maximum
    # approval for any candidate rating c below r.

    # If r equals zero, then candidate c was not able to pass any qualifying threshold 
    # and was sorted at the end of the ranking by total approval.

    while (len(remaining) > 0):
        TG += T[r]            # Note that the T[0] list is all zeros, by construction above.
        checkQ(r,TG,MX[r],remaining,ratings,ranking)

        if len(remaining) > 0:
            r -= 1

    # ----------------------------------------------------------------------
    # Use IBIFA (Relevant Ratings) to sort the Condorcet sets:
    # ----------------------------------------------------------------------
    # invert the rankings list:  rrindex[i] is the placement in the rankings for candidate i
    rrindex = [r for r in ranking]
    for i, r in enumerate(ranking):
        rrindex[r] = i

    appr_level = ratings[ranking[0]][0]

    A2 = B[appr_level:].sum(axis=0)

    Approval = T[appr_level:].sum(axis=0)

    # ----------------------------------------------------------------------
    # Rank and rating calculations:
    # ----------------------------------------------------------------------
    # Find the Smith set (the set of candidates that defeats all candidates
    # outside the set)
    fullsmith = smith_from_losses(np.where(A.T > A, 1, 0),cands)
    # Find the Smith set using only inferred Approval ratings
    apprsmith = smith_from_losses(np.where(A2.T > A2, 1, 0), cands)

    if (len(fullsmith) == 1):
        winner = list(fullsmith)[0]
    elif (len(apprsmith) == 1):
        winner = list(apprsmith)[0]
    else:
        winner = sorted(cands,key=(lambda c:Approval[c]),reverse=True)[0]
    
    return(winner,tw,ncands,maxscore,ratings,ranking,list(fullsmith),list(apprsmith),Approval,A,A2,MX,MXC)
 
def test_ir3f2(ballots,weight,cnames):
    
    winner,tw,ncands,maxscore,ratings,ranking,fullsmith,apprsmith,approval,A,A2,MX,MXC = ir3f2(ballots,weight)

    cands = np.arange(ncands)

    print("\nFull Pairwise Array:")
    for row in A:
        print(row)

    print("\nApproved Pairwise Array:")
    for row in A2:
        print(row)

    print("\nApproval:")
    print(approval)

    print("\nIBIFA Ratings:")
    for c in ranking:
        r, t = ratings[c]
        mx = MX[r,c]
        mxc = MXC[r,c]
        print("\t{}: ({}, {}), {}: {}".format(cnames[c],r,t,cnames[mxc],mx))

    print("\nIBIFA ranking:")
    print("\t{}".format('>'.join([cnames[c] for c in ranking])))
        
    print("\nFull Smith:\n")
    print("\t{}".format(','.join([cnames[c] for c in fullsmith])))

    print("\nApproval Smith:\n")
    print("\t{}".format(','.join([cnames[c] for c in apprsmith])))

    print("\nApproval rankings:\n")
    print("\t{}".format('>'.join([cnames[c] for c in sorted(cands,key=(lambda c:approval[c]),reverse=True)])))

    print("\nWinner =", cnames[winner])
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

    test_ir3f2(ballots, weight, cnames)

if __name__ == "__main__":
    main()
