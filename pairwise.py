#!/usr/bin/env python
#
# Finding Smith Set:
# Start with an empty set (this will hold the solution) and empty queue.
# Find the Copelands winner [candidate with the fewest losses]. (If
# there is a tie, just add all tied candidates.) Add them to the set
# and queue.
# Loop until queue is empty:
# a. Dequeue the first candidate in the queue. Call that candidate C.
# b. Find all the candidates that beat or tie C.
# c. For each such candidate, if they have not already been added to the
#    set, add them to the set and the queue.
# d. Once the queue is empty, the resulting set will be the Smith set.

import numpy as np
from collections import deque
__doc__ = """\
pairwise library --
pairwise.equal_rated_whole:  Count tied ballots as equal whole vote
                             for each candidate against the other
pairwise.equal_rated_none:   Count tied ballots as zero votes
                             for each candidate against the other
pairwise.smith_from_ballots: Return the Smith Set
"""

def equal_rated_whole(ballots, weight):
    "Equal Rated Whole pairwise"
    # Tied votes above the minimum ballot score are counted
    # as whole votes for and against (Equal-Rated Whole = ERW). Note
    # that the pairwise array's resulting diagonal is effectively the
    # total approval score for each candidate
    numcands = np.shape(ballots)[1]
    maxscore = ballots.max()
    totalweight = weight.sum()
    a = np.zeros((numcands,numcands))
    maxscorep1 = maxscore + 1
    for ballot, w in zip(ballots,weight):
        vmin = int(ballot.min() + 1.5)
        for v in range(vmin,maxscorep1):
            a += np.multiply.outer(np.where(ballot==v,w,0),
                                   np.where(ballot<=v,1,0))
            return(totalweight,numcands,a)

def equal_rated_none(ballots, weight):
    # Important: Tied votes above the minimum ballot score are counted
    # as whole votes for and against (Equal-Rated Whole = ERW). Note
    # that the pairwise array's resulting diagonal is effectively the
    # total approval score for each candidate

    numballots, numcands = np.shape(ballots)
    maxscore = ballots.max()
    totalweight = weight.sum()
    A = np.zeros((numcands,numcands))
    TAT = np.zeros((numcands,numcands))
    maxscorep1 = maxscore + 1
    for ballot, w in zip(ballots,weight):
        vmax = int(ballot.max() + 0.5)
        vmin = int(ballot.min() + 1.5)
        for v in range(vmin,maxscorep1):
            A += np.multiply.outer(np.where(ballot==v,w,0),
                                   np.where(ballot<v ,1,0))
        TAT += np.multiply.outer(np.where(ballot==vmax,w,0),
                                 np.where(ballot==vmax,1,0))

    np.fill_diagonal(TAT,0)
    return(totalweight,numcands,A,TAT)

def smith_from_losses(losses):
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

def smith_from_ballots(ballots, weight):
    "Return Smith Set from a bunch of weighted ballots"
    tw, ncand, A = equal_rated_none(ballots, weight)
    cands = np.arange(ncand)
    smith = smith_from_losses(np.where(A.T > A,1,0))
    return(tw,ncands,A,smith)
