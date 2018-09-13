# sorted-margins-elimination
# Sorted Margins Elimination
## Minimum Losing Votes, Equal-Rated Whole
### A Condorcet-completion method that is resistant to the Chicken Dilemma

Approval Sorted Margins was introduced by Forest Simmons as a symmetric
modification to Definitive Majority Choice (AKA Ranked Approval Voting).

https://electowiki.miraheze.org/wiki/Approval_Sorted_Margins
https://wiki.electorama.com/wiki/Approval_Sorted_Margins
https://electowiki.miraheze.org/wiki/Marginal_Ranked_Approval_Voting
https://wiki.electorama.com/wiki/Marginal_Ranked_Approval_Voting

Chris Benham chose a different metric for the seed sort, Minimum Losing Vote
(AKA MinLV), with equal-rated non-zero ratings scored as whole votes for each
candidate in the pairwise match-up.

http://lists.electorama.com/pipermail/election-methods-electorama.com/2016-October/000599.html

Benham claims

> This meets Smith, Plurality, Mono-raise, Mono-switch-plump, Non-drastic Defense.

Since the minimum Losing Votes metric can change, depending on what candidates
are included, Benham also introduced a step to eliminate the lowest rated
candidate after symmetric marginal pairwise sorting, which automatically gives
the method Local Independence from Irrelevant Alternatives (LIIA), since
successively eliminating the lowest pairwise sorted candidate ensures that the
higher ranking result will be the same when lower ranked candidates are
removed.

Benham gives three examples, provided here as `examples/benham1.csv`,
`examples/benham2.csv`, and `examples/benham3.csv`, which demonstrate that
attempts to win by burying a more popular similar candidate will fail.

Interestingly, this resistancy to burying/chicken dilemma strategizing reduces
the likelihood of a Favorite Betrayal Criterion failure.  Such failures can be
defeated by raising a lower-rated compromise candidate to equal top rank.
