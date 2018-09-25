# sorted-margins-elimination
# Sorted Margins Elimination
## Minimum Losing Votes, Equal-Rated Whole
### A Condorcet-completion method that is resistant to the Chicken Dilemma

Approval Sorted Margins was introduced by Forest Simmons as a symmetric
modification to Definitive Majority Choice (AKA Ranked Approval Voting).

https://electowiki.org/wiki/Approval_Sorted_Margins

https://wiki.electorama.com/wiki/Approval_Sorted_Margins

https://electowiki.org/wiki/Marginal_Ranked_Approval_Voting

https://wiki.electorama.com/wiki/Marginal_Ranked_Approval_Voting

Chris Benham chose a different metric for the seed sort, Minimum Losing Vote
(AKA MinLV), with equal-rated above-bottom ratings scored as whole votes for each
candidate in the pairwise match-up.

http://lists.electorama.com/pipermail/election-methods-electorama.com/2016-October/000599.html

Benham originally claimed

> This meets Smith, Plurality, Mono-raise, Mono-switch-plump, Non-drastic Defense.

but has since dropped the claim of Mono-switch-plump, as it appears to be
incompatible with the Condorcet criterion.

Benham notes that without elimination, Simmon's sorted margins using Minimum
Losing Votes might not be clone-independent, since the MinLV metric can change
depending on which candidates are included.

`smeminlv.py` can be used to run a single winner election.  Typical usage is

   `smeminlv.py -i examples/benham1.csv`

`qrarv.py` can be used to run a multiwinner PR election, in the style of
Bucklin Reweighted Voting, but with SME used as the method to find the winner
of each seat.  Typical usage is

   `qrarv.py -m 9 -q droop -i examples/june2011.csv`

Both scripts require that the CSV inputfile be specified using the `-i` option.

Benham gives three examples, provided here as `examples/benham1.csv`,
`examples/benham2.csv`, and `examples/benham3.csv` (input files provided in
the `examples` subdirectory of this repository), which demonstrate that
attempts to win by burying a more popular similar candidate will either fail
to stop that candidate from winning, or will backfire by electing the
candidate from the opposition.
