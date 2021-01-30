# companova

R implementation of COMPANOVA

The purpose of companova is to execute an ANOVA on the
results of an experiment with comparative judgements as
proposed by Scheffe, H.  Journal of the statistical
association 1952, 47, 381-400: An Analysis Of Variance For
Paired Comparisons

Usage: companova::run(<input>)
    where "input" refers to a plaintext file, having 8 columns separated
    by spaces or tabs. Each row corresponds to a judgement. For an odd
    row where object O(i) is compared with object O(j) is followed by a
    row where O(j) is compared with O(i). The first cell of a row contains
    the identifier i (or j) while the other 7 contain counts of thelisteners judgements on a seven
    point scale.

For more details, see Rietveld, Toni (2021) Subjective Measurement Techniques for Speech and Language Pathology. London: Routledge