*   For each target site, encode one purine/one pyramidine base change for every position.  Choose N of these
*   For each target site, encode 2 random base changes.  Do this N times.

# For each candidate
*   Compute the START threshold by adding up the costs of the changed bases
*   Do phred iterations going from 40 + MAX_START down to 1 + MIN_START, with an initial specificity of START for each candidate
*   Record the *relative specificity* for the condidate when the threshold is crossed.


WEAKNESS is the sum of phred scores of changed bases.  It is the brokenness of the initial match.

The SPECIFICITY as previously measured asserted that the target would have to be *at least that weak* a match for a location other than the specified location.  I.e. a SPECIFICITY of 39 indicates that the WEAKNESS of that guide for its second best target was greater than 39, because bowtie found only one match (the original) when allowing 39 phred-points worth of mismatch

A SPECIFICITY of 0, then, says that the WEAKNESS of the guide for the second best target was no greater than its weakness for the intended target.

In the case of deliberately weakened targets, then, the SPECIFICITY should be the WEAKNESS of the second target (the lowest threshold at which both the first and second best targets were permitted) minus the original target's weakness.
