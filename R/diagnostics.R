### =========================================================================
### Filter Diagnostics
### -------------------------------------------------------------------------
###
### To diagnose filter performance, several pieces of information
### would be helpful:
### * Results for each filter (as evalSeparately would return, but convenient)
### * Statistics considered by each filter. How to modularize this?
###   - Filter object has a function that returns the stats
###   - Filter function returns stats as an attribute
###   - Should wait to see whether we ever need this
### * Concordance with some ground truth, including FN/TN rates, +/- exclusive
###
### There may also be some helpful plots.

