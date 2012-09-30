### =========================================================================
### Binomial Likelihood Ratio (p.lower/p.error) Test
### -------------------------------------------------------------------------
### EXPERIMENTAL/private
### 

setClass("VariantTest")

setClass("BinomialLRTest",
         representation(p.lower = "numeric", p.error = "numeric"),
         contains = "VariantTest")


         
