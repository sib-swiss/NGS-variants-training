genotype_likelihood <- function(m,g,e,ref,alt){
  (((m-g)*e+g*(1-e))^alt * ((m-g)*(1-e)+g*e)^ref)/(m^(ref+alt))
}

# For g = 0 (i.e. 0 reference alleles)
-10*log10(genotype_likelihood(m = 2, g= 0, e = 0.01, ref = 22, alt = 4))
# [1] 80.96026 
-10*log10(genotype_likelihood(m = 2, g= 1, e = 0.01, ref = 22, alt = 4))
# [1] 78.2678
-10*log10(genotype_likelihood(m = 2, g= 2, e = 0.01, ref = 22, alt = 4))
# [1] 440.1746
