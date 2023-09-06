genotype_likelihood <- function(m, g, e, ref, alt) {
  (((m - g) * e + g * (1 - e))^alt * ((m - g) * (1 - e) + g * e)^ref) / (m^(ref + alt))
}

genotype_likelihood(m = 2, g = 1, e = 0.01, ref = 22, alt = 4)

gl <- sapply(c(0, 1, 2), function(x) {
  genotype_likelihood(
    m = 2,
    g = x,
    e = 0.01,
    ref = 22,
    alt = 4
  )
})

-10 * log10(gl)
