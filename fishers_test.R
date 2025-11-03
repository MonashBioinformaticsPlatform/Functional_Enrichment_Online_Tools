# https://bioinformatics3.erc.monash.edu/rsconnect/content/241/

# Total differentially expressed genes = 300

# Total background genes = 16,000

# Genes with pathway term = 500

# Overlap (DE genes with the term) = 30

# Then:

a <- 30 #                                    # (differentially expressed genes annotated with the term)

b <- 300 - 30 # = 270 #                        # (differentially expressed genes not annotated with the term)

c <- 500 - 30 # = 470 #                        # (genes in the term but not differentially expressed)

d <- 16000 - ((300 + 500) - 30) # = 15230 #    # (genes neither differentially expressed nor in the term)

data <- matrix(c(30, 270, 470, 15230), nrow = 2, byrow = TRUE)
ft <- fisher.test(data, alternative = "greater")
ft
ft$p.value #  2.047319e-08


# Given 10000 terms being tested, what is the adjusted p-value?
ft$p.value * 10000 <- 0.0002047319
