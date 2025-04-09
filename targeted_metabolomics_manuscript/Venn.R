install.packages("venneuler")     # Install & load venneuler package
library("venneuler")
# Create new plotting page
library("grid")
grid.newpage()
# Draw pairwise venn diagram (3)
plot(venneuler(c(A = 3243,B = 0,C = 26,"A&B" = 0,"B&C" = 0,"A&C" = 0,"A&B&C" = 288)))
# Draw pairwise venn diagram (2)
#plot(venneuler(c(A = 434,B = 0,"A&B" = 198)))
