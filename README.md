# biotools
Tools for Biometry and Applied Statistics in Agricultural Science (R package)

Tools designed to perform and work with cluster analysis (including Tocher's algorithm), discriminant analysis and path analysis (standard and under collinearity), as well as some useful miscellaneous tools for dealing with sample size and optimum plot size calculations. Mantel's permutation test can be found in this package. A new approach for calculating its power is implemented. biotools also contains the new tests for genetic covariance components. An approach for predicting spatial gene diversity is implemented.

# Instalation

You can install and load the released version of biotools from CRAN with:

```r
install.packages("biotools")

library(biotools)
```

# Example of Tocher's clustering

```r
data(garlicdist)

(toc <- tocher(garlicdist))

          Tocher's Clustering 

Call: tocher.dist(d = garlicdist)

Cluster algorithm: original 
Number of objects: 17 
Number of clusters: 6 
Most contrasting clusters: cluster 3 and cluster 5, with 
   average intercluster distance: 11.78786

$`cluster 1`
[1] 8  9  12 4  10 2  7  15

$`cluster 2`
[1] 1  6  14

$`cluster 3`
[1] 11 13

$`cluster 4`
[1] 3 5

$`cluster 5`
[1] 16

$`cluster 6`
[1] 17

toc$distClust # cluster distances
          cluster 1 cluster 2 cluster 3 cluster 4 cluster 5 cluster 6
cluster 1  1.745434  4.333530  3.264753  7.070493  8.816863  3.045773
cluster 2  4.333530  1.930265  7.525301  4.156222  3.476651  3.560654
cluster 3  3.264753  7.525301  2.317785  8.019206 11.787861  6.596850
cluster 4  7.070493  4.156222  8.019206  2.324152  4.043741  8.484307
cluster 5  8.816863  3.476651 11.787861  4.043741  0.000000  5.441962
cluster 6  3.045773  3.560654  6.596850  8.484307  5.441962  0.000000

cof <- cophenetic(toc)

cor(cof, garlicdist)
[1] 0.9086886
```
