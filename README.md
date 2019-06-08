# biotools

Tools designed to perform and work with cluster analysis (including Tocher's algorithm), discriminant analysis and path analysis (standard and under collinearity), as well as some useful miscellaneous tools for dealing with sample size and optimum plot size calculations. Mantel's permutation test can be found in this package. A new approach for calculating its power is implemented. biotools also contains the new tests for genetic covariance components. An approach for predicting spatial gene diversity is implemented. Some of the main tools are described in the next sections.

# Importance of variables: Singh's criterion

We have implemented in the function singh() the method proposed by Singh (1981) for determining the importance of variables based on the squared generalized Mahalanobis distance. In his approach, the importance of the $j$-th variable ($j = 1, 2, ..., p$) on the calculation of the distance matrix can be obtained by:

![](https://latex.codecogs.com/gif.latex?S_%7B.j%7D%20%3D%20%5Csum_%7Bi%3D1%7D%5E%7Bn-1%7D%20%5Csum_%7Bi%27%3Ei%7D%5E%7Bn%7D%20%28x_%7Bij%7D%20-%20x_%7Bi%27j%7D%29%20%28%7B%5Cbf%20x%7D_i%20-%20%7B%5Cbf%20x%7D_%7Bi%27%7D%29%5ET%20%5Cmat%7B%5CSigma%7D_%7B.j%7D%5E%7B-1%7D)

\begin{equation}
	S_{.j} = \sum_{i=1}^{n-1} \sum_{i'>i}^{n} (x_{ij} - x_{i'j}) ({\bf x}_i - {\bf x}_{i'})^T \mat{\Sigma}_{.j}^{-1} 
\end{equation}

\noindent where $x_{ij}$ and $x_{i'j}$ are the observation taken at the $i$-th and $i'$-th objects (individuals) for the $j$-th variable, ${\bf x}_i$ is the $p$-variate vector of the $i$-th object and $\mat{\Sigma}_{.j}^{-1}$ is the $j$-th column of the inverse of the covariance matrix.

The output of singh() is basically a matrix containing every $S_{.j}$ and its relative value.

# Tocher's algorithm

biotools contains the method suggested by K.D. Tocher (Rao 1952) for clustering objects, which consists of an optimization method that follows the algorithm:

\begin{enumerate}
	\item (Input) Take a distance matrix, let us say ${\bf d}$, of size $n$.
	\item Define a clustering criterion, i.e., a limit of distance at which to evaluate the inclusion of objects in a cluster in formation. Usually, this criterion is defined as follows: let ${\bf m}$ be a vector of size $n$ containing the lowest distances involving each object, then take $\theta = \max{({\bf m})}$ to be the clustering criterion.
	\item Locate the pair of objects with the lowest distance in ${\bf d}$ to be the starting cluster.
	\item Compute the average distance ($d_{k(j)}$) of the starting cluster, say $k$, to a certain object, say $j$, that shows the lowest pairwise distance with each object in the starting cluster.
	\item Compare the statistic evaluated at the step 4 with $\theta$. Then, add the new object to the starting cluster if $d_{k(j)} \leq \theta$. Otherwise, go to step 3 not considering the objects already clustered and start a new cluster.
\end{enumerate}

The process continues until the last remaining object is evaluated and either included in the last cluster formed or allocated to its own separated cluster. The function {\ttfamily tocher()} performs optimization clustering and returns an object of class 'tocher', which contains the clusters formed, a numeric vector indicating the cluster of each object, a matrix of cluster distances and also the input, a 'dist' object.

# Cluster distances

After obtaining the clusters, it might be useful to know how divergent they are from each other. In this context, cluster distances are calculated from the original distance matrix through the function distClust(). An intracluster distance is calculated by averaging all pairwise distances among objects in the cluster concerned. Likewise, the distance between two clusters is calculated by averaging all pairwise distances among objects in these clusters.

# Tools for evaluating the clustering outcome

Clustering validation is widely made for hierarchical and iterative methods. Some measures of internal validation for a Tocher's clustering outcome currently implemented on biotools.

# Cophenetic correlation

We have implemented the approach presented by Silva (2013), taking the cluster distances in order to build a cophenetic matrix for clustering performed through the Tocher's method. Their approach consists of taking the cophenetic distance among objects located in the same cluster as the intracluster distance and the cophenetic distance between objects of different clusters as the intercluster distance. Then, the Pearson's correlation between the elements of the original and cophenetic matrix can be taken as a cophenetic correlation. The function to be called is cophenetic(), whose input is an object of class 'tocher' and its output an object of class 'dist'.

# Box's M-test

When clusters are formed, one might be interested in verifying if covariance matrices of the clusters can be considered homogeneous, especially if one intends to perform a linear discriminant analysis using a pooled matrix. In this case, the well known Box's M-test can be applied. The function {\ttfamily boxM()} performs the test using an approximation of the $\chi_{\nu}^2$ distribution, where $\nu = \frac{1}{2} p(p+1)(k-1)$, and $k$ is the number of covariance matrices. Users should be aware that all clusters must have a positive definite covariance matrix.

# Discriminant analysis based on Mahalanobis distance

A simple and efficient object classification rule is based on Mahalanobis distance. It consists of calculating the squared generalized Mahalanobis distances of each multivariate observation to the centre of each cluster. biotools performs this sort of discriminant analysis through {\ttfamily D2.disc()}, as follows: consider $D_{ij, j'}^2$ as the Mahalanobis distance from the $i$-th ($i = 1, 2, ..., n$) observation in the $j$-th ($j = 1, 2, ..., k$) cluster to the probable $j'$-th cluster. Now consider $C_j$ the random variable that represents the cluster at which this observation lies. The predicted class $\hat{C}_{j'}$ for this concerned observation is the one such that $j' \Rightarrow \min_{j' = 1}^{k} ( D_{ij, j'}^2 )$.

D2.disc() outputs the Mahalanobis distances from each observation to the centre of each cluster, the pooled covariance matrix for clusters and the confusion matrix, which contains the number of correct classifications in the diagonal cells and misclassifications off-diagonal. In addition, a column misclass indicates (with an asterisk) where there was disagreement between the original classification (grouping) and the predicted one (pred).

# Instalation

You can install and load the released version of biotools from GitHub with:

```r
devtools::install_github("arsilva87/biotools")

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

cof <- cophenetic(toc)  # cophenetic matrix

cor(cof, garlicdist)  # cophenetic correlation coefficient
[1] 0.9086886
```
