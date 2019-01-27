---
title: 'biotools: a package for optimization cluster analysis'
authors:
- affiliation: 1
  name: Anderson R da Silva
  orcid: 0000-0003-2518-542X
date: "27 January 2019"
bibliography: paper.bib
tags:
- Tocher algorithm
- clustering techniques
- genetic diversity
affiliations:
- index: 1
  name: Statistics and geoprocessing lab., Instituto Federal Goiano, 75790-000, Urutai, GO, Brazil.
---

# Summary

Tocher's optimization method allows one to establish mutually exclusive clusters, with no need to define the number of clusters. It has been widely used in studies [@Nascimento12, @Ngugi13, @Singh14, @Yatung14] of genetic/phenotypic diversity that are based on cluster analysis. Furthermore, Tocher's method can be used to determine the number of clusters in dendrograms.

Clusters are established according to an objective function that adopts an optimization criterion, which minimizes the average intra-cluster distance and maximizes the average inter-cluster distances [@Silva13]. **biotools** contains the method suggested by K.D. Tocher [@Rao52] for clustering objects, based on the algorithm:

1. (Input) Take a distance matrix, let us say ${\bf d}$, of size $n$.
2. Define a clustering criterion, that is, a limit of distance at which to evaluate the inclusion of objects in a cluster in formation. Usually, this criterion is defined as follows: let ${\bf m}$ be a vector of size $n$ containing the smallest distances involving each object, then take $\theta = \max{({\bf m})}$ to be the clustering criterion.
3. Locate the pair of objects with the smallest distance in ${\bf d}$ to be the starting cluster.
4. Compute the average distance ($d_{k(j)}$) of the starting cluster, say $k$, to a certain object, say $j$, that shows the smallest pairwise distance with each object in the starting cluster.
5. Compare the statistic evaluated at the step 4 with $\theta$. Then, add the new object to the starting cluster if $d_{k(j)} \leq \theta$. Otherwise, go to step 3 not considering the objects already clustered and start a new cluster.

The process continues until the last remaining object is evaluated and either included in the last cluster formed or allocated to a new cluster. The function `tocher()` performs optimization clustering and returns an object of class `tocher`, which contains the clusters formed, a numeric vector indicating the cluster of each object, a matrix of cluster distances and also the input - a class `dist` object.

After obtaining the clusters, it might be useful to know how divergent they are. In this context, cluster distances are calculated from the original distance matrix through the function `distClust()`. An intra-cluster distance is calculated by averaging all pairwise distances among objects in the concerning cluster. Likewise, the distance between two clusters is calculated by averaging all pairwise distances among objects in these clusters.

There are several R packages for cluster analysis, with different outcomes and objective functions. However, none of them offers an implementation of the Tocher's algorithm. The R package **biotools** contains an implementation of Tocher's algorithms, the original and the modified or "sequential" algorithm [@Vasconcelos07],  as well as tools for evaluating the quality of clustering outcome. For this last part, **biotools** supplies some new and standard techniques such as: `cophenetic()` (for class `tocher`) - a specific cophenetic correlation coefficient [@Silva13], `boxM()` - the Box's M-test for evaluating the equality of the cluster covariance matrices and `D2.disc()` - a discriminant analysis based on Mahalanobis distances. In addition, the function `singh()` can be used for determining the importance of variables based on the squared generalized Mahalanobis distance [@Singh81].

# Acknowledgements

To the Conselho Nacional de Desenvolvimento Cientifico e Tecnologico (CNPq, grant number: 307334/2018-0). And to Professor Wojtek Krzanowski for revising the paper.

# References
