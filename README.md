OPTICS Clustering
===

The original OPTICS algorithm is due to [Sander et al][1], and is designed to improve on DBSCAN by taking into account the variable density of the data. OPTICS computes a dendogram based on the reachability of points. The clusters have to be extracted from the reachability, and I use the 'automatic' algorithm, also by [Sander et al][2]

To speed things up I use an [rtree spatial index](http://toblerity.github.io/rtree/) for the nearest neighbour queries.

I'll add some more tests- for the moment there is an ipython notebook which generates test data (or takes it from file), and plots the clusters.


[1]: Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jörg Sander (1999). "OPTICS: Ordering Points To Identify the Clustering Structure". ACM SIGMOD international conference on Management of data. ACM Press. pp. 49–60.

[2]: Sander, Jörg, Xuejie Qin, Zhiyong Lu, Nan Niu, and Alex Kovarsky. "Automatic extraction of clusters from hierarchical clustering representations." In Advances in Knowledge Discovery and Data Mining, pp. 75-87. Springer Berlin Heidelberg, 2003.

[3]: http://toblerity.github.io/rtree/