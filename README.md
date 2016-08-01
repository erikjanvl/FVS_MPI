# FVS_MPI

Short description:
The algorithm uses a branching algorithm that branches on the vertex of highest possible degree. It uses the Bafna et al algorithm to find a good upper and lower bound to guide the branching. The Bafna lower bounds are non-standard. In each branching step, some of the reduction rules of the Bodlaender-Van Dijk cubic kernel are used, in particular deletion of nodes of degree 0, 1, contraction of nodes of degree 2, and the strongly forced vertex rule (which also relies on Bafna). We also added a rule that gets rid of nodes of degree 2 with one neighbor.

Implementation:
We use our own graph class, which is a combination of an adjacency matrix and an adjacency list. A particular strength is that deletion of nodes can be reversed in time O(degree), which is of great help during our branching algorithm.

Getting it to work:
Our implementation requires lemon (http://lemon.cs.elte.hu/). Assuming lemon is installed and part of the default search path, the whole thing compiles with:
c++ -o fvs.o -std=c++11 -stdlib=libc++ main.cpp reductions.cpp -DLEMON_ONLY_TEMPLATES