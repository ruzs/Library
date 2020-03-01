#include <bits/stdc++.h>

/*--------------------------------

	minimum spanning tree

input	: undirected graph's edge list
proc	: making spanning tree with minimum weight
output	: minimum sum of spanning tree weight

----------------------------------*/

// Kruskal's algorithm
struct dsu {
	int a[100005];
	dsu() { memset(a, -1, sizeof a); }
	int find(int u) {
		return a[u] < 0 ? u : a[u] = find(a[u]);
	}
	int merge(int u, int v) {
		u = find(u);
		v = find(v);
		if (u != v) {
			a[u] += a[v];
			a[v] = u;
		}
		return u;
	}
};
int mst(vector<tuple<int, int, int>>& e) {
	dsu ds;
	int r = 0;
	sort(e.begin(), e.end());
	for(auto & [w, a, b] : e) {
		if (ds.find(a) == ds.find(b)) continue;
		r += w;
		ds.merge(a, b);
	}
	return r;
}


// Prim's algorithm

