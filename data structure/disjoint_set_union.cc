#include <bits/stdc++.h>



// USAGE :
// constructive 	- disjoint_set ds(size_t);
// 1. find(u) 		- find u's tag;
// 2. merge(u, v) 	- merge u and v

struct disjoint_set {
	vector<int> a;
	disjoint_set(int n) { a.resize(n, -1); }

	inline int find(int u) {
		return a[u] < 0 ? u : a[u] = find(a[u]);
	}
	inline int merge(int u, int v) {
		u = find(u);
		v = find(v);
		return u == v ? -1 : a[u] += a[v], a[v] = u;
	}
};