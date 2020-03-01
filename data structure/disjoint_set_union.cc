#include<bits/stdc++.h>



struct disjoint_set {
	vector<int> a;
	disjoint_set(int n) {
		a.assign(n, -1);
	}
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