#include<bits/stdc++.h>


template <typename T>
class graph {
	struct edge {
		int u, v;
		T w;
	};
	vector<edge> e;
	vector<vector<int>> g;

	virtual int add(int a, int b, int c) = 0;

	int add(int a, int b, T c = 1) {
		int id = (int)e.size();
		g[a].push_back(id);
		g[b].push_back(id);
		e.emplace_back(a, b, c);
		return id;
	}
};