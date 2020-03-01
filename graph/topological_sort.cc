#include <bits/stdc++.h>
#define N 100000

/*--------------------------

	topological sort
	
input	: directed acyclic graph
output	: topological order list

---------------------------*/

// Tarjan's algorithm

auto tpl_sort(vector<vector<int>>& g) {
	bool v[N]{}, cycle = 0;
	vector<int> tpl;
	function<void(int)> dfs = [&](int a) {
		v[a] = 1;
		for(int b : g[a]) {
			if (v[b]) cycle = 1;
			dfs(b);
		}
		tpl.push_back(a);
	};

	for(int i =0; i<g.size(); ++i) if (!v[i]) dfs(i);
	if (cycle) tpl.clear();
	reverse(tpl.begin(), tpl.end());
	return tpl;
}

// Kahn's algorithm

auto tpl_sort(vector<vector<int>>& g) {
	int ind[N]{};
	vector<int> tpl;
	queue<int> q;

	for(auto & i : g) for(int j : i) ind[j]++;
	for(int i =0; i<g.size(); ++i) if (!ind[i]) q.push(i);
	while(q.size()) {
		int a = q.front(); q.pop();
		tpl.push_back(a);
		for(int b : g[a]) {
			ind[b]--;
			if (!ind[b]) {
				q.push(b);
			}
		}
	}
	return tpl;
}