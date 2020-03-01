#include <bits/stdc++.h>

/*------------------------------------------
	Bellman-Ford algorithm
	
input	: graph ( ~10000 adjacent list )
output	: shortest dist from source
+ find negative cycle

------------------------------------------*/

auto bellman_ford(vector<vector<pair<int, int>>>& g, int src) {
	int n = g.size();
	vector<int> d(n, 1e9), minusCycle = {-1};
	d[src] = 0;

	for(int i =1; i<=n; ++i) {
		bool update = 0;
		for(int a =1; a<=n; ++a) {
			if (d[a] == 1e9) continue;
			for(auto& [b, c] : g[a]) {
				if (d[b] > d[a] + c) {
					d[b] = d[a] + c;
					update = 1;
					if (i == n) return minusCycle;
				}
			}
		}
		if (!update) break;
	}
	return d;
}