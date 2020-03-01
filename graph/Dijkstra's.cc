#include <bits/stdc++.h>

/*------------------------------------------

	Dijkstra's algorithm
	
input	: graph that doesn't include negative weight ( ~200000 adjacent list )
output	: shortest dist from source

------------------------------------------*/

auto dijkstra(vector<vector<pair<int, int>>>& g, int src) {
	priority_queue<pair<int, int>> pq;
	vector<int> d(g.size(), 1e9);

	pq.push({ 0, src });
	d[src] = 0;

	while(pq.size()) {
		int c = -pq.top().first;
		int a = pq.top().second; pq.pop();
		if (c > d[a]) continue;
		for(auto& [b, w] : g[a]) {
			int nc = c + w;
			if (d[b] > nc) {
				d[b] = nc;
				pq.push({ -nc, b });
			}
		}
	}
	return d;
}