#include <bits/stdc++.h>

/*------------------------------------

	Hopcroft-Karp algorithm
	
input	: bipartite graph
output	: maximum matching

------------------------------------*/


int bipartite_matching(vector<vector<int>>& g) {
	int n = g.size();
	vector<int> mat(n, 0), use(n, 0), lv(n);
	auto bfs = [&]() {
		queue<int> q;
		for(int i =1; i<=n; ++i) {
			if (!use[i]) lv[i] = 0, q.push(i);
			else lv[i] = -1;
		}
		while(q.size()) {
			int a = q.front(); q.pop();
			for(int b : g[a])
				if (mat[b] && lv[mat[b]] == -1)
					lv[mat[b]] = lv[a] + 1,
					q.push(mat[b]);
		}
	};
	function<bool(int)> bimatch = [&](int a) {
		for(int b : g[a]) {
			if (!mat[b] || lv[a] < lv[mat[b]] && bimatch(mat[b])) {
				use[a] = 1;
				mat[b] = a;
				return 1;
			}
		}
		return 0;
	};
	int ans = 0, fl;
	do {
		bfs();
		fl = 0;
		for(int i =1; i<=n; ++i)
			if (!use[i] && bimatch(i)) fl++;
		ans += fl;
	} while(fl);
	return ans;
}