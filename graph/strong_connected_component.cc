#include <bits/stdc++.h>

/*------------------------------------------

	strong connected component

input	: directed graph
proc	: make strong connected component

------------------------------------------*/

// Tarjan's algorithm

struct SCC {
	int c, n;
	vector<int> v, scc_n, st;
	vector<vector<int>> scc;

	void init(vector<vector<int>>& g) {
		c = 0;
		n = g.size();
		v.resize(n, 0);
		scc_n.resize(n, 0);
		st.clear();
		scc.clear();
		function<int(int)> dfs = [&](int a) {
			v[a] = ++c;
			int mn = c;
			st.push_back(a);
			for (int b : g[a]) {
				if (!v[b])
					mn = min(mn, dfs(b));
				else if (!scc_n[b])
					mn = min(mn, v[b]);
			}
			if (mn == v[a]) {
				scc.push_back(vector<int>());
				while (st.size()) {
					int tmp = st.back();
					st.pop_back();
					scc_n[tmp] = scc.size();
					scc.back().push_back(tmp);
					if (tmp == a) break;
				}
				sort(scc.back().begin(), scc.back().end());
			}
			return mn;
		};
		for(int i =0; i<n; ++i) if (!v[i]) dfs(i);
	}
};