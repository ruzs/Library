#include <bits/stdc++.h>

auto suffixArray(const string& s) {
	int n = s.size(), d = 1;
	vector<int> sa(n), g(n);
	auto cmp = [&](int i, int j) {
		if (g[i] != g[j]) return g[i] < g[j];
		i += d;
		j += d;
		return i < n && j < n ? g[i] < g[j] : i > j;
	};
	for (int i =0; i<n; ++i) {
		sa[i] = i;
		g[i] = s[i];
	}
	while (1) {
		sort(sa.begin(), sa.end(), cmp);
		vector<int> ng(n, 0);
		for (int i =1; i<n; ++i)
			ng[sa[i]] = ng[sa[i - 1]] + cmp(sa[i - 1], sa[i]);
		if (ng[sa[n - 1]] == n - 1) break;
		g = ng;
		d <<= 1;
	}
	return sa;
}