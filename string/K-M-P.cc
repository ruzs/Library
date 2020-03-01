#include <bits/stdc++.h>

auto getPi(string& p) {
	vector<int> r(p.size());
	for (int i = 1, j = 0; i < p.size(); ++i) {
		while (j && p[i] != p[j]) j = r[j - 1];
		if (p[i] == p[j]) r[i] = ++j;
	}
	return r;
}
auto KMP(string& h, string& t) {
	vector<int> r, pi = getPi(t);
	for (int i = 0, j = 0; i < h.size(); ++i) {
		while (j && h[i] != t[j]) j = pi[j - 1];
		if (h[i] == t[j]) {
			if (j == t.size() - 1) {
				r.push_back(i - t.size() + 1);
				j = pi[j];
			}
			else j++;
		}
	}
	return r;
}