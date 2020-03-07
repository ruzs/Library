#include <bits/stdc++.h>


template <typename T>
struct fenwick {
	vector<T> a;
	fenwick(int n) { a.resize(n + 1); }

	void add(int p, T x) {
		do a[p] += x;
		while ((p += p & -p) < a.size());
	}
	T sum(int p) {
		T r{};
		do r += a[p];
		while (p -= p & -p);
		return r;
	}
};