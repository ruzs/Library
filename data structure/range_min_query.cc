#include <bits/stdc++.h>

// USAGE

// 1. range_minimum_query
// rmq<int> rm;
// rmq<double> rm;
// else...

// 2. range_maximum_query
// rmq<int, greater<int>> rm;

// 3. Index tracking
// rmq<pair<value, index>> rm;
// idx = rm.query(s, e).second;

// O(n + Q lg n), preprocess: O(n), for each query: O(lg n), for each update: O(lg n)
template <class T, class Comp = less<T>>
struct range_min {
	vector<T> a;
	int n;
	range_min(const vector<T>& _a) {
		n = 1 << 32 - __builtin_clz(_a.size());
		a.resize(n << 1);
		copy(_a.begin(), _a.end(), a.begin() + n);
		for (int i = n - 1; i; --i) a[i] = min(a[i << 1], a[i << 1 | 1], Comp());
	}
	void update(int p, T x) {
		a[p += n] = x;
		while (p >>= 1) a[p] = min(a[p << 1], a[p << 1 | 1], Comp());
	}
	T query(int s, int e) {
		T r = min(a[s += n], a[e += n], Comp());
		while (s <= e) {
			if (s % 2 == 1) r = min(r, a[s++], Comp());
			if (e % 2 == 0) r = min(r, a[e--], Comp());
			s >>= 1;
			e >>= 1;
		}
		return r;
	}
};

// if there a log of query...
// O(n lg n + Q), preprocess: O(n lg n), for each query: O(1)
template <class T, class Comp = less<T>>
struct range_min {
	vector<vector<T>> a;
	range_min(const vector<T>& _a) {
		int d = 32 - __builtin_clz(_a.size());
		a.resize(d, vector<T>(1 << d));
		copy(_a.begin(), _a.end(), a[0].begin());
		for (int i = 0; i < d - 1; ++i)
			for (int j = 0; j < n; ++j)
				a[i + 1][j] = min(a[i][j], a[i][min(n - 1, j + (1 << i))], Comp());
	}
	T query(int s, int e) {
		int d = 31 - __builtin_clz(e - s + 1);
		return min(a[d][s], a[d][e - (1 << d) + 1], Comp());
	}
};