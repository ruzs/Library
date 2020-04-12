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
struct rmq {
	vector<T> a;
	int n;
	rmq(const vector<T> &_a) {
		n = 1 << 32 - __builtin_clz(_a.size());
		a.resize(n << 1);
		copy(_a.begin(), _a.end(), a.begin() + n);
		for (int i = n - 1; i; --i)
			a[i] = min(a[i << 1], a[i << 1 | 1], Comp());
	}
	void update(int p, T x) {
		a[p += n] = x;
		while (p >>= 1)
			a[p] = min(a[p << 1], a[p << 1 | 1], Comp());
	}
	T query(int s, int e) {
		T r = min(a[s += n], a[e += n], Comp());
		while (s <= e) {
			if (s & 1)
				r = min(r, a[s++], Comp());
			if (!(e & 1))
				r = min(r, a[e--], Comp());
			s >>= 1;
			e >>= 1;
		}
		return r;
	}
};

// if there a log of query...
// O(n lg n + Q), preprocess: O(n lg n), for each query: O(1)
template <class T, class Cmp = less<T>>
struct rmq {
	vector<array<int, 20>> sp;
	vector<int> lg;

	void assign(vector<int>& a) {
		int n = a.size();
		sp.resize(n);
		lg.resize(n + 1);
		for(int i =0; i<n; ++i) sp[i][0] = a[i];
		for(int i =2; i<=n; ++i) lg[i] = lg[i>>1] + 1;
		for(int j =1; j<20; ++j)
			for(int i =0; i+(1<<j)<=n; ++i)
				sp[i][j] = min(sp[i][j - 1], sp[i + (1<<j-1)][j - 1], Cmp());
	}
	int query(int l, int r) {
		int j = lg[r - l + 1];
		return min(sp[l][j], sp[r - (1 << j) + 1][j], Cmp());
	}
};


template <class T, class Cmp = less<T>> struct rmq {
	vector<array<int, 20>> sp;
	void assign(int* a, int n) {
		sp.resize(n);
		for (int i = 0; i < n; ++i) sp[i][0] = a[i];
		for (int j = 1; j < 20; ++j) for (int i = 0; i + (1 << j) <= n; ++i)
			sp[i][j] = min(sp[i][j - 1], sp[i + (1 << j - 1)][j - 1], Cmp());
	}
	void assign(const vector<int>& a) { assign(a.data(), a.size()); }
	int query(int l, int r) {
		int j = 31 - __builtin_clz(r - l + 1);
		return min(sp[l][j], sp[r - (1 << j) + 1][j], Cmp());
	}
};
