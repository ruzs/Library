#include<bits/stdc++.h>



template <typename T>
struct fenwick {
	vector<T> a;
	int n;
	fenwick(int _n) : n(_n) {
		a.resize(n);
	}
	void add(int p, T x) {
		while(p < n) {
			a[p] += x;
			p += p & -p;
		}
	}
	T sum(int p) {
		T r{};
		while(p) {
			r += a[p];
			p -= p & -p;
		}
		return r;
	}
};