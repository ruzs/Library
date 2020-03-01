#include <bits/stdc++.h>

template <typename T>
struct segment {
	#define LEFT n<<1,l,l+r>>1
	#define RIGHT n<<1|1,(l+r>>1)+1,r

	vector<T> tree, lazy;

	segment(int n) {
		n = 1 << 33 - __builtin_clz(n);
		tree.resize(n);
		lazy.resize(n);
	}
	segment(const vector<T>& a) {
		segment(a.size());
		modify(1, 0, a.size() - 1, a);
	}
	void modify(int n, int l, int r, const vector<int>& a) {
		if (l == r) {
			tree[n] = a[l];
			return;
		}
		modify(LEFT, a);
		modify(RIGHT, a);
		tree[n] = tree[n << 1] + tree[n << 1 | 1];
	}
	inline void apply(int n, int l, int r) {
		tree[n] += lazy[n] * (r - l + 1);
		if (l != r) {
			lazy[n << 1] += lazy[n];
			lazy[n<<1|1] += lazy[n];
		}
		lazy[n] = 0;
	}
	void update(int n, int l, int r, int s, int e, T x) {
		if (lazy[n]) apply(n, l, r);
		if (e < l || r < s) return;
		if (s <= l && r <= e) {
			lazy[n] += x;
			apply(n, l, r);
			return;
		}
		update(LEFT, s, e, x);
		update(RIGHT, s, e, x);
		tree[n] = tree[n<<1] + tree[n<<1|1];
	}
	int query(int n, int l, int r, int s, int e) {
		if (e < l || r < s) return 0;
		if (lazy[n]) apply(n, l, r);
		if (s <= l && r <= e) return tree[n];
		return query(LEFT, s, e) + query(RIGHT, s, e);
	}
};



template <typename T>
struct segment {
	vector<T> a;
	int n;
	segment(int _n) {
		n = 1 << 32 - __builtin_clz(_n);
		a.resize(n << 1);
	}
	segment(const vector<T>& _a) {
		segment(a.size());
		copy(_a.begin(), _a.end(), a.begin() + n);
		for(int i =n-1; i; --i) a[i] = a[i<<1] + a[i<<1|1];
	}
	void update(int p, T x) {
		a[p += n] = x;
		while(p >>= 1) a[p] = a[p<<1] + a[p<<1|1];
	}
	T query(int s, int e) {
		T r{};
		while(s <= e) {
			if (s % 2 == 1) r += a[s++];
			if (e % 2 == 0) r += a[e--];
			s >>= 1; e >>= 1;
		}
		return r;
	}
};