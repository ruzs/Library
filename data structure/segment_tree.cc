#include <bits/stdc++.h>

class segment {
#define LEFT i<<1,l,l+r>>1
#define RIGHT i<<1|1,(l+r>>1)+1,r
    struct data {
		// declare var
        template <typename U>
        void add(const U& v) {
            // define assign rules
        }
    };
	data uni(const data& l, const data& r) {
		// define unite rules
	}
    void propagate(int i, int l, int r) {
        if (lz[i]) {
            st[i].add(lz[i]);
            if (l != r) {
                lz[i << 1] += lz[i];
                lz[i<<1|1] += lz[i]; 
            }
            lz[i] = 0;
        }
    }
    vector<data> st;
    vector<int> sz, lz;
	static int n;
public:
    segment(int N) {
		n = 1 << 32 - __builtin_clz(N);
        st.resize(n << 1);
		sz.resize(n << 1);
        lz.resize(n << 1);
    }
    template <typename U>
    segment(vector<U> & a) {
        segment(a.size());
        for(int i =0; i<a.size(); ++i) {
			st[i + n] = a[i];
			sz[i + n] = 1;
		}
        for(int i =n - 1; i; --i) {
			st[i] = uni(st[i << 1], st[i << 1 | 1]);
			sz[i] = sz[i << 1] + sz[i<<1|1];
		}
    }
    template <typename U>
    void modify(int p, const U & v, int i = 1, int l = 0, int r = n - 1) {
        if (p < l || r < p) return;
        propagate(i, l, r);
        if (l == r) {
			st[i].add(v);
			if (!sz[i]) sz[i] = 1;
			return;
		}
        modify(p, v, LEFT);
        modify(p, v, RIGHT);
        st[i].uni(st[i << 1], st[i<<1|1]);
		sz[i] = sz[i << 1] + sz[i<<1|1];
    }
	data query(int s, int e, int i = 1, int l = 0, int r = n - 1) {
        if (r < s || e < l) return data();
		propagate(i, l, r);
        if (s <= l && r <= e) return st[i];
        return uni(query(s, e, LEFT), query(s, e, RIGHT));
    }
    data all() { return st[1]; }
    data kth(int k, int i = 1) {
        assert(sz[i] <= k);
		if (sz[i] == k) return st[i];
        if (sz[i << 1] >= k) return kth(k, i << 1);
        return kth(k - sz[i << 1], i << 1 | 1);
    }
};