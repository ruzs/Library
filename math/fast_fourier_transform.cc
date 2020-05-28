#include <bits/stdc++.h>

using cd = complex<double>;
const double PI = acos(-1);

void fft(vector<cd>& a, bool inv) {
    int n = a.size();
    for (int i = 1, j = 0, k; i < n; ++i) {
        for (k = n >> 1; j & k; k >>= 1) j ^= k;
        j ^= k;
        if (i < j) swap(a[i], a[j]);
    }
    for (int k = 2; k <= n; k <<= 1) {
        double ang = 2 * PI / k * (inv ? -1 : 1);
        cd wl(cos(ang), sin(ang));
        for (int i = 0; i < n; i += k) {
            cd w(1);
            for (int j = 0; j < k / 2; ++j, w *= wl) {
                cd u = a[i + j], v = a[i + j + k / 2] * w;
                a[i + j] = u + v;
                a[i + j + k / 2] = u - v;
            }
        }
    }
    if (inv) for (cd & x : a) x /= n;
}
template <typename T>
vector<T> multiply(vector<T>& a, vector<T>& b, T mod = 0) {
    int n = 2;
    while (n <= a.size() + b.size()) n <<= 1;
    vector<cd> A(n), B(n), C(n), D(n);

    for (int i = 0; i < a.size(); ++i) A[i] = {a[i] >> 15, a[i] & 32767};
    for (int i = 0; i < b.size(); ++i) B[i] = {b[i] >> 15, b[i] & 32767};
    fft(A, 0);
    fft(B, 0);
    for (int i = 0; i < n; ++i) {
        int j = i ? n - i : i;
        cd f1 = (A[i] + conj(A[j])) * cd(0.5, 0),
           f2 = (A[i] - conj(A[j])) * cd(0, -0.5),
           f3 = (B[i] + conj(B[j])) * cd(0.5, 0),
           f4 = (B[i] - conj(B[j])) * cd(0, -0.5);
        C[i] = f1 * (f3 + f4 * cd(0, 1));
        D[i] = f2 * (f3 + f4 * cd(0, 1));
    }
    fft(C, 1);
    fft(D, 1);
    vector<T> ans(a.size() + b.size() + 1);
    for (int i = 0; i < ans.size(); ++i) {
        T a = (T)round(C[i].real()),
          b = (T)round(C[i].imag()) + (T)round(D[i].real()),
          c = (T)round(D[i].imag());
        if (mod)
            a %= mod, b %= mod, c %= mod;
        ans[i] = (a << 30) + (b << 15) + c;
        if (mod)
            ans[i] = (ans[i] % mod + mod) % mod;
    }
    return ans;
}

/*------------------------------------------------------------------*/

// xor convolution

void fft(vector<cd>& a, bool inv) {
    int n = a.size();
    for (int i = 1, j = 0, k; i < n; i++) {
        for (k = n >> 1; (j ^= k) < k; k >>= 1);
        if (i < j) swap(a[i], a[j]);
    }
    for (int k = 2; k <= n; k <<= 1) {
        int s = n / k;
        for (int i = 0; i < n; i += k)
            for (int j = 0; j < k / 2; j++) {
                cd v = a[i + j + k / 2];
                a[i + j + k / 2] = a[i + j] - v;
                a[i + j] += v;
            }
    }
    if (inv) for (cd& x : a) x /= n;
}

//
typedef complex<double> base;
void fft(vector<base>& a, bool inv) {
    int n = a.size();
    vector<base> r(n / 2);
    for (int i = 1, j = 0; i < n; i++) {
        int bit = (n >> 1);
        while (j >= bit) {
            j -= bit;
            bit >>= 1;
        }
        j += bit;
        if (i < j)
            swap(a[i], a[j]);
    }
    double ang = 2 * acos(-1) / n * (inv ? -1 : 1);
    for (int i = 0; i < n / 2; i++) { r[i] = base(cos(ang * i), sin(ang * i)); }
    /* In NTT, let prr = primitive root. Then,
    int ang = ipow(prr, (mod - 1) / n);
    if(inv) ang = ipow(ang, mod - 2);
    for(int i=0; i<n/2; i++){
        r[i] = (i ? (1ll * r[i-1] * ang % mod) : 1);
    }
    XOR Convolution : set r[*] = 1.
    OR Convolution : set r[*] = 1, and do following:
    if (!inv) {
        a[j + k] = u + v;
        a[j + k + i/2] = u;
    } else {
        a[j + k] = v;
        a[j + k + i/2] = u - v;
    }
    */
    for (int i = 2; i <= n; i <<= 1) {
        for (int j = 0; j < n; j += i) {
            for (int k = 0; k < i / 2; k++) {
                base u = a[j + k], v = a[j + k + i / 2] * r[n / i * k];
                a[j + k] = u + v;
                a[j + k + i / 2] = u - v;
            }
        }
    }
    if (inv)
        for (auto& x : a) x /= n;
}
template <typename T> vector<T> multiply(vector<T>& v, vector<T>& w) {
    vector<base> fv(v.begin(), v.end()), fw(w.begin(), w.end());
    int n = 2;
    while (n < v.size() + w.size()) n <<= 1;
    fv.resize(n);
    fw.resize(n);
    fft(fv, 0);
    fft(fw, 0);
    for (int i = 0; i < n; i++) fv[i] *= fw[i];
    fft(fv, 1);
    vector<T> ret(n);
    for (int i = 0; i < n; i++) ret[i] = (T)round(fv[i].real());
    return ret;
}
vector<long long> multiply(vector<long long>& v, vector<long long>& w,
                           long long mod) {
    int n = 2;
    while (n < v.size() + w.size()) n <<= 1;
    vector<base> v1(n), v2(n), r1(n), r2(n);
    for (int i = 0; i < v.size(); i++) v1[i] = base(v[i] >> 15, v[i] & 32767);
    for (int i = 0; i < w.size(); i++) v2[i] = base(w[i] >> 15, w[i] & 32767);
    fft(v1, 0);
    fft(v2, 0);
    for (int i = 0; i < n; i++) {
        int j = (i ? (n - i) : i);
        base ans1 = (v1[i] + conj(v1[j])) * base(0.5, 0);
        base ans2 = (v1[i] - conj(v1[j])) * base(0, -0.5);
        base ans3 = (v2[i] + conj(v2[j])) * base(0.5, 0);
        base ans4 = (v2[i] - conj(v2[j])) * base(0, -0.5);
        r1[i] = (ans1 * ans3) + (ans1 * ans4) * base(0, 1);
        r2[i] = (ans2 * ans3) + (ans2 * ans4) * base(0, 1);
    }
    fft(r1, 1);
    fft(r2, 1);
    vector<long long> ret(n);
    for (int i = 0; i < n; i++) {
        long long av = (long long)round(r1[i].real());
        long long bv =
            (long long)round(r1[i].imag()) + (long long)round(r2[i].real());
        long long cv = (long long)round(r2[i].imag());
        av %= mod, bv %= mod, cv %= mod;
        ret[i] = (av << 30) + (bv << 15) + cv;
        ret[i] = (ret[i] % mod + mod) % mod;
    }
    return ret;
}

/*-----------------------------------------------------------------*/
// NTT sait2000
constexpr int powmod(long long a, long long b, int mod) {
    a %= mod;
    a += mod;
    a %= mod;
    int r = 1 % mod;
    for (; b > 0; b /= 2) {
        if (b % 2)
            r = r * a % mod;
        a = a * a % mod;
    }
    return r;
}

constexpr int invmod(long long a, int mod) { return powmod(a, mod - 2, mod); }

template <int prime, int primitive_root> struct FFT {
    static const int P = prime;
    static const int PR = primitive_root;
    static const int PRI = invmod(PR, P);
    static const int MSZ = (P - 1) & (1 - P);
    static const int W = powmod(PR, (P - 1) / MSZ, P);
    static const int WI = powmod(PRI, (P - 1) / MSZ, P);

    static void fft(int sz, int* a, bool inv) {
        int w = powmod(inv ? WI : W, MSZ / sz, P);
        for (int ssz = sz / 2; ssz; ssz /= 2) {
            for (int i = 0; i < sz; i += ssz * 2) {
                long long wt = 1;
                for (int j = 0; j < ssz; ++j) {
                    int& lft = a[i + j];
                    int& rgh = a[i + j + ssz];
                    int lftold = lft;
                    lft = (lft + rgh) % P;
                    rgh = (lftold + (P - rgh)) * wt % P;
                    wt = wt * w % P;
                }
            }
            w = 1LL * w * w % P;
        }
        if (inv) {
            long long invsz = invmod(sz, P);
            for (int i = 0; i < sz; ++i) a[i] = a[i] * invsz % P;
        }
        for (int j = 0, i = 1; i < sz; ++i) {
            int dg = sz / 2;
            for (; j & dg; dg >>= 1) j ^= dg;
            j ^= dg;
            if (i < j) {
                int tmp = a[i];
                a[i] = a[j];
                a[j] = tmp;
            }
        }
    }

    static vector<int> convolution(const vector<int>& a, const vector<int>& b) {
        int sz = a.size() + b.size();
        int fftsz = sz;
        while (fftsz & (fftsz - 1)) fftsz += fftsz & -fftsz;
        vector<int> av(fftsz);
        vector<int> bv(fftsz);
        for (int i = 0, ilim = a.size(); i < ilim; ++i) av[i] = a[i];
        for (int i = 0, ilim = b.size(); i < ilim; ++i) bv[i] = b[i];
        fft(fftsz, av.data(), false);
        fft(fftsz, bv.data(), false);
        for (int i = 0; i < fftsz; ++i) av[i] = 1LL * av[i] * bv[i] % P;
        fft(fftsz, av.data(), true);
        return av;
    }
};

typedef FFT<1012924417, 5> FFT1;
