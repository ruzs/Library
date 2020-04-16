#include <bits/stdc++.h>

long long f(long long x, long long c, long long mod) {
    return ((__int128_t)x * x % mod + c) % mod;
}
long long rho(long long n, long long x0 = 2, long long c = 1) {
    long long x = x0;
    long long y = x0;
    long long g = 1;
	int cnt = 0;
    while (g == 1) {
        x = f(x, c, n);
        y = f(y, c, n);
        y = f(y, c, n);
        g = __gcd(abs(x - y), n);
		if (++cnt > 50000) return -1;
    }
    return g;
}