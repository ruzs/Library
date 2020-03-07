#include <bits/stdc++.h>

struct prime {
	vector<int> p;
	bitset<5000005> np;
	prime(int max_n) : p({2}) {
		int i, j, k;
		for (i = 3; i * i <= max_n; i += 2) {
			if (np[i >> 1]) continue;
			p.push_back(i);
			for (j = i * i, k = i + i; j <= max_n; j += k)
				np[j >> 1] = 1;
		}
		while (i <= max_n) {
			if (!np[i >> 1]) p.push_back(i);
			i += 2;
		}
	}
	auto begin() { return p.begin(); }
	auto end() { return p.end(); }

	bool isPrime(int x) {
		return x % 2 == 0 ? x == 2 : !np[x >> 1];
	}
};


//sqrt(x) decision function
auto isPrime = [](auto x) {
	if (x % 2 == 0) return x == 2;
	for(decltype(x) i=3; i*i<=x; i+=2) {
		if (x % i == 0) return 0;
	}
	return 1;
};

//prime factorization
auto factorization = [](auto x) {
	map<decltype(x), int> r;
	for(decltype(x) i =2; i*i<=x; ++i) {
		int c =0;
		while(x % i == 0) c++, x /= i;
		if (c) r[i] = c;
	}
	if (x > 1) r[x] = 1;
	return r;
};

//divisors
auto divisors = [](auto x) {
	set<decltype(x)> r;
	for(decltype(x) i =2; i*i<=x; ++i) {
		if (x % i == 0) {
			r.insert(i);
			r.insert(x / i);
		}
	}
	return r;
};

//Euler phi function, phi(x) = x보다 작은 서로소의 개수
int getEulerPhi(int x) {
	int phi = x;
	for(int i =2; i*i<=x; ++i) {
		if (x % i == 0) phi = phi / i * (i - 1);
		while(x % i == 0) x /= i;
	}
	if (x > 1) phi = phi / x * (x - 1);
	return phi;
}