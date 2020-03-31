#include <bits/stdc++.h>



auto KMP = [](string & s, string & t) {
	int n = s.size();
	vector<int> pi(n), r;
	for(int i =1, j =0; i<n; ++i) {
		while(j && s[i] != s[j]) j = pi[j - 1];
		if (s[i] == s[j]) pi[i] = ++j;
	}
	for(int i =0, j =0; i<t.size(); ++i) {
		while(j && t[i] != s[j]) j = pi[j - 1];
		if (t[i] == s[j] && ++j == n) r.push_back(i - n + 1);
	}
	return r;
};