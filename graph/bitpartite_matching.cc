#include <bits/stdc++.h>
#define N 1000

int mat[N];
vector<int> g[N];

bool vis[N];
bool bimatch(int a) {
	vis[a] = 1;
	for(int b : g[a]) {
		if (!mat[b] || !vis[mat[b]] && bimatch(mat[b])) {
			mat[b] = a;
			return 1;
		}
	}
	return 0;
}