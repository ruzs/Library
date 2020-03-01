#include <bits/stdc++.h>

/*------------------------------------------

	Floyd-Warshall algorithm
	
input	: graph ( 200 x 200 adjacent matrix )
proc 	: making all pair shortest path

------------------------------------------*/

void floyed_warshall(vector<vector<int>>& g) {
	int n = g.size();
	for(int k =0; k<n; ++k) for(int i =0; i<n; ++i) for(int j =0; j<n; ++j) {
		g[i][j] = min(g[i][j], g[i][k] + g[k][j]);
}