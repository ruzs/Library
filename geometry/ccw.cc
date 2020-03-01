using point = complex<long long>;
int ccw(point a, point b, point c) {
	long long r = (conj(b - a) * (c - b)).imag();
	return r ? r < 0 ? -1 : 1 : 0;
}
bool cross(point a, point b, point c, point d) {
	return (ccw(a, c, b) * ccw(a, d, b)) <= 0;
}