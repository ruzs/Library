#define X real()
#define Y imag()
using point = complex<long long>;

struct line {
	// ax + by + c = 0
	long long a, b, c, sx = -1e18, ex = 1e18;
	line() {}
	line(const point& u, const point& v) {
		a = u.X - v.X;
		b = v.Y - u.Y;
		c = (v.X - u.X) * u.X - (v.Y - u.Y) * u.Y;
		tie(sx, ex) = minmax(u.X, v.X);
	}
	double diag() {
		return 1.0 * -a / b;
	}
	bool cross(const line& oth) {
		if (a * oth.b == oth.a * b) return 0;
		double cy = double(oth.a * c - a * oth.c) / double(a * oth.b - oth.a * b);
		double cx = double(-b * cy - c) / a;
		return sx < cx && cx < ex;
	}
};