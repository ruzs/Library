template <const int mod>
class modular {
	int x;
public:
	modular() {}
	template <typename U>
	modular(const U& v) : x(v % mod) { if (x < 0) x += mod; }

	template <typename U>
	modular pow(U e) {
		modular a = x, r = 1;
		for (; e > 0; a *= a, e >>= 1) if (e & 1) r *= a;
		return r;
	}
	modular inv() { return pow(mod - 2); }

	template <typename U> auto& operator+=(const U& b) { if ((x += b) >= mod) x -= mod; return *this; }
	template <typename U> auto& operator-=(const U& b) { if ((x -= b) < 0) x += mod; return *this; }
	template <typename U> auto& operator*=(const U& b) { x = b * x - static_cast<long long>(static_cast<long double>(b) * x / mod) * mod; return *this; }
	template <typename U> auto& operator/=(const U& b) { return *this *= modular(b).inv(); }
	auto  operator+ () { return *this; }
	auto  operator- () { return modular(!x ? 0 : mod - x); }
	auto& operator++() { x = x + 1 == mod ? 0 : x + 1; return *this; }
	auto& operator--() { x = !x ? mod - 1 : x - 1; return *this; }
	auto  operator++(int) { return modular(x++); }
	auto  operator--(int) { return modular(x--); }
	auto& operator+=(const modular& rhs) { return *this += rhs.x; }
	auto& operator-=(const modular& rhs) { return *this -= rhs.x; }
	auto& operator*=(const modular& rhs) { return *this *= rhs.x; }
	auto& operator/=(const modular& rhs) { return *this /= rhs.x; }
	auto  operator==(const modular& rhs) { return x == rhs.x; }
	auto  operator!=(const modular& rhs) { return x != rhs.x; }

	template <typename U> friend auto operator+ (const modular& lhs, const U& b) { return modular(lhs) += b; }
	template <typename U> friend auto operator- (const modular& lhs, const U& b) { return modular(lhs) -= b; }
	template <typename U> friend auto operator* (const modular& lhs, const U& b) { return modular(lhs) *= b; }
	template <typename U> friend auto operator/ (const modular& lhs, const U& b) { return modular(lhs) /= b; }
	template <typename U> friend auto operator==(const modular& lhs, const U& b) { return lhs.x == b; }
	template <typename U> friend auto operator!=(const modular& lhs, const U& b) { return lhs.x != b; }
	friend string to_string(modular& m) { return to_string(m.x); }
	friend auto& operator>>(istream& in, modular& m) { cin >> m.x; return in; }
	friend auto& operator<<(ostream& out, modular& m) { cout << m.x; return out; }
};

using mint = modular<998244353>;
