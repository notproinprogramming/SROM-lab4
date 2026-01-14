#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <x86intrin.h>
#define M 173 // m = 173, 2m+1 = 347 (просте), 347 = 3(mod4) -> ОНБ другого типу

using LambdaMatrix = std::vector<std::vector<bool>>;

LambdaMatrix computeLambda() {
    constexpr int p = 2 * M + 1;
    LambdaMatrix Lambda(M, std::vector<bool>(M, false));

    std::vector<int> pow2(M);
    int current = 1;
    for (int i = 0; i < M; ++i) {
        pow2[i] = current;
        current = (current * 2) % p;
    }

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            int A = pow2[i];
            int B = pow2[j];

            if ((A + B) % p == 1 || (A - B + p) % p == 1 ||
                (B - A + p) % p == 1 || (p - (A + B) % p) % p == 1) {
                Lambda[i][j] = true;
            }
        }
    }
    return Lambda;
}

LambdaMatrix Lambda = computeLambda();

inline uint64_t rdtsc() { return __rdtsc(); }

template <typename Func>
uint64_t measure_cycles(Func f, int iterations = 1000) {
    uint64_t start = rdtsc();
    for (int i = 0; i < iterations; ++i) {
        f();
    }
    uint64_t end = rdtsc();
    return (end - start) / iterations;
}

static char hexChar(int v) {
    if (v < 10)
        return '0' + v;
    return 'a' + (v - 10);
}

static int hexValue(char c) {
    if ('0' <= c && c <= '9')
        return c - '0';
    if ('a' <= c && c <= 'f')
        return c - 'a' + 10;
    if ('A' <= c && c <= 'F')
        return c - 'A' + 10;
    return 0;
}

class Galois {
    private:
        std::vector<bool> bits;
        static int get_index(int i) {
            constexpr int p = 2 * M + 1;

            int r = i % p;
            if (r < 0)
                r += p;

            if (r == 0)
                r = p - 1;

            if (r > M)
                r = p - r;

            return r - 1;
        }

    public:
        Galois()
            : bits(M, false) {}

        static Galois zero() { return Galois(); }

        static Galois one() {
            Galois result;
            for (int i = 0; i < M; ++i)
                result.bits[i] = true;
            return result;
        }

        Galois(std::initializer_list<bool> values)
            : Galois() {
            int i = 0;
            for (bool val : values) {
                if (i >= M)
                    break;
                bits[i] = val;
                i++;
            }
        }

        Galois operator+(const Galois &other) const {
            Galois result;
            for (int i = 0; i < M; ++i) {
                result.bits[i] = bits[i] ^ other.bits[i];
            }
            return result;
        }

        Galois square() const {
            Galois result;
            for (int i = 0; i < M; ++i) {
                result.bits[(i + 1) % M] = bits[i];
            }
            return result;
        }

        Galois multiply(const Galois &other) const {
            Galois result;

            for (int k = 0; k < M; ++k) {
                bool sum_bit = false;
                for (int i = 0; i < M; ++i) {
                    if (this->bits[(i + k) % M]) {
                        for (int j = 0; j < M; ++j) {
                            if (Lambda[i][j] && other.bits[(j + k) % M]) {
                                sum_bit = !sum_bit;
                            }
                        }
                    }
                }
                result.bits[k] = sum_bit;
            }
            return result;
        }

        Galois operator*(const Galois &other) const {
            return this->multiply(other);
        }
        void Print(const std::string &name) const {
            std::string hex;
            for (int i = 0; i < M; i += 4) {
                int v = 0;
                for (int b = 0; b < 4; ++b) {
                    if (i + b < M && bits[i + b])
                        v |= (1 << b);
                }
                hex.push_back(hexChar(v));
            }
            while (hex.size() > 1 && hex.back() == '0')
                hex.pop_back();
            std::reverse(hex.begin(), hex.end());
            std::cout << name << " = " << hex << std::endl;
        }

        static Galois fromHex(const std::string &hex) {
            Galois res;
            int bitPos = 0;
            for (int i = hex.size() - 1; i >= 0; --i) {
                int v = hexValue(hex[i]);
                for (int b = 0; b < 4; ++b) {
                    if (bitPos < M) {
                        res.bits[bitPos] = (v >> b) & 1;
                        bitPos++;
                    }
                }
            }
            return res;
        }
        // Піднесення до степеня
        Galois operator^(const std::string &hexExponent) const {
            Galois result = Galois::one();
            Galois base = *this;
            for (char c : hexExponent) {
                int v = hexValue(c);
                for (int b = 3; b >= 0; --b) {
                    result = result.square();
                    if ((v >> b) & 1) {
                        result = result * base;
                    }
                }
            }
            return result;
        }

        Galois operator^(int pow) const {
            if (pow == 0)
                return Galois::one();
            Galois res = Galois::one();
            Galois base = *this;
            while (pow > 0) {
                if (pow & 1)
                    res = res * base;
                base = base.square();
                pow >>= 1;
            }
            return res;
        }

        // Слід (Trace) в ОНБ - це просто XOR усіх бітів
        bool Trace() const {
            bool tr = false;
            for (int i = 0; i < M; ++i) {
                tr = tr ^ bits[i];
            }
            return tr;
        }

        // Інверсія за алгоритмом Іто-Цудзіі
        Galois Inverse() const {
            Galois res = *this;
            Galois one_val = Galois::one();
            Galois temp = *this;
            for (int i = 1; i < M; ++i) {
                temp = temp.square() * res;
            }
            return temp.square();
        }
};

int main() {
    Galois one = Galois::one();
    one.Print("1");
    Galois A = Galois::fromHex("0B5BEB4A7D26A6452F37734773C1EC148B875533996A");
    Galois B = Galois::fromHex("0F30EFBEDEDDC63D3A7E355B8177B3582F9A1E83C6A9");
    std::string N = "0C7ECEA2D10F386F4A9F1A476B55A06E95CD68F5B197";

    A.Print("A");
    B.Print("B");
    std::cout << "N = " + N << std::endl;

    Galois sum = A + B;
    sum.Print("A + B");

    Galois mul = A * B;
    mul.Print("A * B");

    Galois sq = A.square();
    sq.Print("A^2");

    std::cout << "Trace(A) = " << A.Trace() << std::endl;

    Galois exponention = A ^ N;
    exponention.Print("A^N");

    /* Galois inv = A.Inverse();
     inv.Print("Inverse(A)");


     Galois check = A * inv;
     check.Print("A * Inverse(A) (Should be all 1s)");

     int Num = 1000;
     uint64_t sqr_cycles =
         measure_cycles([&]() { volatile Galois r = A.square(); }, Num);
     uint64_t mul_cycles =
         measure_cycles([&]() { volatile Galois r = A * B; }, 100);

     std::cout << "\nAverage CPU cycles (ONB):\n";
     std::cout << "Square (Shift): " << sqr_cycles << "\n";
     std::cout << "Multiply:       " << mul_cycles << "\n";
 */
    return 0;
}
