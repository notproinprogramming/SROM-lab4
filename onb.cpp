#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <x86intrin.h>
#define M 173 // m = 173, 2m+1 = 347 (просте), 347 = 3(mod4) -> ОНБ другого типу

enum class OutputFormat { Binary, Hex };

OutputFormat g_outputFormat = OutputFormat::Binary;

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

void PrintLambda(const LambdaMatrix &Lambda) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            std::cout << (Lambda[i][j] ? '1' : '0');
        }
        std::cout << '\n';
    }
}

void SaveLambdaToFile(const LambdaMatrix &Lambda, const std::string &filename) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            out << (Lambda[i][j] ? '1' : '0');
        }
        out << '\n';
    }

    out.close();
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

        /*
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
      */

        void PrintHex(const std::string &name) const {
            std::string hex;
            hex.reserve((M + 3) / 4);

            int pos = 0;

            int first_bits = M % 4;
            if (first_bits == 0)
                first_bits = 4;

            int v = 0;
            for (int b = 0; b < first_bits; ++b) {
                if (bits[pos++]) {
                    v |= 1 << (first_bits - 1 - b);
                }
            }
            hex.push_back(hexChar(v));

            while (pos < M) {
                v = 0;
                for (int b = 0; b < 4; ++b) {
                    if (bits[pos++]) {
                        v |= 1 << (3 - b);
                    }
                }
                hex.push_back(hexChar(v));
            }

            std::cout << name << " = " << hex << std::endl;
        }

        void PrintBin(const std::string &name) const {
            std::cout << name << " =\n";
            for (int i = 0; i < M; ++i) {
                std::cout << (bits[i] ? '1' : '0');
            }
            std::cout << "\n\n";
        }

        void Print(const std::string &name) const {
            if (g_outputFormat == OutputFormat::Binary) {
                PrintBin(name);
            } else {
                PrintHex(name);
            }
        }
        /*
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
        */
        static Galois fromHex(const std::string &hex) {
            Galois res;

            int bitPos = M - 1;

            for (int i = 0; i < (int)hex.size(); ++i) {
                int v = hexValue(hex[i]);

                for (int b = 0; b <= 3; ++b) {
                    if (bitPos >= 0) {
                        res.bits[bitPos] = (v << b) & 1;
                        bitPos--;
                    }
                }
            }

            return res;
        }

        static Galois fromBinaryString(const std::string &s) {
            Galois res;

            int len = std::min((int)s.size(), M);
            for (int i = 0; i < len; ++i) {
                if (s[i] == '0')
                    res.bits[i] = false;
                else if (s[i] == '1')
                    res.bits[i] = true;
                else
                    throw std::runtime_error(
                        "Invalid character in binary string");
            }

            return res;
        }

        static std::vector<bool> exponentBits(const std::string &s) {
            std::vector<bool> bits;

            if (g_outputFormat == OutputFormat::Binary) {
                for (char c : s) {
                    if (c == '0')
                        bits.push_back(false);
                    else if (c == '1')
                        bits.push_back(true);
                    else
                        throw std::runtime_error("Invalid binary exponent");
                }
            } else {
                for (char c : s) {
                    int v = hexValue(c);
                    for (int b = 3; b >= 0; --b) {
                        bits.push_back((v >> b) & 1);
                    }
                }
            }

            return bits;
        }

        Galois operator^(const std::string &exp) const {
            Galois result = Galois::one();
            Galois base = *this;

            std::vector<bool> bits = exponentBits(exp);

            for (bool bit : bits) {
                result = result.square();
                if (bit) {
                    result = result * base;
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

        bool Trace() const {
            bool tr = false;
            for (int i = 0; i < M; ++i) {
                tr = tr ^ bits[i];
            }
            return tr;
        }

        bool operator==(const Galois &other) const {
            for (int i = 0; i < M; ++i) {
                if (bits[i] != other.bits[i])
                    return false;
            }
            return true;
        }

        Galois Inverse() const {
            Galois copy = *this;
            if (copy == Galois::zero()) {
                throw std::runtime_error("Inverse of zero");
            }

            Galois t = copy;

            Galois tmp = copy;

            for (int i = 1; i < M - 1; ++i) {
                tmp = tmp.square();
                t = t * tmp;
            }

            return t.square();
        }
};

int main() {

    std::ofstream fout("result.txt");
    if (!fout) {
        std::cerr << "Cannot open result.txt\n";
        return 1;
    }

    std::streambuf *oldCout = std::cout.rdbuf();
    std::cout.rdbuf(fout.rdbuf());

    Galois one = Galois::one();
    one.Print("1");
    /* Galois A =
     Galois::fromHex("0B5BEB4A7D26A6452F37734773C1EC148B875533996A"); Galois B =
     Galois::fromHex("0F30EFBEDEDDC63D3A7E355B8177B3582F9A1E83C6A9");
     std::string N = "0C7ECEA2D10F386F4A9F1A476B55A06E95CD68F5B197";
 */
    Galois A = Galois::fromBinaryString(
        "1000100001111101110000000100001110101111000101011110000001101101010111"
        "1110011101111111010111110111010000010010010101000111010100010010110101"
        "000110110110101001000010110010110");
    Galois B = Galois::fromBinaryString(
        "0111001011100001001100001100001100100001111111000010010111001101111011"
        "1100011111101011001001101010110001010111110111001100100110100011110101"
        "111111011000111110100000110100000");
    std::string N = "0011101011010101110100011100011001111001100110001000000110"
                    "1000001001101101100101010011101101001100110111010101101100"
                    "110111011001100100000011000011110101000110110010010110001";

    A.Print("A");
    B.Print("B");
    std::cout << "N = " + N
              << "\n "
                 "===================================================\n";

    Galois sum = A + B;
    sum.Print("A + B");

    Galois mul = A * B;
    mul.Print("A * B");

    Galois sq = A.square();
    sq.Print("A^2");

    std::cout << "Trace(A) = " << A.Trace() << std::endl;

    Galois exponention = A ^ N;
    exponention.Print("A^N");

    Galois inv = A.Inverse();
    inv.Print("Inverse(A)");

    std::cout
        << "===================================================\nCheck:\n";
    Galois check = A * inv;
    check.Print("A * Inverse(A) (Should be all 1s)");

    Galois C = mul;
    check = (A + B) * C;
    check.Print("(A+B)*C");
    check = B * C + A * C;
    check.Print("B*C+A*C");

    int Num = 1000;
    uint64_t sqr_cycles =
        measure_cycles([&]() { volatile Galois r = A.square(); }, Num);
    uint64_t mul_cycles =
        measure_cycles([&]() { volatile Galois r = A * B; }, 100);

    std::cout << "\nAverage CPU cycles (ONB):\n";
    std::cout << "Square (Shift): " << sqr_cycles << "\n";
    std::cout << "Multiply:       " << mul_cycles << "\n";

    std::cout.rdbuf(oldCout);
    fout.close();
    return 0;
}
