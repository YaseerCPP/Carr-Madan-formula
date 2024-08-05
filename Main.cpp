#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

const double PI = 3.14159265358979323846;

// Black-Scholes characteristic function
std::complex<double> characteristicFunction(double v, double S0, double r, double sigma, double T) {
    std::complex<double> i(0, 1);
    std::complex<double> phi = exp(i * v * (log(S0) + (r - 0.5 * sigma * sigma) * T) - 0.5 * sigma * sigma * v * v * T);
    return phi;
}

// Carr-Madan formula for European call option pricing
double carrMadanCallPrice(double K, double S0, double r, double sigma, double T, double alpha = 1.5, int N = 1000) {
    double deltaV = 0.1;
    double callPrice = 0.0;

    for (int k = 1; k <= N; ++k) {
        double v = k * deltaV;
        std::complex<double> i(0, 1);
        std::complex<double> numerator = exp(-i * v * log(K)) * characteristicFunction(v - i * (alpha + 1), S0, r, sigma, T);
        std::complex<double> denominator = alpha * alpha + alpha - v * v + i * (2.0 * alpha + 1.0) * v;
        std::complex<double> integrand = exp(-r * T) / PI * numerator / denominator;
        callPrice += deltaV * real(integrand);
    }

    return callPrice;
}

int main() {
    double S0 = 100.0;     // Initial price of the underlying asset
    double K = 100.0;      // Strike price
    double r = 0.05;       // Risk-free interest rate
    double sigma = 0.2;    // Volatility
    double T = 1.0;        // Time to maturity

    double callPrice = carrMadanCallPrice(K, S0, r, sigma, T);

    std::cout << "European Call Option Price: " << callPrice << std::endl;

    return 0;
}
