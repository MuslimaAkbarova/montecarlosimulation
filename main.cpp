#include <iostream>
#include <random>
#include <cmath>
#include <chrono>

using namespace std;

double calculateCallPriceMC(double S, double K, double r, double sigma, double T, int numSimulations)
{
    double sumCallPayoffs = 0.0;

    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> distribution(0.0, 1.0);

    double dt = T / static_cast<double>(numSimulations);

    for (int i = 0; i < numSimulations; ++i)
    {
        double epsilon = distribution(gen);
        double S_T = S * exp((r - 0.5 * sigma * sigma) * T + sigma * sqrt(T) * epsilon);

        double callPayoff = std::max(S_T - K, 0.0);
        sumCallPayoffs += callPayoff;
    }

    double averageCallPayoff = sumCallPayoffs / static_cast<double>(numSimulations);
    double callPrice = exp(-r * T) * averageCallPayoff;

    return callPrice;
}

double calculatePutPriceMC(double S, double K, double r, double sigma, double T, int numSimulations)
{
    double sumPutPayoffs = 0.0;

    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> distribution(0.0, 1.0);

    double dt = T / static_cast<double>(numSimulations);

    for (int i = 0; i < numSimulations; ++i)
    {
        double epsilon = distribution(gen);
        double S_T = S * exp((r - 0.5 * sigma * sigma) * T + sigma * sqrt(T) * epsilon);

        double putPayoff = std::max(K - S_T, 0.0);
        sumPutPayoffs += putPayoff;
    }

    double averagePutPayoff = sumPutPayoffs / static_cast<double>(numSimulations);
    double putPrice = exp(-r * T) * averagePutPayoff;

    return putPrice;
}

int main()
{
    // Input parameters
    double S;       // Current price of the underlying asset
    double K;       // Strike price of the option
    double r;      // Risk-free interest rate
    double sigma;   // Volatility of the underlying asset
    double T;         // Time to expiration (in years)
    int numSimulations = 1000000;  // Number of Monte Carlo simulations

    cout << "Current price of the underlying asset:" <<endl;
    cin >> S;
    cout << "Strike price of the option:" <<endl;
    cin >> K;
    cout << "Risk-free interest rate:" <<endl;
    cin >> r;
    cout << "Volatility of the underlying asset:" <<endl;
    cin >> sigma;
    cout << "Time to expiration (in years):" <<endl;
    cin >> T;

    auto startTime = chrono::high_resolution_clock::now();
    // Calculate call price using Monte Carlo simulation
    double callPrice = calculateCallPriceMC(S, K, r, sigma, T, numSimulations);

    // Calculate put price using Monte Carlo simulation
    double putPrice = calculatePutPriceMC(S, K, r, sigma, T, numSimulations);

    auto endTime = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(endTime - startTime);
    // Output the results
    cout << "Call Price: " << callPrice << endl;
    cout << "Put Price: " << putPrice << endl;
    cout << "Execution Time: " << duration.count() << " milliseconds" << endl;
    return 0;
}
