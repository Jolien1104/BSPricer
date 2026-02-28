using MathNet.Numerics.Distributions;


namespace Pricer.Numerics;

public enum OptionType
{
    Call,
    Put
}

public class BlackScholes
{
    private OptionType option;
    private double r; // Risk-free interest rate
    private double T; // Time to maturity
    private double sigma; // Volatility
    private double K; // Strike price
    private double S; // Underlying asset price
    private double q; // Dividend yield

    public BlackScholes(
        OptionType optionType,
        double riskFreeRate,
        double timeToMaturity,
        double volatility,
        double strike,
        double underlyingPrice,
        double dividendYield = 0.0)
    {
        option = optionType;
        r = riskFreeRate;
        T = timeToMaturity;
        sigma = volatility;
        K = strike;
        S = underlyingPrice;
        q = dividendYield;
    }

    // Cumulative normal distribution
    private static double N(double x)
        => Normal.CDF(0.0, 1.0, x);

    // Normal probability density function
    private static double n(double x)
        => Normal.PDF(0.0, 1.0, x);


    // Black–Scholes price
    public double Price()
    {
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * sigma * sigma) * T) 
                    / (sigma * Math.Sqrt(T));
        double d2 = d1 - sigma * Math.Sqrt(T);

        if (option == OptionType.Call)
            return S * Math.Exp(-q * T) * N(d1) - K * Math.Exp(-r * T) * N(d2);
        else
            return K * Math.Exp(-r * T) * N(-d2) - S * Math.Exp(-q * T) * N(-d1);
    }

    public double Vega()
        => S * Math.Exp(-q * T) * Math.Sqrt(T) * n(
            (Math.Log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * Math.Sqrt(T))
        );

}

// Test price compared with party at the mooonlight option calculator
// Test put call parity for consistency


public static class ImpliedVolatilityCalculator
{
    public static double Compute(
        OptionType optionType,
        double marketPrice,
        double interestRate,
        double timeToMaturity,
        double strike,
        double underlyingPrice,
        double initialGuess = 0.2,
        double tolerance = 1e-8,
        int maxIterations = 100)
    {
        double sigma = initialGuess;

        for (int i = 0; i < maxIterations; i++)
        {
            var bs = new BlackScholes(optionType, interestRate, timeToMaturity, sigma, strike, underlyingPrice);
            double price = bs.Price();
            double vega = bs.Vega();

            double diff = price - marketPrice;
            if (Math.Abs(diff) < tolerance) return sigma;

            sigma -= diff / vega; // Newton-Raphson stap
        }

        return sigma; // beste schatting na maxIterations
    }

    public static double BachelierImpliedVolATM(
        double optionPrice,
        double underlyingPrice,
        double timeToMaturity,
        double interestRate)
    {
        return (optionPrice / underlyingPrice) * Math.Sqrt(2 * Math.PI / timeToMaturity);
    }
}