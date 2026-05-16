namespace Pricer.Numerics;

public static class BinomialTree
{
    public static double EuropeanCall(double s, double k, double r, double sigma, double t, int n)
        => Price(s, k, r, sigma, t, n, OptionType.Call, isAmerican: false);

    public static double EuropeanPut(double s, double k, double r, double sigma, double t, int n)
        => Price(s, k, r, sigma, t, n, OptionType.Put, isAmerican: false);

    public static double AmericanCall(double s, double k, double r, double sigma, double t, int n)
        => Price(s, k, r, sigma, t, n, OptionType.Call, isAmerican: true);

    public static double AmericanPut(double s, double k, double r, double sigma, double t, int n)
        => Price(s, k, r, sigma, t, n, OptionType.Put, isAmerican: true);

    public static IReadOnlyList<ConvergencePoint> ConvergenceProfile(
        double s, double k, double r, double sigma, double t,
        IEnumerable<int> stepCounts, bool isAmerican = false)
    {
        double bsRef = new BlackScholes(OptionType.Call, r, t, sigma, k, s).Price();
        return stepCounts
            .Select(n => new ConvergencePoint(
                n,
                Price(s, k, r, sigma, t, n, OptionType.Call, isAmerican),
                bsRef))
            .ToList();
    }

    private static double Price(
        double s, double k, double r, double sigma, double t,
        int n, OptionType type, bool isAmerican)
    {
        if (n <= 0) throw new ArgumentOutOfRangeException(nameof(n));

        double dt   = t / n;
        double u    = Math.Exp(sigma * Math.Sqrt(dt));
        double d    = 1.0 / u;
        double disc = Math.Exp(-r * dt);
        double p    = (Math.Exp(r * dt) - d) / (u - d);

        double[] v = new double[n + 1];
        for (int j = 0; j <= n; j++)
            v[j] = Payoff(s * Math.Pow(u, j) * Math.Pow(d, n - j), k, type);

        for (int i = n - 1; i >= 0; i--)
        {
            for (int j = 0; j <= i; j++)
            {
                double continuation = disc * (p * v[j + 1] + (1 - p) * v[j]);
                v[j] = isAmerican
                    ? Math.Max(continuation, Payoff(s * Math.Pow(u, j) * Math.Pow(d, i - j), k, type))
                    : continuation;
            }
        }

        return v[0];
    }

    private static double Payoff(double spot, double k, OptionType type) =>
        type == OptionType.Call ? Math.Max(spot - k, 0.0) : Math.Max(k - spot, 0.0);
}

public record ConvergencePoint(int Steps, double BinomialPrice, double BlackScholesPrice)
{
    public double AbsoluteError => Math.Abs(BinomialPrice - BlackScholesPrice);
    public double RelativeError => BlackScholesPrice > 1e-10 ? AbsoluteError / BlackScholesPrice : 0;
}