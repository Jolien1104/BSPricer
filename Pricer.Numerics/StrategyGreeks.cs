using MathNet.Numerics.Distributions;

namespace Pricer.Numerics;

public enum Strategy
{
    Call,
    Put,
    CallSpread,
    Straddle,
    Butterfly
}

public static class StrategyGreeks
{
    private const double Wing = 20.0;

    // ── Normale verdeling hulpfuncties ────────────────────────────────────
    private static double N(double x) => Normal.CDF(0.0, 1.0, x);
    private static double n(double x) => Normal.PDF(0.0, 1.0, x);

    // ── d1 en d2 ──────────────────────────────────────────────────────────
    private static double D1(double S, double K, double r, double sigma, double tau)
        => (Math.Log(S / K) + (r + 0.5 * sigma * sigma) * tau) / (sigma * Math.Sqrt(tau));

    private static double D2(double S, double K, double r, double sigma, double tau)
        => D1(S, K, r, sigma, tau) - sigma * Math.Sqrt(tau);

    // ── Analytische Greeks ────────────────────────────────────────────────

    private static double CallPrice(double S, double K, double r, double sigma, double tau)
    {
        if (tau <= 0) return Math.Max(S - K, 0);
        return S * N(D1(S, K, r, sigma, tau)) - K * Math.Exp(-r * tau) * N(D2(S, K, r, sigma, tau));
    }

    private static double PutPrice(double S, double K, double r, double sigma, double tau)
    {
        if (tau <= 0) return Math.Max(K - S, 0);
        return CallPrice(S, K, r, sigma, tau) - S + K * Math.Exp(-r * tau);
    }

    private static double CallDelta(double S, double K, double r, double sigma, double tau)
    {
        if (tau <= 0) return S > K ? 1.0 : 0.0;
        return N(D1(S, K, r, sigma, tau));
    }

    private static double PutDelta(double S, double K, double r, double sigma, double tau)
        => CallDelta(S, K, r, sigma, tau) - 1.0;

    // Gamma is hetzelfde voor call en put
    private static double Gamma(double S, double K, double r, double sigma, double tau)
    {
        if (tau <= 0) return 0;
        return n(D1(S, K, r, sigma, tau)) / (sigma * S * Math.Sqrt(tau));
    }

    // Theta call: -½·σ·S·n(d1)/√τ - r·K·e^{-rτ}·N(d2)
    private static double CallTheta(double S, double K, double r, double sigma, double tau)
    {
        if (tau <= 0) return 0;
        double d1 = D1(S, K, r, sigma, tau);
        double d2 = D2(S, K, r, sigma, tau);
        return -0.5 * sigma * S * n(d1) / Math.Sqrt(tau)
               - r * K * Math.Exp(-r * tau) * N(d2);
    }

    // Theta put: -½·σ·S·n(d1)/√τ + r·K·e^{-rτ}·N(-d2)
    private static double PutTheta(double S, double K, double r, double sigma, double tau)
    {
        if (tau <= 0) return 0;
        double d1 = D1(S, K, r, sigma, tau);
        double d2 = D2(S, K, r, sigma, tau);
        return -0.5 * sigma * S * n(d1) / Math.Sqrt(tau)
               + r * K * Math.Exp(-r * tau) * N(-d2);
    }

    // Vega is hetzelfde voor call en put: S·√τ·n(d1)
    private static double Vega(double S, double K, double r, double sigma, double tau)
    {
        if (tau <= 0) return 0;
        return S * Math.Sqrt(tau) * n(D1(S, K, r, sigma, tau));
    }

    // ── Strategy aggregator ───────────────────────────────────────────────

    public static (double Price, double Delta, double Gamma, double Theta, double Vega)
        Compute(Strategy strategy, double S, double K, double r, double sigma, double tau)
    {
        return strategy switch
        {
            Strategy.Call => (
                CallPrice(S, K, r, sigma, tau),
                CallDelta(S, K, r, sigma, tau),
                Gamma(S, K, r, sigma, tau),
                CallTheta(S, K, r, sigma, tau),
                Vega(S, K, r, sigma, tau)
            ),

            Strategy.Put => (
                PutPrice(S, K, r, sigma, tau),
                PutDelta(S, K, r, sigma, tau),
                Gamma(S, K, r, sigma, tau),
                PutTheta(S, K, r, sigma, tau),
                Vega(S, K, r, sigma, tau)
            ),

            Strategy.CallSpread => (
                CallPrice(S, K, r, sigma, tau)  - CallPrice(S, K + Wing, r, sigma, tau),
                CallDelta(S, K, r, sigma, tau)  - CallDelta(S, K + Wing, r, sigma, tau),
                Gamma(S, K, r, sigma, tau)      - Gamma(S, K + Wing, r, sigma, tau),
                CallTheta(S, K, r, sigma, tau)  - CallTheta(S, K + Wing, r, sigma, tau),
                Vega(S, K, r, sigma, tau)       - Vega(S, K + Wing, r, sigma, tau)
            ),

            Strategy.Straddle => (
                CallPrice(S, K, r, sigma, tau) + PutPrice(S, K, r, sigma, tau),
                CallDelta(S, K, r, sigma, tau) + PutDelta(S, K, r, sigma, tau),
                2 * Gamma(S, K, r, sigma, tau),
                CallTheta(S, K, r, sigma, tau) + PutTheta(S, K, r, sigma, tau),
                2 * Vega(S, K, r, sigma, tau)
            ),

            Strategy.Butterfly => (
                CallPrice(S, K - Wing, r, sigma, tau) - 2 * CallPrice(S, K, r, sigma, tau) + CallPrice(S, K + Wing, r, sigma, tau),
                CallDelta(S, K - Wing, r, sigma, tau) - 2 * CallDelta(S, K, r, sigma, tau) + CallDelta(S, K + Wing, r, sigma, tau),
                Gamma(S, K - Wing, r, sigma, tau)     - 2 * Gamma(S, K, r, sigma, tau)     + Gamma(S, K + Wing, r, sigma, tau),
                CallTheta(S, K - Wing, r, sigma, tau) - 2 * CallTheta(S, K, r, sigma, tau) + CallTheta(S, K + Wing, r, sigma, tau),
                Vega(S, K - Wing, r, sigma, tau)      - 2 * Vega(S, K, r, sigma, tau)      + Vega(S, K + Wing, r, sigma, tau)
            ),

            _ => throw new ArgumentOutOfRangeException(nameof(strategy))
        };
    }
}
