namespace Pricer.Numerics;

/// <summary>
/// Generates two correlated Brownian motion paths W1 and W2 using
/// Cholesky decomposition of the 2×2 correlation matrix.
///
/// Correlation structure:
///   W1 = Z1
///   W2 = ρ·Z1 + √(1−ρ²)·Z2
/// where Z1, Z2 are independent standard Brownian motions.
///
/// Diffusive scaling:
///   Each increment is scaled by (Δt)^(α/2) so that:
///     α = 1  → standard Brownian motion  (Var(ΔX) = Δt)
///     α > 1  → freezing                  (Var(ΔX) → 0)
///     α &lt; 1  → explosion                 (Var(ΔX) → ∞)
/// </summary>
public static class CorrelatedBrownianGenerator
{
    /// <summary>
    /// Generate two correlated paths of length <paramref name="steps"/> on [0, 1].
    /// </summary>
    /// <param name="steps">Number of time steps n (Δt = 1/n).</param>
    /// <param name="correlation">Pearson correlation ρ ∈ [−1, 1].</param>
    /// <param name="alpha">Diffusive scaling exponent α (default 1.0 = BM).</param>
    /// <param name="rng">Optional shared Random instance.</param>
    /// <returns>
    /// A tuple (W1, W2) where each array has length <paramref name="steps"/> + 1,
    /// starting at 0 and containing cumulative path values.
    /// </returns>
    public static (double[] W1, double[] W2) Generate(
        int    steps,
        double correlation,
        double alpha = 1.0,
        Random? rng  = null)
    {
        rng ??= new Random();

        // Clamp ρ to a numerically safe range
        double rho  = Math.Clamp(correlation, -1.0 + 1e-9, 1.0 - 1e-9);
        double sqrt1MinusRho2 = Math.Sqrt(1.0 - rho * rho);

        double dt    = 1.0 / steps;
        double scale = Math.Pow(dt, alpha / 2.0);   // (Δt)^(α/2)

        var w1 = new double[steps + 1];
        var w2 = new double[steps + 1];

        w1[0] = 0.0;
        w2[0] = 0.0;

        for (int i = 1; i <= steps; i++)
        {
            // Two independent N(0,1) draws via Box-Muller
            double z1 = SampleNormal(rng);
            double z2 = SampleNormal(rng);

            // Cholesky: decorrelate into W1 and W2
            double dw1 = scale * z1;
            double dw2 = scale * (rho * z1 + sqrt1MinusRho2 * z2);

            w1[i] = w1[i - 1] + dw1;
            w2[i] = w2[i - 1] + dw2;
        }

        return (w1, w2);
    }

    // Box-Muller transform: produces a standard normal sample
    private static double SampleNormal(Random rng)
    {
        double u1 = 1.0 - rng.NextDouble();   // avoid log(0)
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
