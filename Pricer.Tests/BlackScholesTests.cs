namespace Pricer.Tests;

using Pricer.Numerics;
using Xunit;



public class BlackScholesTests
{
    [Fact]
    public void CallPrice_ShouldMatchKnownValue()
    {
        var bs = new BlackScholes(OptionType.Call, 0.05, 1.0, 0.2, 100, 100);
        Assert.Equal(10.4506, bs.Price(), precision: 3);
    }

    [Fact]
    public void PutCallParity_ShouldHold()
    {
        var call = new BlackScholes(OptionType.Call, 0.05, 1.0, 0.2, 100, 100);
        var put  = new BlackScholes(OptionType.Put,  0.05, 1.0, 0.2, 100, 100);

        double lhs = call.Price() - put.Price();
        double rhs = 100 * Math.Exp(0) - 100 * Math.Exp(-0.05 * 1.0);

        Assert.Equal(lhs, rhs, precision: 6);
    }

    [Fact]
    public void ImpliedVol_ShouldRecoverOriginalSigma()
    {
        double sigma = 0.2;
        var bs = new BlackScholes(OptionType.Call, 0.05, 1.0, sigma, 100, 100);
        double price = bs.Price();

        double iv = ImpliedVolatilityCalculator.Compute(OptionType.Call, price, 0.05, 1.0, 100, 100);

        Assert.Equal(sigma, iv, precision: 6);
    }
}

