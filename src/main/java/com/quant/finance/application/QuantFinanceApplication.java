package com.quant.finance.application;

import com.quant.finance.derivatives.CreditDefaultSwap;
import com.quant.finance.derivatives.Option;
import com.quant.finance.instruments.Bond;
import com.quant.finance.models.MonteCarloPricer;
import com.quant.finance.models.OptionPricingModel;
import com.quant.finance.models.YieldCurve;
import com.quant.finance.models.YieldCurveModels;
import com.quant.finance.risk.CreditRiskModels;
import com.quant.finance.risk.PortfolioRiskManager;
import com.quant.finance.statistics.PortfolioStatistics;
import com.quant.finance.statistics.TimeSeriesAnalysis;
import com.quant.finance.statistics.VolatilityModels;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

/**
 * Example that demonstrates the capabilities of the quantitative finance library
 */
public class QuantFinanceApplication {

    public static void main(String[] args) {
        System.out.println("Quantitative Finance Application");
        System.out.println("============================");
        
        LocalDate today = LocalDate.now();
        
        // Example 1: Bond pricing and analytics
        System.out.println("\n1. Bond Pricing and Analytics");
        System.out.println("----------------------------");
        
        // Create a 5-year bond with 5% coupon rate paid semi-annually
        LocalDate bondIssueDate = today.minusYears(1);
        LocalDate bondMaturityDate = bondIssueDate.plusYears(6);
        Bond corporateBond = new Bond("BOND001", "Corporate Bond 5%", bondIssueDate, 
                                     1000.0, 0.05, 2, bondMaturityDate);
        
        // Value the bond using a 6% yield
        Map<String, Object> bondMarketData = new HashMap<>();
        bondMarketData.put("yieldRate", 0.06);
        
        double bondPrice = corporateBond.presentValue(today, bondMarketData);
        System.out.printf("Bond Price: $%.2f%n", bondPrice);
        
        // Calculate bond risk metrics
        Map<String, Double> bondRisks = corporateBond.calculateRisks(today, bondMarketData);
        System.out.printf("Macaulay Duration: %.2f years%n", bondRisks.get("macaulayDuration"));
        System.out.printf("Modified Duration: %.2f%n", bondRisks.get("modifiedDuration"));
        System.out.printf("Price Change for 1%% Yield Increase: %.2f%%%n", 
                         bondRisks.get("priceChangeFor1PercentYieldChange") * 100);
        System.out.printf("Convexity: %.4f%n", bondRisks.get("convexity"));
        
        // Calculate yield to maturity for a given market price
        double ytm = corporateBond.getYield(today, bondPrice);
        System.out.printf("Yield to Maturity: %.2f%%%n", ytm * 100);
        
        // Example 2: Option pricing with Black-Scholes
        System.out.println("\n2. Option Pricing with Black-Scholes");
        System.out.println("----------------------------------");
        
        // Create a European call option on a stock
        LocalDate optionIssueDate = today.minusMonths(2);
        LocalDate optionExpiryDate = today.plusMonths(4);
        Option callOption = new Option("OPT001", "AAPL Call $150", optionIssueDate, 
                                      Option.OptionType.CALL, Option.ExerciseStyle.EUROPEAN,
                                      150.0, optionExpiryDate, "AAPL");
        
        // Market data for option valuation
        Map<String, Object> optionMarketData = new HashMap<>();
        optionMarketData.put("spotPrice", 145.0); // Current stock price
        optionMarketData.put("riskFreeRate", 0.02); // 2% risk-free rate
        optionMarketData.put("volatility", 0.25); // 25% volatility
        
        // Value the option
        double callPrice = callOption.presentValue(today, optionMarketData);
        System.out.printf("Call Option Price: $%.2f%n", callPrice);
        
        // Calculate option Greeks
        Map<String, Double> greeks = callOption.calculateRisks(today, optionMarketData);
        System.out.printf("Delta: %.4f%n", greeks.get("delta"));
        System.out.printf("Gamma: %.6f%n", greeks.get("gamma"));
        System.out.printf("Theta: %.6f%n", greeks.get("theta"));
        System.out.printf("Vega: %.6f%n", greeks.get("vega"));
        System.out.printf("Rho: %.6f%n", greeks.get("rho"));
        
        // Example 3: Monte Carlo Option Pricing
        System.out.println("\n3. Monte Carlo Option Pricing");
        System.out.println("----------------------------");
        
        // Create a Monte Carlo pricer
        MonteCarloPricer mcPricer = new MonteCarloPricer(10000, 252, true, true, 42);
        
        // Price the same option using Monte Carlo
        Map<String, Double> mcResults = mcPricer.priceOption(callOption, today, 
                                                            (double) optionMarketData.get("spotPrice"),
                                                            (double) optionMarketData.get("volatility"),
                                                            (double) optionMarketData.get("riskFreeRate"),
                                                            0.0); // No dividends
        
        System.out.printf("Monte Carlo Price: $%.2f%n", mcResults.get("price"));
        System.out.printf("Standard Error: $%.4f%n", mcResults.get("standardError"));
        System.out.printf("95%% Confidence Interval: $%.2f to $%.2f%n", 
                         mcResults.get("confidenceLower95"), mcResults.get("confidenceUpper95"));
        
        // Run a sensitivity analysis on volatility
        List<Double> volatilityRange = DoubleStream.iterate(0.1, v -> v + 0.05)
                                                  .limit(9)
                                                  .boxed()
                                                  .collect(Collectors.toList());
        
        Map<Double, Double> volatilitySensitivity = mcPricer.runSensitivityAnalysis(
                callOption, today, 
                (double) optionMarketData.get("spotPrice"),
                (double) optionMarketData.get("volatility"),
                (double) optionMarketData.get("riskFreeRate"),
                0.0, // No dividends 
                "volatility", volatilityRange);
        
        System.out.println("\nVolatility Sensitivity:");
        volatilitySensitivity.forEach((vol, price) -> 
            System.out.printf("Volatility: %.2f, Price: $%.2f%n", vol, price));
        
        // Example 4: Yield Curve Modeling
        System.out.println("\n4. Yield Curve Modeling");
        System.out.println("----------------------");
        
        // Create a sample yield curve
        Map<Double, Double> tenorRates = new HashMap<>();
        tenorRates.put(0.25, 0.0120); // 3 months
        tenorRates.put(0.5, 0.0145);  // 6 months
        tenorRates.put(1.0, 0.0175);  // 1 year
        tenorRates.put(2.0, 0.0205);  // 2 years
        tenorRates.put(3.0, 0.0230);  // 3 years
        tenorRates.put(5.0, 0.0260);  // 5 years
        tenorRates.put(7.0, 0.0280);  // 7 years
        tenorRates.put(10.0, 0.0295); // 10 years
        tenorRates.put(20.0, 0.0310); // 20 years
        tenorRates.put(30.0, 0.0315); // 30 years
        
        YieldCurve yieldCurve = new YieldCurve(today, tenorRates);
        
        // Calculate various rates from the curve
        System.out.println("Zero Rates at Different Tenors:");
        double[] checkTenors = {0.5, 1.0, 2.0, 5.0, 10.0, 30.0};
        for (double tenor : checkTenors) {
            System.out.printf("%.1f year rate: %.2f%%%n", tenor, yieldCurve.getRate(tenor) * 100);
        }
        
        // Calculate forward rates
        System.out.println("\nForward Rates:");
        System.out.printf("Forward 1y1y: %.2f%%%n", yieldCurve.getForwardRate(1.0, 2.0) * 100);
        System.out.printf("Forward 2y3y: %.2f%%%n", yieldCurve.getForwardRate(2.0, 5.0) * 100);
        
        // Price a zero-coupon bond using the yield curve
        double zeroPrice = yieldCurve.priceZeroCouponBond(5.0, 1000.0);
        System.out.printf("\nPrice of 5-year Zero-Coupon Bond: $%.2f%n", zeroPrice);
        
        // Example 5: Portfolio Statistics
        System.out.println("\n5. Portfolio Statistics");
        System.out.println("----------------------");
        
        // Sample stock price history
        List<Double> stockPrices = Arrays.asList(
            100.0, 102.5, 101.8, 103.2, 106.7, 105.3, 107.8, 109.2, 
            108.5, 110.3, 111.8, 114.2, 115.1, 113.7, 116.4, 118.2, 120.5
        );
        
        // Sample benchmark price history
        List<Double> benchmarkPrices = Arrays.asList(
            1000.0, 1010.2, 1005.3, 1012.4, 1025.1, 1022.7, 1035.2, 1042.5, 
            1037.6, 1045.3, 1050.8, 1060.2, 1065.4, 1058.9, 1070.3, 1075.2, 1085.7
        );
        
        // Calculate returns
        double[] stockReturns = PortfolioStatistics.calculateLogReturns(stockPrices);
        double[] benchmarkReturns = PortfolioStatistics.calculateLogReturns(benchmarkPrices);
        
        // Calculate annualized return (assuming daily prices, 252 trading days)
        double annualizedReturn = PortfolioStatistics.calculateAnnualizedReturn(stockPrices, 252);
        System.out.printf("Annualized Return: %.2f%%%n", annualizedReturn * 100);
        
        // Calculate volatility
        double annualizedVol = PortfolioStatistics.calculateAnnualizedVolatility(stockReturns, 252);
        System.out.printf("Annualized Volatility: %.2f%%%n", annualizedVol * 100);
        
        // Calculate Sharpe ratio (assuming 2% risk-free rate)
        double sharpeRatio = PortfolioStatistics.calculateSharpeRatio(stockReturns, 0.02, 252);
        System.out.printf("Sharpe Ratio: %.2f%n", sharpeRatio);
        
        // Calculate maximum drawdown
        double maxDrawdown = PortfolioStatistics.calculateMaxDrawdown(stockPrices);
        System.out.printf("Maximum Drawdown: %.2f%%%n", maxDrawdown * 100);
        
        // Calculate beta
        double beta = PortfolioStatistics.calculateBeta(stockReturns, benchmarkReturns);
        System.out.printf("Beta: %.2f%n", beta);
        
        // Calculate alpha
        double alpha = PortfolioStatistics.calculateAlpha(stockReturns, benchmarkReturns, 0.02, 252);
        System.out.printf("Alpha (annualized): %.2f%%%n", alpha * 100);
        
        // Calculate Value at Risk
        double var95 = PortfolioStatistics.calculateHistoricalVaR(stockReturns, 0.95, 100000.0);
        System.out.printf("95%% 1-day VaR for $100,000 portfolio: $%.2f%n", var95);
        
        // Calculate Conditional Value at Risk (Expected Shortfall)
        double cvar95 = PortfolioStatistics.calculateHistoricalCVaR(stockReturns, 0.95, 100000.0);
        System.out.printf("95%% 1-day CVaR for $100,000 portfolio: $%.2f%n", cvar95);
        
        // Example 6: Time Series Analysis and Volatility Modeling
        System.out.println("\n6. Time Series Analysis and Volatility Modeling");
        System.out.println("-------------------------------------------");
        
        // Synthetic price series for demonstration
        // Create a non-stationary price series with trend and momentum
        double[] syntheticPrices = new double[200];
        syntheticPrices[0] = 100.0;
        for (int i = 1; i < syntheticPrices.length; i++) {
            double trend = 0.05; // upward trend
            double momentum = 0.2; // slight momentum effect
            double vol = 0.015; // volatility
            double randomShock = Math.random() * 2 * vol - vol;
            syntheticPrices[i] = syntheticPrices[i-1] * (1 + trend/252 + momentum * randomShock + randomShock);
        }
        
        // Convert to returns
        double[] syntheticReturns = new double[syntheticPrices.length - 1];
        for (int i = 0; i < syntheticReturns.length; i++) {
            syntheticReturns[i] = Math.log(syntheticPrices[i + 1] / syntheticPrices[i]);
        }
        
        // ADF Test for stationarity
        double adfStat = TimeSeriesAnalysis.adfTest(syntheticPrices, 1);
        System.out.printf("ADF Test Statistic for Prices: %.4f%n", adfStat);
        
        // ADF Test on returns (should be stationary)
        double adfStatReturns = TimeSeriesAnalysis.adfTest(syntheticReturns, 1);
        System.out.printf("ADF Test Statistic for Returns: %.4f%n", adfStatReturns);
        System.out.println("Values below -2.86 suggest stationarity at 5% significance level");
        
        // Hurst Exponent Analysis
        double hurstExponent = TimeSeriesAnalysis.calculateHurstExponent(syntheticReturns);
        System.out.printf("Hurst Exponent of Returns: %.4f%n", hurstExponent);
        System.out.println("Values > 0.5 suggest trending behavior, < 0.5 suggest mean reversion");
        
        // Autocorrelation analysis
        int maxLag = 5;
        double[] acf = TimeSeriesAnalysis.calculateACF(syntheticReturns, maxLag);
        System.out.println("\nAutocorrelation Function (ACF) of Returns:");
        for (int i = 0; i < acf.length; i++) {
            System.out.printf("Lag %d: %.4f%n", i + 1, acf[i]);
        }
        
        // PACF
        double[] pacf = TimeSeriesAnalysis.calculatePACF(syntheticReturns, maxLag);
        System.out.println("\nPartial Autocorrelation Function (PACF) of Returns:");
        for (int i = 0; i < pacf.length; i++) {
            System.out.printf("Lag %d: %.4f%n", i + 1, pacf[i]);
        }
        
        // Volatility modeling with EWMA
        double[] ewmaVol = TimeSeriesAnalysis.calculateEWMAVolatility(syntheticReturns, 0.94);
        System.out.printf("\nEWMA Volatility (last value): %.4f%n", ewmaVol[ewmaVol.length - 1]);
        
        // Volatility modeling with GARCH
        double[] garchVol = TimeSeriesAnalysis.calculateGARCHVolatility(syntheticReturns, 0.1, 0.85, 0.00001);
        System.out.printf("GARCH Volatility (last value): %.4f%n", garchVol[garchVol.length - 1]);
        
        // More advanced volatility modeling using VolatilityModels class
        VolatilityModels historicalModel = VolatilityModels.createModel(syntheticReturns, "historical");
        double forecastHistVol = historicalModel.forecastVolatility();
        System.out.printf("Historical Volatility Forecast: %.4f%n", forecastHistVol);
        
        VolatilityModels garchModel = VolatilityModels.createModel(syntheticReturns, "garch");
        double forecastGarchVol = garchModel.forecastVolatility();
        System.out.printf("GARCH Volatility Forecast: %.4f%n", forecastGarchVol);
        
        // Ljung-Box Test
        double lbStat = TimeSeriesAnalysis.ljungBoxTest(syntheticReturns, 10, 0);
        System.out.printf("\nLjung-Box Test Statistic (Q10): %.4f%n", lbStat);
        System.out.println("High values suggest presence of autocorrelation in the series");
        
        // Calculate and display squared returns' autocorrelation for volatility clustering
        double[] squaredReturns = Arrays.stream(syntheticReturns).map(r -> r * r).toArray();
        double[] acfSquared = TimeSeriesAnalysis.calculateACF(squaredReturns, maxLag);
        System.out.println("\nAutocorrelation of Squared Returns (Volatility Clustering):");
        for (int i = 0; i < acfSquared.length; i++) {
            System.out.printf("Lag %d: %.4f%n", i + 1, acfSquared[i]);
        }
        
        // Example 7: Advanced Option Pricing Models
        System.out.println("\n7. Advanced Option Pricing Models");
        System.out.println("-------------------------------");
        
        // Create options for different exercise styles
        Option europeanCall = new Option("OPT002", "AAPL European Call $150", optionIssueDate, 
                                       Option.OptionType.CALL, Option.ExerciseStyle.EUROPEAN,
                                       150.0, optionExpiryDate, "AAPL");
        
        Option americanCall = new Option("OPT003", "AAPL American Call $150", optionIssueDate, 
                                       Option.OptionType.CALL, Option.ExerciseStyle.AMERICAN,
                                       150.0, optionExpiryDate, "AAPL");
        
        // Using alternative option pricing models
        OptionPricingModel bsModel = OptionPricingModel.createModel("black-scholes");
        OptionPricingModel binomialModel = OptionPricingModel.createModel("binomial");
        OptionPricingModel mcModel = OptionPricingModel.createModel("montecarlo");
        
        // Price European option with different models
        double bsPrice = bsModel.calculatePrice(europeanCall, today, optionMarketData);
        double binomialPrice = binomialModel.calculatePrice(europeanCall, today, optionMarketData);
        double mcModelPrice = mcModel.calculatePrice(europeanCall, today, optionMarketData);
        
        System.out.println("European Call Option Prices:");
        System.out.printf("Black-Scholes: $%.2f%n", bsPrice);
        System.out.printf("Binomial Tree: $%.2f%n", binomialPrice);
        System.out.printf("Monte Carlo: $%.2f%n", mcModelPrice);
        
        // Price American option with Binomial model (Black-Scholes can't price American options)
        double americanPrice = binomialModel.calculatePrice(americanCall, today, optionMarketData);
        System.out.printf("\nAmerican Call Option (Binomial Tree): $%.2f%n", americanPrice);
        
        // Compare Greeks between models
        Map<String, Double> bsGreeks = bsModel.calculateGreeks(europeanCall, today, optionMarketData);
        Map<String, Double> binomialGreeks = binomialModel.calculateGreeks(europeanCall, today, optionMarketData);
        
        System.out.println("\nComparison of Greeks:");
        System.out.println("Black-Scholes vs Binomial Tree");
        System.out.printf("Delta: %.4f vs %.4f%n", bsGreeks.get("delta"), binomialGreeks.get("delta"));
        System.out.printf("Gamma: %.6f vs %.6f%n", bsGreeks.get("gamma"), binomialGreeks.get("gamma"));
        System.out.printf("Theta: %.6f vs %.6f%n", bsGreeks.get("theta"), binomialGreeks.get("theta"));
        
        // Example 8: Alternative Yield Curve Models
        System.out.println("\n8. Alternative Yield Curve Models");
        System.out.println("-------------------------------");
        
        // Create yield curve models with different interpolation methods
        YieldCurveModels linearCurve = YieldCurveModels.createModel(today, tenorRates, "linear");
        YieldCurveModels logLinearCurve = YieldCurveModels.createModel(today, tenorRates, "loglinear");
        YieldCurveModels cubicSplineCurve = YieldCurveModels.createModel(today, tenorRates, "cubicspline");
        YieldCurveModels nelsonSiegelCurve = YieldCurveModels.createModel(today, tenorRates, "nelsonsiegel");
        
        // Compare rates at a non-grid point (e.g., 4-year rate)
        double fourYearTenor = 4.0;
        System.out.printf("4-year Zero Rate Comparison:%n");
        System.out.printf("Linear Interpolation: %.2f%%%n", linearCurve.getRate(fourYearTenor) * 100);
        System.out.printf("Log-Linear Interpolation: %.2f%%%n", logLinearCurve.getRate(fourYearTenor) * 100);
        System.out.printf("Cubic Spline Interpolation: %.2f%%%n", cubicSplineCurve.getRate(fourYearTenor) * 100);
        System.out.printf("Nelson-Siegel Model: %.2f%%%n", nelsonSiegelCurve.getRate(fourYearTenor) * 100);
        
        // Compare forward rates
        System.out.printf("\nForward Rate 3y2y Comparison:%n");
        System.out.printf("Linear Interpolation: %.2f%%%n", linearCurve.getForwardRate(3.0, 5.0) * 100);
        System.out.printf("Log-Linear Interpolation: %.2f%%%n", logLinearCurve.getForwardRate(3.0, 5.0) * 100);
        System.out.printf("Cubic Spline Interpolation: %.2f%%%n", cubicSplineCurve.getForwardRate(3.0, 5.0) * 100);
        System.out.printf("Nelson-Siegel Model: %.2f%%%n", nelsonSiegelCurve.getForwardRate(3.0, 5.0) * 100);
        
        // Example 9: Credit Default Swap Pricing
        System.out.println("\n9. Credit Default Swap Pricing");
        System.out.println("-----------------------------");
        
        // Create a CDS contract
        LocalDate cdsEffectiveDate = today;
        LocalDate cdsMaturityDate = today.plusYears(5);
        
        CreditDefaultSwap cds = new CreditDefaultSwap(
            "CDS001", "XYZ Corp 5Y CDS", today, 
            10000000.0, // $10M notional
            200.0, // 200 bps spread
            4, // quarterly payments
            cdsEffectiveDate,
            cdsMaturityDate,
            "XYZ Corporation",
            0.4 // 40% recovery rate
        );
        
        // Market data for CDS valuation
        Map<String, Object> cdsMarketData = new HashMap<>();
        cdsMarketData.put("riskFreeRate", 0.025);
        cdsMarketData.put("hazardRate", 0.033);
        
        // Calculate present value
        double cdsPV = cds.presentValue(today, cdsMarketData);
        System.out.printf("CDS Present Value: $%.2f%n", cdsPV);
        
        // Calculate risk metrics
        Map<String, Double> cdsRisks = cds.calculateRisks(today, cdsMarketData);
        System.out.printf("CS01 (1bp increase): $%.2f%n", cdsRisks.get("cs01"));
        System.out.printf("IR01 (1bp increase): $%.2f%n", cdsRisks.get("ir01"));
        System.out.printf("Risky Duration: %.2f years%n", cdsRisks.get("riskyDuration"));
        
        // Calculate implied spread
        double impliedSpread = cds.calculateImpliedSpread(today, cdsMarketData);
        System.out.printf("Implied CDS Spread: %.2f bps%n", impliedSpread);
        
        // Example 10: Credit Risk Models
        System.out.println("\n10. Credit Risk Models");
        System.out.println("---------------------");
        
        // Create different credit risk models
        CreditRiskModels mertonModel = CreditRiskModels.createModel("merton");
        CreditRiskModels hazardRateModel = CreditRiskModels.createModel("hazardrate");
        CreditRiskModels creditMetricsModel = CreditRiskModels.createModel("creditmetrics");
        
        // Parameters for Merton model
        Map<String, Object> mertonParams = new HashMap<>();
        mertonParams.put("assetValue", 120.0); // Asset value
        mertonParams.put("assetVolatility", 0.3); // Asset volatility
        mertonParams.put("debtValue", 80.0); // Face value of debt
        mertonParams.put("riskFreeRate", 0.02); // Risk-free rate
        
        // Calculate default probability
        double pd1y = mertonModel.calculateProbabilityOfDefault(1.0, mertonParams);
        System.out.printf("Merton Model 1-year PD: %.2f%%%n", pd1y * 100);
        
        // Parameters for hazard rate model
        Map<String, Object> hazardRateParams = new HashMap<>();
        hazardRateParams.put("hazardRate", 0.04); // Constant hazard rate
        
        // Calculate default probability with hazard rate model
        double pd5y = hazardRateModel.calculateProbabilityOfDefault(5.0, hazardRateParams);
        System.out.printf("Hazard Rate Model 5-year PD: %.2f%%%n", pd5y * 100);
        
        // Parameters for CreditMetrics model
        Map<String, Object> creditMetricsParams = new HashMap<>();
        creditMetricsParams.put("initialRating", "BBB");
        
        // Calculate default probability with CreditMetrics
        double pdCreditMetrics = creditMetricsModel.calculateProbabilityOfDefault(1.0, creditMetricsParams);
        System.out.printf("CreditMetrics 1-year PD for BBB rating: %.2f%%%n", pdCreditMetrics * 100);
        
        // Example 11: Portfolio Risk Management
        System.out.println("\n11. Portfolio Risk Management");
        System.out.println("----------------------------");
        
        // Create a portfolio risk manager
        PortfolioRiskManager riskManager = new PortfolioRiskManager();
        
        // Add positions to the portfolio
        riskManager.addPosition(corporateBond, 100); // 100 bonds
        riskManager.addPosition(cds, 1); // 1 CDS contract
        
        // Calculate instrument risks
        riskManager.calculateInstrumentRisks(today, bondMarketData);
        
        // Calculate portfolio DV01
        double portfolioDV01 = riskManager.calculatePortfolioDV01();
        System.out.printf("Portfolio DV01: $%.2f%n", portfolioDV01);
        
        // Define historical returns for VaR calculation
        Map<String, double[]> historicalReturns = new HashMap<>();
        historicalReturns.put(corporateBond.getId(), stockReturns); // Using stock returns as proxy
        historicalReturns.put(cds.getId(), benchmarkReturns); // Using benchmark returns as proxy
        
        // Calculate historical VaR
        double portfolioVaR = riskManager.calculateHistoricalVaR(historicalReturns, 0.95, 1);
        System.out.printf("95%% 1-day Historical VaR: $%.2f%n", portfolioVaR);
        
        // Calculate parametric VaR
        Map<String, Double> volatilities = new HashMap<>();
        volatilities.put(corporateBond.getId(), 0.05);
        volatilities.put(cds.getId(), 0.10);
        
        Map<String, Map<String, Double>> correlationMatrix = new HashMap<>();
        Map<String, Double> correlations1 = new HashMap<>();
        correlations1.put(corporateBond.getId(), 1.0);
        correlations1.put(cds.getId(), 0.3);
        
        Map<String, Double> correlations2 = new HashMap<>();
        correlations2.put(corporateBond.getId(), 0.3);
        correlations2.put(cds.getId(), 1.0);
        
        correlationMatrix.put(corporateBond.getId(), correlations1);
        correlationMatrix.put(cds.getId(), correlations2);
        
        double parametricVaR = riskManager.calculateParametricVaR(volatilities, correlationMatrix, 0.95, 1);
        System.out.printf("95%% 1-day Parametric VaR: $%.2f%n", parametricVaR);
        
        // Define stress scenarios
        Map<String, Map<String, Object>> stressScenarios = new HashMap<>();
        
        // Interest rate shock scenario
        Map<String, Object> irShockScenario = new HashMap<>();
        irShockScenario.put("riskFreeRate", 0.01); // +100bps shock
        stressScenarios.put("Interest Rate +100bps", irShockScenario);
        
        // Credit spread widening scenario
        Map<String, Object> spreadScenario = new HashMap<>();
        spreadScenario.put("hazardRate", 0.03); // Increased hazard rate
        stressScenarios.put("Credit Spread Widening", spreadScenario);
        
        // Calculate stress test results
        Map<String, Double> stressResults = riskManager.calculateStressTests(today, bondMarketData, stressScenarios);
        
        System.out.println("\nStress Test Results:");
        stressResults.forEach((scenario, impact) -> 
            System.out.printf("%s: $%.2f%n", scenario, impact));
        
        // Risk contribution analysis
        Map<String, Double> riskContributions = riskManager.calculateRiskContribution("dv01");
        
        System.out.println("\nRisk Contributions (DV01):");
        riskContributions.forEach((instrumentId, contribution) -> 
            System.out.printf("%s: %.2f%%%n", instrumentId, contribution * 100));
    }
} 