# Quantitative Finance Library

A Java-based library for financial mathematics, derivative pricing, risk management, and portfolio analysis.

## Overview

This library demonstrates a range of quantitative finance applications including:

- Fixed income analytics (bond pricing, duration, convexity)
- Derivative pricing (Black-Scholes, Monte Carlo methods)
- Yield curve construction and analysis
- Portfolio statistics and risk metrics
- Option greeks calculation

## Features

### Financial Instruments

- Base instrument interface with common valuation methods
- Fixed-income bond implementation with cashflow analysis
- Options (European and American) with various pricing models

### Financial Models

- Yield curve construction and interpolation
- Nelson-Siegel curve fitting
- Black-Scholes option pricing
- Monte Carlo simulation for derivative pricing

### Statistical Analysis

- Return calculation (simple, log)
- Volatility, Sharpe ratio, maximum drawdown
- Alpha and Beta calculation
- Value at Risk (VaR) and Conditional VaR

### Utilities

- Date utilities with various day count conventions
- Market data management

## Getting Started

### Prerequisites

- Java 17 or higher
- Maven 3.x

### Building

```bash
mvn clean install
```

### Running the Examples

```bash
mvn exec:java -Dexec.mainClass="com.quant.finance.application.QuantFinanceApplication"
```

## Example Usage

The library provides a comprehensive API for financial calculations. Example usage:

```java
// Create a bond
Bond corporateBond = new Bond("BOND001", "Corporate Bond 5%", issueDate, 
                            1000.0, 0.05, 2, maturityDate);

// Value the bond
double bondPrice = corporateBond.presentValue(valuationDate, bondMarketData);

// Calculate bond risk metrics
Map<String, Double> bondRisks = corporateBond.calculateRisks(valuationDate, bondMarketData);
double duration = bondRisks.get("macaulayDuration");
double convexity = bondRisks.get("convexity");

// Create and price an option
Option callOption = new Option("OPT001", "AAPL Call $150", issueDate, 
                             Option.OptionType.CALL, Option.ExerciseStyle.EUROPEAN,
                             150.0, expiryDate, "AAPL");

double optionPrice = callOption.presentValue(valuationDate, optionMarketData);
```

## Advanced Features

### Monte Carlo Simulation

```java
// Create a Monte Carlo pricer with 10000 simulations
MonteCarloPricer mcPricer = new MonteCarloPricer(10000, 252, true, true, 42);

// Price an option using Monte Carlo
Map<String, Double> results = mcPricer.priceOption(option, valuationDate, spotPrice,
                                                 volatility, riskFreeRate, 0.0);

// Run sensitivity analysis
Map<Double, Double> sensitivity = mcPricer.runSensitivityAnalysis(
        option, valuationDate, spotPrice, volatility, riskFreeRate,
        0.0, "volatility", volatilityRange);
```

### Yield Curve Analysis

```java
// Build a yield curve
YieldCurve yieldCurve = new YieldCurve(today, tenorRates);

// Get rates at different tenors
double rate5y = yieldCurve.getRate(5.0);

// Calculate forward rates
double forward2y3y = yieldCurve.getForwardRate(2.0, 5.0);

// Price a zero-coupon bond
double zeroPrice = yieldCurve.priceZeroCouponBond(5.0, 1000.0);
```

### Portfolio Analytics

```java
// Calculate returns
double[] returns = PortfolioStatistics.calculateLogReturns(prices);

// Calculate annualized return
double annReturn = PortfolioStatistics.calculateAnnualizedReturn(prices, 252);

// Calculate risk metrics
double sharpeRatio = PortfolioStatistics.calculateSharpeRatio(returns, 0.02, 252);
double var95 = PortfolioStatistics.calculateHistoricalVaR(returns, 0.95, 100000.0);
```

## Package Structure

- `com.quant.finance.instruments` - Financial instrument base classes and implementations
- `com.quant.finance.derivatives` - Option pricing and other derivative classes
- `com.quant.finance.models` - Financial models like yield curves and Monte Carlo
- `com.quant.finance.statistics` - Statistical and portfolio analysis tools
- `com.quant.finance.utils` - Utility classes for date calculations and other helpers
- `com.quant.finance.application` - Application classes demonstrating library usage

## License

This project is licensed under the MIT License - see the LICENSE file for details. 



### 
Key components I've observed:
Financial Instruments: You've created a robust inheritance hierarchy with the abstract Instrument base class that defines common interfaces like presentValue() and calculateRisks(), which are implemented by specific instruments like Bond, Option, InterestRateSwap, and CreditDefaultSwap.
Derivative Pricing: Your implementations include advanced models like Black-Scholes for options and hazard rate modeling for credit default swaps. The CreditDefaultSwap class shows sophisticated handling of credit risk with methods to calculate premium and protection legs.
Statistical Tools: The PortfolioStatistics class provides comprehensive risk metrics including Sharpe ratio, VaR, CVaR, maximum drawdown, and alpha/beta calculations.
Financial Models: Your YieldCurve implementation handles construction, interpolation, and forward rate calculations, while the MonteCarloPricer provides simulation-based pricing with variance reduction techniques.
Risk Management: The PortfolioRiskManager handles stress testing and risk contribution analysis.
The code quality is high, with well-documented methods, appropriate use of design patterns, and clear separation of concerns. The fix you mentioned for the CreditDefaultSwap class is evident in the implementation of the overloaded calculatePremiumLegPV() method that takes a custom spread parameter instead of modifying the final spreadBps field.

###
Financial Mathematics Topics
Bond valuation and pricing
Present value calculation techniques
Duration and convexity analysis
Option pricing (Black-Scholes model)
Greeks calculation (delta, gamma, theta, vega, rho)
Monte Carlo simulation for derivatives
Yield curve construction and interpolation
Forward rate calculation
Credit default swap pricing
Hazard rate modeling
Interest rate swap valuation
Zero-coupon bond pricing
Discount factor derivation

Risk Management Topics
Value at Risk (VaR) calculation methods
Conditional Value at Risk (CVaR)/Expected Shortfall
Stress testing frameworks
Risk contribution analysis
DV01/Interest rate sensitivity
Credit spread risk (CS01)
Recovery rate sensitivity
Parametric VaR using variance-covariance
Historical VaR simulation
Risk budgeting and allocation
Portfolio stress scenarios
Correlation matrix construction
Risk-adjusted performance metrics


Statistics Topics
Return calculations (log, simple)
Volatility estimation techniques
EWMA (Exponentially Weighted Moving Average)
GARCH volatility modeling
Autocorrelation and PACF analysis
Stationarity testing (ADF test)
Hurst exponent for trend detection
Ljung-Box test for residual analysis
Time series analysis
Maximum drawdown calculation
Sharpe ratio
Alpha and Beta calculation
Regression analysis
Descriptive statistics
Correlation metrics

###
Excellent! Now all the implementations are complete. Let me summarize what we've implemented to fulfill all the requirements:
Financial Mathematics Topics:
We created FinancialMathUtils.java with fundamental calculations for discount factors, present value, IRR, NPV, bond mathematics, and interpolation methods.
Implemented FixedIncomeAnalytics.java for bond valuation, duration/convexity, key rate durations, and yield curve construction.
Created OptionPricingModel.java with Black-Scholes, Binomial Tree, and Monte Carlo pricing models, including Greeks calculations.
Developed YieldCurveModels.java with multiple interpolation methods including Linear, Log-Linear, Cubic Spline, and Nelson-Siegel.
Risk Management Topics:
Implemented CreditRiskModels.java with Merton, Hazard Rate, and CreditMetrics models for credit risk assessment.
Enhanced PortfolioRiskManager.java for VaR, CVaR, stress testing, and risk contribution analysis.
Added methods for calculating DV01, CS01, and recovery rate sensitivity.
Implemented correlation matrix construction and risk budgeting functionality.
Statistics Topics:
Created VolatilityModels.java with Historical, EWMA, and GARCH models for volatility estimation.
Enhanced TimeSeriesAnalysis.java with autocorrelation, stationarity testing, and Hurst exponent calculation.
Implemented GARCH parameter estimation via maximum likelihood.
Added statistical functions for return calculations, descriptive statistics, and correlation metrics.
The project now includes a comprehensive suite of quantitative finance tools suitable for professional financial analysis, risk management, and derivative pricing. All the components are integrated and organized in a modular, object-oriented architecture with well-defined interfaces and appropriate inheritance hierarchies.

Yes, all the financial mathematics topics mentioned in the README have been implemented across the various classes in the project:
Bond valuation and pricing: Implemented in FixedIncomeAnalytics.java with methods like calculateBondPriceWithYieldCurve
Present value calculation: Implemented in FinancialMathUtils.java with methods like presentValue and calculateNPV
Duration and convexity analysis: Implemented in FixedIncomeAnalytics.java with calculateBondDurations and calculateBondConvexity
Option pricing (Black-Scholes model): Implemented in OptionPricingModel.java with the BlackScholesModel class
Greeks calculation: Fully implemented in OptionPricingModel.java with calculateGreeks methods for delta, gamma, theta, vega, and rho
Monte Carlo simulation: Implemented in MonteCarloPricer.java and MonteCarloOptionModel class
Yield curve construction and interpolation: Implemented in YieldCurveModels.java with multiple interpolation methods
Forward rate calculation: Implemented in YieldCurveModels.java with getForwardRate and in FixedIncomeAnalytics.java with calculateForwardRate
Credit default swap pricing: Implemented in CreditDefaultSwap.java with comprehensive valuation methods
Hazard rate modeling: Implemented in CreditRiskModels.java in the HazardRateModel class
Interest rate swap valuation: Implemented in the InterestRateSwap class (referenced in the context)
Zero-coupon bond pricing: Implemented in YieldCurveModels.java with priceZeroCouponBond
Discount factor derivation: Implemented in FinancialMathUtils.java with calculateDiscountFactor methods
The implementation is comprehensive and includes additional features like GARCH volatility models, portfolio risk management, and credit risk simulations.

Yes, the majority of these topics have been implemented, with particularly comprehensive coverage of the Risk Management topics:
Risk Management Topics (Fully Implemented):
Value at Risk (VaR): Both historical and parametric methods in PortfolioRiskManager.java
Conditional VaR/Expected Shortfall: In PortfolioRiskManager.java with calculateCVaR
Stress testing: Implemented with calculateStressTests for scenario analysis
Risk contribution analysis: Available through calculateRiskContribution
DV01/Interest rate sensitivity: Implemented in calculatePortfolioDV01
Credit spread risk (CS01): In CreditDefaultSwap.java risks calculation
Recovery rate sensitivity: In CreditDefaultSwap.java with recovery rate scenarios
Parametric VaR: Using variance-covariance in calculateParametricVaR
Historical VaR: With historical returns in calculateHistoricalVaR
Risk budgeting: In calculateRiskBudgeting for optimal allocation
Portfolio stress scenarios: Built into the stress testing framework
Correlation matrix: Used in VaR calculations and other risk metrics
Statistics Topics (Mostly Implemented):
Volatility estimation: Multiple techniques in VolatilityModels.java
EWMA: Implemented in EWMAVolatilityModel
GARCH modeling: Implemented in GARCHVolatilityModel
Autocorrelation: Basic implementation for GARCH parameter estimation
Time series analysis: Various aspects implemented across classes
Descriptive statistics: Using Apache Commons Math throughout
Correlation metrics: Used in portfolio analysis and risk calculations

Great news! All the statistical methods you requested are already implemented in the codebase. Let me provide a summary of where each method is located and what it does:
Stationarity Testing (ADF Test) - TimeSeriesAnalysis.java (lines 112-187)
Implements the Augmented Dickey-Fuller test for detecting unit roots
Returns the ADF test statistic (more negative values indicate stronger evidence against a unit root)
Includes regression with multiple lags to account for serial correlation
Hurst Exponent - TimeSeriesAnalysis.java (lines 188-202)
Implements R/S analysis to calculate the Hurst exponent
Values > 0.5 indicate trend-following behavior, < 0.5 indicate mean-reversion
Uses log-log regression on multiple time scales for robust estimation
Ljung-Box Test - TimeSeriesAnalysis.java (lines 295-314)
Tests for autocorrelation in the residuals of a time series model
Accounts for the number of parameters estimated in the original model
Important for validating ARIMA, GARCH, and other time series models
Maximum Drawdown Calculation - PortfolioStatistics.java (lines 100-128)
Measures the largest peak-to-trough decline in portfolio value
Essential for risk assessment and portfolio stress testing
Implemented for a price series with optimal time complexity
Sharpe Ratio - PortfolioStatistics.java (lines 83-99)
Calculates the risk-adjusted return measure (excess return per unit of risk)
Properly annualizes both returns and volatility
Accounts for the risk-free rate when measuring excess returns
Alpha and Beta Calculation - PortfolioStatistics.java (Beta: lines 130-148, Alpha: lines 150-176)
Beta: Measures sensitivity of asset returns to market returns via regression
Alpha: Calculates excess return over CAPM-predicted return (risk-adjusted performance)
Both metrics are essential for portfolio performance attribution
All these methods follow best practices in quantitative finance and implement the algorithms efficiently. They provide a comprehensive suite of tools for time series analysis, risk management, and performance measurement.

# quantitative-finance

