package com.quant.finance.statistics;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Statistical utility class for financial portfolio analysis
 */
public class PortfolioStatistics {

    /**
     * Calculate annualized returns from a sequence of prices
     * 
     * @param prices List of sequential asset prices
     * @param periodsPerYear Number of periods in a year (e.g., 252 for daily trading days)
     * @return Annualized return as a decimal (e.g., 0.08 for 8%)
     */
    public static double calculateAnnualizedReturn(List<Double> prices, int periodsPerYear) {
        if (prices.size() < 2) {
            throw new IllegalArgumentException("Need at least two prices to calculate returns");
        }
        
        double firstPrice = prices.get(0);
        double lastPrice = prices.get(prices.size() - 1);
        int periods = prices.size() - 1;
        
        // Calculate total return
        double totalReturn = lastPrice / firstPrice - 1;
        
        // Annualize the return
        double annualizedReturn = Math.pow(1 + totalReturn, periodsPerYear / (double) periods) - 1;
        
        return annualizedReturn;
    }
    
    /**
     * Calculate a series of logarithmic returns from a list of prices
     * 
     * @param prices List of sequential asset prices
     * @return List of logarithmic returns
     */
    public static double[] calculateLogReturns(List<Double> prices) {
        if (prices.size() < 2) {
            throw new IllegalArgumentException("Need at least two prices to calculate returns");
        }
        
        double[] logReturns = new double[prices.size() - 1];
        
        for (int i = 1; i < prices.size(); i++) {
            double currentPrice = prices.get(i);
            double previousPrice = prices.get(i - 1);
            
            logReturns[i - 1] = Math.log(currentPrice / previousPrice);
        }
        
        return logReturns;
    }
    
    /**
     * Calculate the annualized volatility (standard deviation) of returns
     * 
     * @param returns Array of periodic returns (not prices)
     * @param periodsPerYear Number of periods in a year (e.g., 252 for daily trading days)
     * @return Annualized volatility
     */
    public static double calculateAnnualizedVolatility(double[] returns, int periodsPerYear) {
        DescriptiveStatistics stats = new DescriptiveStatistics(returns);
        double stdDev = stats.getStandardDeviation();
        
        // Annualized volatility
        return stdDev * Math.sqrt(periodsPerYear);
    }
    
    /**
     * Calculate the Sharpe ratio of an investment
     * 
     * @param returns Array of periodic returns
     * @param riskFreeRate Risk-free rate for the period as a decimal (e.g., 0.02 for 2%)
     * @param periodsPerYear Number of periods in a year
     * @return Sharpe ratio
     */
    public static double calculateSharpeRatio(double[] returns, double riskFreeRate, int periodsPerYear) {
        DescriptiveStatistics stats = new DescriptiveStatistics(returns);
        double meanReturn = stats.getMean();
        double stdDev = stats.getStandardDeviation();
        
        // Annualize the excess return
        double excessReturn = (meanReturn - (riskFreeRate / periodsPerYear));
        double annualizedExcessReturn = excessReturn * periodsPerYear;
        
        // Annualize the volatility
        double annualizedVolatility = stdDev * Math.sqrt(periodsPerYear);
        
        // Calculate Sharpe ratio
        return annualizedExcessReturn / annualizedVolatility;
    }
    
    /**
     * Calculate maximum drawdown from a price series
     * 
     * @param prices List of sequential asset prices
     * @return Maximum drawdown as a positive decimal (e.g., 0.2 for 20% drawdown)
     */
    public static double calculateMaxDrawdown(List<Double> prices) {
        if (prices.isEmpty()) {
            return 0.0;
        }
        
        double maxDrawdown = 0.0;
        double peak = prices.get(0);
        
        for (Double price : prices) {
            if (price > peak) {
                peak = price;
            }
            
            double drawdown = (peak - price) / peak;
            if (drawdown > maxDrawdown) {
                maxDrawdown = drawdown;
            }
        }
        
        return maxDrawdown;
    }
    
    /**
     * Calculate Beta of an asset against a benchmark
     * 
     * @param assetReturns Array of asset returns
     * @param benchmarkReturns Array of benchmark returns (same length as assetReturns)
     * @return Beta value
     */
    public static double calculateBeta(double[] assetReturns, double[] benchmarkReturns) {
        if (assetReturns.length != benchmarkReturns.length) {
            throw new IllegalArgumentException("Asset and benchmark return series must be the same length");
        }
        
        SimpleRegression regression = new SimpleRegression();
        
        for (int i = 0; i < assetReturns.length; i++) {
            regression.addData(benchmarkReturns[i], assetReturns[i]);
        }
        
        return regression.getSlope();
    }
    
    /**
     * Calculate Alpha of an asset (using CAPM)
     * 
     * @param assetReturns Array of asset returns
     * @param benchmarkReturns Array of benchmark returns
     * @param riskFreeRate Risk-free rate for the period
     * @param periodsPerYear Number of periods in a year
     * @return Alpha value (annualized)
     */
    public static double calculateAlpha(double[] assetReturns, double[] benchmarkReturns, 
                                        double riskFreeRate, int periodsPerYear) {
        double beta = calculateBeta(assetReturns, benchmarkReturns);
        
        DescriptiveStatistics assetStats = new DescriptiveStatistics(assetReturns);
        DescriptiveStatistics benchmarkStats = new DescriptiveStatistics(benchmarkReturns);
        
        double meanAssetReturn = assetStats.getMean();
        double meanBenchmarkReturn = benchmarkStats.getMean();
        
        // Calculate period risk-free rate
        double periodRiskFreeRate = riskFreeRate / periodsPerYear;
        
        // Calculate alpha for the period
        double alpha = meanAssetReturn - (periodRiskFreeRate + beta * (meanBenchmarkReturn - periodRiskFreeRate));
        
        // Annualize alpha
        return alpha * periodsPerYear;
    }
    
    /**
     * Calculate the correlation matrix for a list of asset return series
     * 
     * @param returnsList List of return arrays for multiple assets
     * @param assetNames List of asset names (must match the order of returnsList)
     * @return Correlation matrix as a Map with keys formatted as "asset1,asset2"
     */
    public static Map<String, Double> calculateCorrelationMatrix(List<double[]> returnsList, List<String> assetNames) {
        if (returnsList.size() != assetNames.size()) {
            throw new IllegalArgumentException("Number of return series must match number of asset names");
        }
        
        Map<String, Double> correlationMatrix = new HashMap<>();
        PearsonsCorrelation correlation = new PearsonsCorrelation();
        
        for (int i = 0; i < returnsList.size(); i++) {
            for (int j = i; j < returnsList.size(); j++) {
                double correlationValue = correlation.correlation(returnsList.get(i), returnsList.get(j));
                String key = assetNames.get(i) + "," + assetNames.get(j);
                correlationMatrix.put(key, correlationValue);
                
                // Add the symmetric pair if different assets
                if (i != j) {
                    String symmetricKey = assetNames.get(j) + "," + assetNames.get(i);
                    correlationMatrix.put(symmetricKey, correlationValue);
                }
            }
        }
        
        return correlationMatrix;
    }
    
    /**
     * Calculate the Value at Risk (VaR) using historical method
     * 
     * @param returns Array of historical returns
     * @param confidenceLevel Confidence level (e.g., 0.95 for 95%)
     * @param investmentValue Current value of the investment
     * @return Value at Risk (a positive number representing potential loss)
     */
    public static double calculateHistoricalVaR(double[] returns, double confidenceLevel, double investmentValue) {
        // Sort returns in ascending order
        double[] sortedReturns = Arrays.copyOf(returns, returns.length);
        Arrays.sort(sortedReturns);
        
        // Find the index corresponding to the confidence level
        int index = (int) Math.floor((1 - confidenceLevel) * sortedReturns.length);
        
        // Get the return at that index
        double varReturn = sortedReturns[index];
        
        // Convert to a loss amount (VaR is expressed as a positive number)
        return -varReturn * investmentValue;
    }
    
    /**
     * Calculate the Conditional Value at Risk (CVaR) using historical method
     * Also known as Expected Shortfall
     * 
     * @param returns Array of historical returns
     * @param confidenceLevel Confidence level (e.g., 0.95 for 95%)
     * @param investmentValue Current value of the investment
     * @return Conditional Value at Risk (a positive number representing potential loss)
     */
    public static double calculateHistoricalCVaR(double[] returns, double confidenceLevel, double investmentValue) {
        // Sort returns in ascending order
        double[] sortedReturns = Arrays.copyOf(returns, returns.length);
        Arrays.sort(sortedReturns);
        
        // Find the index corresponding to the confidence level
        int index = (int) Math.floor((1 - confidenceLevel) * sortedReturns.length);
        
        // Calculate the average of returns below the VaR threshold
        double sum = 0.0;
        for (int i = 0; i < index; i++) {
            sum += sortedReturns[i];
        }
        
        double averageWorstReturns = sum / index;
        
        // Convert to a loss amount (CVaR is expressed as a positive number)
        return -averageWorstReturns * investmentValue;
    }
} 