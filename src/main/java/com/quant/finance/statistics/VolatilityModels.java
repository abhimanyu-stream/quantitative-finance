package com.quant.finance.statistics;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * Abstract base class for volatility models
 */
public abstract class VolatilityModels {
    
    protected double[] returns;
    
    /**
     * Constructor
     * 
     * @param returns Array of return data
     */
    public VolatilityModels(double[] returns) {
        this.returns = Arrays.copyOf(returns, returns.length);
    }
    
    /**
     * Calculate volatility estimates
     * 
     * @return Array of volatility estimates
     */
    public abstract double[] calculateVolatility();
    
    /**
     * Calculate forecast volatility for the next period
     * 
     * @return Forecasted volatility
     */
    public abstract double forecastVolatility();
    
    /**
     * Create a volatility model of the specified type
     * 
     * @param returns Array of return data
     * @param modelType Type of volatility model
     * @return VolatilityModels instance
     */
    public static VolatilityModels createModel(double[] returns, String modelType) {
        switch (modelType.toLowerCase()) {
            case "historical":
                return new HistoricalVolatilityModel(returns);
            case "ewma":
                return new EWMAVolatilityModel(returns);
            case "garch":
                return new GARCHVolatilityModel(returns);
            default:
                throw new IllegalArgumentException("Unknown volatility model type: " + modelType);
        }
    }
}

/**
 * Historical volatility model
 * Uses a rolling window of returns to calculate volatility
 */
class HistoricalVolatilityModel extends VolatilityModels {
    private int windowSize = 30; // Default window size
    
    public HistoricalVolatilityModel(double[] returns) {
        super(returns);
    }
    
    public HistoricalVolatilityModel(double[] returns, int windowSize) {
        super(returns);
        this.windowSize = windowSize;
    }
    
    @Override
    public double[] calculateVolatility() {
        if (returns.length < windowSize) {
            throw new IllegalArgumentException("Not enough data points for window size");
        }
        
        double[] volatility = new double[returns.length];
        
        // For the first (windowSize-1) data points, use the volatility of the first window
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = 0; i < windowSize; i++) {
            stats.addValue(returns[i]);
        }
        
        double firstWindowVol = Math.sqrt(stats.getVariance());
        for (int i = 0; i < windowSize; i++) {
            volatility[i] = firstWindowVol;
        }
        
        // Calculate rolling window volatility
        for (int i = windowSize; i < returns.length; i++) {
            stats.addValue(returns[i]);
            stats.removeMostRecentValue();
            volatility[i] = Math.sqrt(stats.getVariance());
        }
        
        return volatility;
    }
    
    @Override
    public double forecastVolatility() {
        // For historical volatility, the forecast is simply the most recent volatility
        if (returns.length < windowSize) {
            throw new IllegalArgumentException("Not enough data points for window size");
        }
        
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = returns.length - windowSize; i < returns.length; i++) {
            stats.addValue(returns[i]);
        }
        
        return Math.sqrt(stats.getVariance());
    }
}

/**
 * Exponentially Weighted Moving Average (EWMA) volatility model
 */
class EWMAVolatilityModel extends VolatilityModels {
    private double lambda = 0.94; // Default decay factor for daily returns
    
    public EWMAVolatilityModel(double[] returns) {
        super(returns);
    }
    
    public EWMAVolatilityModel(double[] returns, double lambda) {
        super(returns);
        this.lambda = lambda;
    }
    
    @Override
    public double[] calculateVolatility() {
        int n = returns.length;
        double[] variance = new double[n];
        
        // Initialize with sample variance of first few observations
        DescriptiveStatistics stats = new DescriptiveStatistics(Arrays.copyOfRange(returns, 0, Math.min(10, n)));
        variance[0] = stats.getVariance();
        
        // Calculate EWMA variance
        for (int i = 1; i < n; i++) {
            variance[i] = lambda * variance[i - 1] + (1 - lambda) * returns[i - 1] * returns[i - 1];
        }
        
        // Convert variance to volatility (standard deviation)
        double[] volatility = new double[n];
        for (int i = 0; i < n; i++) {
            volatility[i] = Math.sqrt(variance[i]);
        }
        
        return volatility;
    }
    
    @Override
    public double forecastVolatility() {
        // For EWMA, forecast is based on the most recent variance estimate
        double[] volatility = calculateVolatility();
        return volatility[volatility.length - 1];
    }
    
    /**
     * Set the lambda decay factor
     * 
     * @param lambda Decay factor (e.g., 0.94 for daily returns, 0.97 for monthly)
     */
    public void setLambda(double lambda) {
        if (lambda <= 0 || lambda >= 1) {
            throw new IllegalArgumentException("Lambda must be between 0 and 1");
        }
        this.lambda = lambda;
    }
}

/**
 * GARCH(1,1) volatility model
 */
class GARCHVolatilityModel extends VolatilityModels {
    private double omega = 0.000001;  // Long-term variance parameter
    private double alpha = 0.1;       // ARCH parameter (weight of past returns)
    private double beta = 0.85;       // GARCH parameter (persistence)
    
    public GARCHVolatilityModel(double[] returns) {
        super(returns);
        estimateParameters();
    }
    
    public GARCHVolatilityModel(double[] returns, double omega, double alpha, double beta) {
        super(returns);
        this.omega = omega;
        this.alpha = alpha;
        this.beta = beta;
        
        if (alpha + beta >= 1.0) {
            throw new IllegalArgumentException("alpha + beta must be less than 1 for stationarity");
        }
    }
    
    /**
     * Estimate GARCH parameters using a simplified approach
     */
    private void estimateParameters() {
        // Use default parameters if less than 100 data points
        if (returns.length < 100) {
            return;
        }
        
        // Simplified parameter estimation approach
        // Instead of using the optimization framework, we'll use reasonable defaults
        // and then adjust based on the data characteristics
        
        // Calculate sample variance
        DescriptiveStatistics stats = new DescriptiveStatistics(returns);
        double sampleVariance = stats.getVariance();
        
        // Calculate first-order autocorrelation of squared returns
        double[] squaredReturns = new double[returns.length];
        for (int i = 0; i < returns.length; i++) {
            squaredReturns[i] = returns[i] * returns[i];
        }
        
        double meanSquaredReturns = stats.getMean();
        double autocorr = calculateAutocorrelation(squaredReturns, 1);
        
        // Adjust parameters based on autocorrelation
        // Higher autocorrelation suggests more persistent volatility (higher beta)
        if (autocorr > 0.2) {
            beta = 0.85 + 0.1 * (autocorr - 0.2); // Increase beta for higher autocorrelation
            beta = Math.min(beta, 0.98); // Cap at 0.98
        } else {
            beta = 0.7 + 0.75 * autocorr; // Lower beta for lower autocorrelation
        }
        
        // Set alpha based on beta, ensuring stationarity
        alpha = Math.min(0.15, 0.98 - beta);
        
        // Set omega to target the unconditional variance (sample variance)
        omega = sampleVariance * (1.0 - alpha - beta);
        
        // Ensure omega is positive
        omega = Math.max(omega, 1e-6);
    }
    
    /**
     * Calculate autocorrelation with given lag
     */
    private double calculateAutocorrelation(double[] series, int lag) {
        if (lag >= series.length) {
            return 0.0;
        }
        
        // Calculate mean
        double mean = 0.0;
        for (double value : series) {
            mean += value;
        }
        mean /= series.length;
        
        // Calculate autocorrelation
        double numerator = 0.0;
        double denominator = 0.0;
        
        for (int i = lag; i < series.length; i++) {
            numerator += (series[i] - mean) * (series[i - lag] - mean);
        }
        
        for (int i = 0; i < series.length; i++) {
            denominator += Math.pow(series[i] - mean, 2);
        }
        
        return numerator / denominator;
    }
    
    /**
     * Calculate log-likelihood for GARCH model
     */
    private double calculateLogLikelihood(double omega, double alpha, double beta) {
        int n = returns.length;
        double[] variance = new double[n];
        double logLikelihood = 0.0;
        
        // Initialize variance
        double unconditionalVariance = omega / (1 - alpha - beta);
        variance[0] = unconditionalVariance;
        
        // Calculate GARCH variance and log-likelihood
        for (int i = 1; i < n; i++) {
            // Update variance
            variance[i] = omega + alpha * returns[i - 1] * returns[i - 1] + beta * variance[i - 1];
            
            // Add term to log-likelihood
            logLikelihood += -0.5 * (Math.log(2 * Math.PI) + Math.log(variance[i]) + 
                                    returns[i] * returns[i] / variance[i]);
        }
        
        return logLikelihood;
    }
    
    @Override
    public double[] calculateVolatility() {
        int n = returns.length;
        double[] variance = new double[n];
        
        // Initialize with unconditional variance
        double unconditionalVariance = omega / (1 - alpha - beta);
        variance[0] = unconditionalVariance;
        
        // Calculate GARCH variance recursively
        for (int i = 1; i < n; i++) {
            variance[i] = omega + alpha * returns[i - 1] * returns[i - 1] + beta * variance[i - 1];
        }
        
        // Convert variance to volatility
        double[] volatility = new double[n];
        for (int i = 0; i < n; i++) {
            volatility[i] = Math.sqrt(variance[i]);
        }
        
        return volatility;
    }
    
    @Override
    public double forecastVolatility() {
        // Calculate current variance
        double[] variance = new double[returns.length];
        
        // Initialize with unconditional variance
        double unconditionalVariance = omega / (1 - alpha - beta);
        variance[0] = unconditionalVariance;
        
        // Calculate GARCH variance recursively
        for (int i = 1; i < returns.length; i++) {
            variance[i] = omega + alpha * returns[i - 1] * returns[i - 1] + beta * variance[i - 1];
        }
        
        // Forecast one-step ahead
        double forecastVariance = omega + alpha * returns[returns.length - 1] * returns[returns.length - 1] + 
                                 beta * variance[variance.length - 1];
        
        return Math.sqrt(forecastVariance);
    }
    
    /**
     * Create a Monte Carlo simulation of future returns based on GARCH model
     * 
     * @param horizon Number of periods to forecast
     * @param numSimulations Number of simulations
     * @return 2D array of simulated returns [simulation][time]
     */
    public double[][] simulateReturns(int horizon, int numSimulations) {
        double[][] simulations = new double[numSimulations][horizon];
        Random random = new Random();
        
        // Calculate current variance
        double[] variance = new double[returns.length];
        double unconditionalVariance = omega / (1 - alpha - beta);
        variance[0] = unconditionalVariance;
        
        for (int i = 1; i < returns.length; i++) {
            variance[i] = omega + alpha * returns[i - 1] * returns[i - 1] + beta * variance[i - 1];
        }
        
        // Last variance is our starting point
        double currentVariance = variance[variance.length - 1];
        
        // For each simulation
        for (int s = 0; s < numSimulations; s++) {
            double simulatedVariance = currentVariance;
            
            // For each time step
            for (int t = 0; t < horizon; t++) {
                // Generate random return from normal distribution
                double z = random.nextGaussian();
                double simulatedReturn = z * Math.sqrt(simulatedVariance);
                simulations[s][t] = simulatedReturn;
                
                // Update variance for next step
                simulatedVariance = omega + alpha * simulatedReturn * simulatedReturn + beta * simulatedVariance;
            }
        }
        
        return simulations;
    }
    
    /**
     * Get the parameters of the GARCH model
     * 
     * @return Map of parameter names to values
     */
    public Map<String, Double> getParameters() {
        Map<String, Double> params = new HashMap<>();
        params.put("omega", omega);
        params.put("alpha", alpha);
        params.put("beta", beta);
        params.put("persistence", alpha + beta);
        params.put("longRunVariance", omega / (1 - alpha - beta));
        return params;
    }
} 