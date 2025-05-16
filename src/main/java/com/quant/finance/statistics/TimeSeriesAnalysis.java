package com.quant.finance.statistics;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Class for financial time series analysis including autocorrelation, stationarity,
 * and other time series specific statistical measures.
 */
public class TimeSeriesAnalysis {
    
    /**
     * Calculate autocorrelation of a time series at a specific lag
     * 
     * @param timeSeries Array of return data
     * @param lag Lag period (e.g., 1 for lag-1 autocorrelation)
     * @return Autocorrelation coefficient
     */
    public static double calculateAutocorrelation(double[] timeSeries, int lag) {
        if (lag >= timeSeries.length) {
            throw new IllegalArgumentException("Lag cannot be greater than or equal to series length");
        }
        
        int n = timeSeries.length;
        double[] series1 = Arrays.copyOfRange(timeSeries, 0, n - lag);
        double[] series2 = Arrays.copyOfRange(timeSeries, lag, n);
        
        PearsonsCorrelation pearson = new PearsonsCorrelation();
        return pearson.correlation(series1, series2);
    }
    
    /**
     * Calculate autocorrelation function (ACF) for multiple lags
     * 
     * @param timeSeries Array of return data
     * @param maxLag Maximum lag to calculate
     * @return Array of autocorrelation coefficients for lags 1 to maxLag
     */
    public static double[] calculateACF(double[] timeSeries, int maxLag) {
        double[] acf = new double[maxLag];
        
        for (int lag = 1; lag <= maxLag; lag++) {
            acf[lag - 1] = calculateAutocorrelation(timeSeries, lag);
        }
        
        return acf;
    }
    
    /**
     * Calculate partial autocorrelation function (PACF)
     * Implementation of the Durbin-Levinson algorithm
     * 
     * @param timeSeries Array of return data
     * @param maxLag Maximum lag to calculate
     * @return Array of partial autocorrelation coefficients for lags 1 to maxLag
     */
    public static double[] calculatePACF(double[] timeSeries, int maxLag) {
        double[] acf = calculateACF(timeSeries, maxLag);
        double[] pacf = new double[maxLag];
        
        // Allocate workspace
        double[][] phi = new double[maxLag][maxLag];
        
        // Initialize PACF at lag 1
        pacf[0] = acf[0];
        phi[0][0] = acf[0];
        
        // Calculate PACF for lags 2 to maxLag using Durbin-Levinson recursion
        for (int k = 1; k < maxLag; k++) {
            // Calculate numerator
            double numerator = acf[k];
            for (int j = 0; j < k; j++) {
                numerator -= phi[k-1][j] * acf[k-j-1];
            }
            
            // Calculate phi_kk (the PACF value at lag k+1)
            phi[k][k] = numerator;
            
            // Calculate denominator
            double denominator = 1.0;
            for (int j = 0; j < k; j++) {
                denominator -= phi[k-1][j] * acf[j];
            }
            
            // Complete the calculation of phi_kk
            phi[k][k] /= denominator;
            pacf[k] = phi[k][k];
            
            // Update the remaining phi values for this lag
            for (int j = 0; j < k; j++) {
                phi[k][j] = phi[k-1][j] - phi[k][k] * phi[k-1][k-j-1];
            }
        }
        
        return pacf;
    }
    
    /**
     * Calculate the Augmented Dickey-Fuller test statistic for stationarity testing
     * Implementation of ADF test without trend
     * 
     * @param timeSeries Array of price data (not returns)
     * @param lags Number of lags to include in the regression
     * @return ADF test statistic (more negative indicates more evidence against unit root)
     */
    public static double adfTest(double[] timeSeries, int lags) {
        int n = timeSeries.length;
        if (n <= lags + 2) {
            throw new IllegalArgumentException("Time series must be longer than lags+2");
        }
        
        // Calculate first differences
        double[] diff = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            diff[i] = timeSeries[i + 1] - timeSeries[i];
        }
        
        // Set up regression predictors and response
        // Response: diff[t]
        // Predictors: timeSeries[t-1], diff[t-1], diff[t-2], ..., diff[t-lags]
        int regressionLength = n - lags - 1;
        double[][] x = new double[regressionLength][lags + 1];
        double[] y = new double[regressionLength];
        
        for (int t = 0; t < regressionLength; t++) {
            // Level term (y_{t-1})
            x[t][0] = timeSeries[t + lags - 1];
            
            // Lagged difference terms
            for (int j = 1; j <= lags; j++) {
                x[t][j] = diff[t + lags - j];
            }
            
            // Response (y_t)
            y[t] = diff[t + lags];
        }
        
        // Perform the regression
        SimpleRegression regression = new SimpleRegression(true);
        
        // Extract the coefficient of the level term (delta) and its standard error
        double delta = 0.0;
        double sumX = 0.0;
        double sumY = 0.0;
        double sumXY = 0.0;
        double sumX2 = 0.0;
        
        for (int i = 0; i < regressionLength; i++) {
            sumX += x[i][0];
            sumY += y[i];
            sumXY += x[i][0] * y[i];
            sumX2 += x[i][0] * x[i][0];
        }
        
        delta = (regressionLength * sumXY - sumX * sumY) / (regressionLength * sumX2 - sumX * sumX);
        
        // Calculate residuals
        double[] residuals = new double[regressionLength];
        for (int i = 0; i < regressionLength; i++) {
            residuals[i] = y[i] - delta * x[i][0];
        }
        
        // Calculate standard error of delta
        double sumResiduals2 = 0.0;
        for (double residual : residuals) {
            sumResiduals2 += residual * residual;
        }
        
        double se = Math.sqrt(sumResiduals2 / (regressionLength - 2)) / Math.sqrt(sumX2 - sumX * sumX / regressionLength);
        
        // Calculate test statistic
        return delta / se;
    }
    
    /**
     * Calculate Hurst exponent using R/S analysis
     * Used to detect long-term memory in time series
     * 
     * @param timeSeries Array of return data
     * @return Hurst exponent (0.5 = random walk, >0.5 = trending, <0.5 = mean reverting)
     */
    public static double calculateHurstExponent(double[] timeSeries) {
        int n = timeSeries.length;
        
        // Need a sufficient number of data points
        if (n < 100) {
            throw new IllegalArgumentException("Time series should have at least 100 data points");
        }
        
        // Calculate various window sizes
        List<Integer> tauValues = new ArrayList<>();
        for (int i = 10; i <= n/10; i += 10) {
            tauValues.add(i);
        }
        
        // Calculate R/S values for each window size
        double[] logTau = new double[tauValues.size()];
        double[] logRS = new double[tauValues.size()];
        
        for (int i = 0; i < tauValues.size(); i++) {
            int tau = tauValues.get(i);
            double rs = calculateRSValue(timeSeries, tau);
            
            logTau[i] = Math.log10(tau);
            logRS[i] = Math.log10(rs);
        }
        
        // Linear regression to find Hurst exponent
        SimpleRegression regression = new SimpleRegression();
        
        for (int i = 0; i < logTau.length; i++) {
            regression.addData(logTau[i], logRS[i]);
        }
        
        return regression.getSlope();
    }
    
    /**
     * Calculate R/S value for a specific window size
     * 
     * @param timeSeries Time series data
     * @param windowSize Size of the window
     * @return R/S value
     */
    private static double calculateRSValue(double[] timeSeries, int windowSize) {
        int n = timeSeries.length;
        int numWindows = n / windowSize;
        
        double sumRS = 0.0;
        
        for (int w = 0; w < numWindows; w++) {
            double[] windowData = Arrays.copyOfRange(timeSeries, w * windowSize, (w + 1) * windowSize);
            
            // Calculate mean
            double mean = 0.0;
            for (double value : windowData) {
                mean += value;
            }
            mean /= windowSize;
            
            // Calculate cumulative deviations
            double[] deviations = new double[windowSize];
            for (int i = 0; i < windowSize; i++) {
                deviations[i] = windowData[i] - mean;
            }
            
            double[] cumulativeDeviations = new double[windowSize];
            cumulativeDeviations[0] = deviations[0];
            for (int i = 1; i < windowSize; i++) {
                cumulativeDeviations[i] = cumulativeDeviations[i - 1] + deviations[i];
            }
            
            // Calculate range (max - min of cumulative deviations)
            double max = cumulativeDeviations[0];
            double min = cumulativeDeviations[0];
            
            for (double deviation : cumulativeDeviations) {
                if (deviation > max) max = deviation;
                if (deviation < min) min = deviation;
            }
            
            double range = max - min;
            
            // Calculate standard deviation
            double sumSquaredDeviations = 0.0;
            for (double deviation : deviations) {
                sumSquaredDeviations += deviation * deviation;
            }
            
            double stdDev = Math.sqrt(sumSquaredDeviations / windowSize);
            
            // Calculate R/S and add to sum
            double rs = range / stdDev;
            sumRS += rs;
        }
        
        // Return average R/S value
        return sumRS / numWindows;
    }
    
    /**
     * Calculate the Ljung-Box test statistic for autocorrelation in residuals
     * 
     * @param residuals Array of residuals (e.g., from an ARIMA model)
     * @param lags Number of lags to include in the test
     * @param numParams Number of parameters estimated in the model
     * @return Ljung-Box test statistic
     */
    public static double ljungBoxTest(double[] residuals, int lags, int numParams) {
        int n = residuals.length;
        double[] autocorr = calculateACF(residuals, lags);
        
        double sumTerm = 0.0;
        for (int i = 0; i < lags; i++) {
            sumTerm += (autocorr[i] * autocorr[i]) / (n - i - 1);
        }
        
        return n * (n + 2) * sumTerm;
    }
    
    /**
     * Calculate volatility using exponentially weighted moving average (EWMA)
     * 
     * @param returns Array of return data
     * @param lambda Decay factor (typical value: 0.94 for daily returns)
     * @return Array of EWMA volatility estimates
     */
    public static double[] calculateEWMAVolatility(double[] returns, double lambda) {
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
    
    /**
     * Calculate GARCH(1,1) volatility forecast
     * Simplified GARCH model with parameters alpha, beta, and omega
     * 
     * @param returns Array of historical returns
     * @param alpha ARCH parameter (sensitivity to recent returns)
     * @param beta GARCH parameter (persistence of volatility)
     * @param omega Long-run variance parameter
     * @return Array of GARCH volatility estimates
     */
    public static double[] calculateGARCHVolatility(double[] returns, double alpha, double beta, double omega) {
        if (alpha + beta >= 1.0) {
            throw new IllegalArgumentException("alpha + beta must be less than 1 for stationarity");
        }
        
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
    
    /**
     * Calculate log-likelihood for GARCH(1,1) model
     * Used for parameter estimation
     * 
     * @param returns Array of returns
     * @param alpha ARCH parameter
     * @param beta GARCH parameter
     * @param omega Long-run variance parameter
     * @return Log-likelihood value
     */
    public static double calculateGARCHLogLikelihood(double[] returns, double alpha, double beta, double omega) {
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
            logLikelihood += -0.5 * (Math.log(2 * Math.PI) + Math.log(variance[i]) + returns[i] * returns[i] / variance[i]);
        }
        
        return logLikelihood;
    }
} 