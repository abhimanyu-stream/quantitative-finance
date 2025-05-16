package com.quant.finance.risk;

import com.quant.finance.utils.DateUtils;
import com.quant.finance.utils.FinancialMathUtils;

import java.time.LocalDate;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * Abstract base class for credit risk models
 */
public abstract class CreditRiskModels {
    
    /**
     * Calculate probability of default over a time horizon
     * 
     * @param timeHorizon Time horizon in years
     * @param parameters Model parameters
     * @return Probability of default
     */
    public abstract double calculateProbabilityOfDefault(double timeHorizon, Map<String, Object> parameters);
    
    /**
     * Calculate loss given default
     * 
     * @param exposureAtDefault Exposure at default
     * @param recoveryRate Recovery rate (as a decimal)
     * @return Loss given default
     */
    public double calculateLossGivenDefault(double exposureAtDefault, double recoveryRate) {
        return exposureAtDefault * (1.0 - recoveryRate);
    }
    
    /**
     * Calculate expected loss
     * 
     * @param exposureAtDefault Exposure at default
     * @param probabilityOfDefault Probability of default
     * @param recoveryRate Recovery rate (as a decimal)
     * @return Expected loss
     */
    public double calculateExpectedLoss(double exposureAtDefault, double probabilityOfDefault, double recoveryRate) {
        double lossGivenDefault = calculateLossGivenDefault(exposureAtDefault, recoveryRate);
        return probabilityOfDefault * lossGivenDefault;
    }
    
    /**
     * Factory method to create a credit risk model of the specified type
     * 
     * @param modelType Type of credit risk model
     * @return CreditRiskModels instance
     */
    public static CreditRiskModels createModel(String modelType) {
        switch (modelType.toLowerCase()) {
            case "merton":
                return new MertonModel();
            case "hazardrate":
                return new HazardRateModel();
            case "creditmetrics":
                return new CreditMetricsModel();
            default:
                throw new IllegalArgumentException("Unknown credit risk model type: " + modelType);
        }
    }
    
    /**
     * Calculate the standard normal cumulative distribution function
     * 
     * @param x Input value
     * @return Standard normal CDF value
     */
    protected double standardNormalCDF(double x) {
        // Use an approximation of the standard normal CDF
        // Abramowitz and Stegun approximation 7.1.26
        double b1 = 0.319381530;
        double b2 = -0.356563782;
        double b3 = 1.781477937;
        double b4 = -1.821255978;
        double b5 = 1.330274429;
        double p = 0.2316419;
        
        double t = 1.0 / (1.0 + p * Math.abs(x));
        double standardNormalPDF = Math.exp(-x * x / 2.0) / Math.sqrt(2.0 * Math.PI);
        
        double cumulativeProbability = 1.0 - standardNormalPDF * 
                                     (b1 * t + 
                                      b2 * t * t + 
                                      b3 * t * t * t + 
                                      b4 * t * t * t * t + 
                                      b5 * t * t * t * t * t);
                                      
        return (x < 0) ? 1.0 - cumulativeProbability : cumulativeProbability;
    }
}

/**
 * Merton model for credit risk (structural model)
 * Based on option pricing theory to model default
 */
class MertonModel extends CreditRiskModels {
    
    @Override
    public double calculateProbabilityOfDefault(double timeHorizon, Map<String, Object> parameters) {
        // Extract parameters
        double assetValue = (double) parameters.getOrDefault("assetValue", 0.0);
        double assetVolatility = (double) parameters.getOrDefault("assetVolatility", 0.0);
        double debtValue = (double) parameters.getOrDefault("debtValue", 0.0);
        double riskFreeRate = (double) parameters.getOrDefault("riskFreeRate", 0.0);
        
        if (assetValue <= 0 || assetVolatility <= 0 || debtValue <= 0) {
            throw new IllegalArgumentException("Asset value, volatility, and debt value must be positive");
        }
        
        // Calculate distance to default
        double d2 = calculateDistanceToDefault(
            assetValue, debtValue, riskFreeRate, assetVolatility, timeHorizon
        );
        
        // Convert to probability of default
        // Using cumulative normal distribution
        return standardNormalCDF(-d2);
    }
    
    /**
     * Calculate distance to default (d2 in the Merton model)
     */
    private double calculateDistanceToDefault(double assetValue, double debtValue, 
                                            double riskFreeRate, double assetVolatility, double timeHorizon) {
        double d1 = (Math.log(assetValue / debtValue) + 
                   (riskFreeRate + 0.5 * assetVolatility * assetVolatility) * timeHorizon) / 
                   (assetVolatility * Math.sqrt(timeHorizon));
                   
        return d1 - assetVolatility * Math.sqrt(timeHorizon);
    }
    
    /**
     * Calculate credit spread based on the Merton model
     * 
     * @param timeHorizon Time horizon in years
     * @param parameters Model parameters
     * @return Credit spread in basis points
     */
    public double calculateCreditSpread(double timeHorizon, Map<String, Object> parameters) {
        // Extract parameters
        double assetValue = (double) parameters.getOrDefault("assetValue", 0.0);
        double assetVolatility = (double) parameters.getOrDefault("assetVolatility", 0.0);
        double debtValue = (double) parameters.getOrDefault("debtValue", 0.0);
        double riskFreeRate = (double) parameters.getOrDefault("riskFreeRate", 0.0);
        
        // Calculate d1 and d2
        double d1 = (Math.log(assetValue / debtValue) + 
                   (riskFreeRate + 0.5 * assetVolatility * assetVolatility) * timeHorizon) / 
                   (assetVolatility * Math.sqrt(timeHorizon));
        double d2 = d1 - assetVolatility * Math.sqrt(timeHorizon);
        
        // Calculate debt value based on Merton model
        double putOption = debtValue * Math.exp(-riskFreeRate * timeHorizon) * standardNormalCDF(-d2) - 
                          assetValue * standardNormalCDF(-d1);
        double riskDebtValue = debtValue * Math.exp(-riskFreeRate * timeHorizon) - putOption;
        
        // Calculate implied yield
        double impliedYield = -Math.log(riskDebtValue / debtValue) / timeHorizon;
        
        // Calculate credit spread
        double creditSpread = impliedYield - riskFreeRate;
        
        // Return in basis points
        return creditSpread * 10000.0;
    }
}

/**
 * Hazard rate model for credit risk (reduced-form model)
 */
class HazardRateModel extends CreditRiskModels {
    
    @Override
    public double calculateProbabilityOfDefault(double timeHorizon, Map<String, Object> parameters) {
        // Extract parameters
        double hazardRate = (double) parameters.getOrDefault("hazardRate", 0.0);
        boolean isTimeVarying = (boolean) parameters.getOrDefault("isTimeVarying", false);
        
        if (hazardRate < 0) {
            throw new IllegalArgumentException("Hazard rate must be non-negative");
        }
        
        if (isTimeVarying) {
            // For time-varying hazard rate, we need a function or array of values
            // Simplified implementation using piecewise constant hazard rates
            @SuppressWarnings("unchecked")
            Map<Double, Double> hazardRateCurve = (Map<Double, Double>) parameters.get("hazardRateCurve");
            
            if (hazardRateCurve == null || hazardRateCurve.isEmpty()) {
                throw new IllegalArgumentException("Hazard rate curve must be provided for time-varying model");
            }
            
            return calculateTimeVaryingPD(timeHorizon, hazardRateCurve);
        } else {
            // For constant hazard rate, use exponential formula
            return 1.0 - Math.exp(-hazardRate * timeHorizon);
        }
    }
    
    /**
     * Calculate probability of default with time-varying hazard rates
     */
    private double calculateTimeVaryingPD(double timeHorizon, Map<Double, Double> hazardRateCurve) {
        // Sort time points
        double[] timePoints = hazardRateCurve.keySet().stream()
                             .mapToDouble(Double::doubleValue)
                             .sorted()
                             .toArray();
        
        // Survival probability starts at 1
        double survivalProbability = 1.0;
        double currentTime = 0.0;
        
        for (int i = 0; i < timePoints.length; i++) {
            double nextTime = timePoints[i];
            
            if (nextTime > timeHorizon) {
                // We've reached the time horizon, calculate final segment
                double hazardRate = hazardRateCurve.get(timePoints[Math.max(0, i - 1)]);
                survivalProbability *= Math.exp(-hazardRate * (timeHorizon - currentTime));
                break;
            }
            
            if (i > 0) {
                // Calculate survival probability for this time segment
                double hazardRate = hazardRateCurve.get(timePoints[i - 1]);
                survivalProbability *= Math.exp(-hazardRate * (nextTime - currentTime));
            }
            
            currentTime = nextTime;
        }
        
        // If we've gone through all time points and still haven't reached the horizon
        if (currentTime < timeHorizon) {
            // Use the last hazard rate for the remaining time
            double hazardRate = hazardRateCurve.get(timePoints[timePoints.length - 1]);
            survivalProbability *= Math.exp(-hazardRate * (timeHorizon - currentTime));
        }
        
        return 1.0 - survivalProbability;
    }
    
    /**
     * Convert CDS spread to hazard rate
     * 
     * @param cdsSpreadBps CDS spread in basis points
     * @param recoveryRate Recovery rate (as a decimal)
     * @return Hazard rate
     */
    public double convertCdsSpreadToHazardRate(double cdsSpreadBps, double recoveryRate) {
        // Simple approximation: hazard rate ≈ spread / (1 - recovery rate)
        double spreadDecimal = cdsSpreadBps / 10000.0;
        return spreadDecimal / (1.0 - recoveryRate);
    }
    
    /**
     * Calculate survival probability curve
     * 
     * @param timeHorizons Array of time horizons
     * @param hazardRate Hazard rate
     * @return Map of time horizons to survival probabilities
     */
    public Map<Double, Double> calculateSurvivalCurve(double[] timeHorizons, double hazardRate) {
        Map<Double, Double> survivalCurve = new HashMap<>();
        
        for (double t : timeHorizons) {
            double survivalProbability = Math.exp(-hazardRate * t);
            survivalCurve.put(t, survivalProbability);
        }
        
        return survivalCurve;
    }
}

/**
 * Credit Metrics model for credit risk and portfolio management
 */
class CreditMetricsModel extends CreditRiskModels {
    
    // Credit rating transition matrix (example - should be calibrated to real data)
    // Rows and columns represent ratings (AAA, AA, A, BBB, BB, B, CCC, Default)
    private double[][] transitionMatrix = {
        {0.9078, 0.0848, 0.0060, 0.0010, 0.0004, 0.0000, 0.0000, 0.0000}, // AAA
        {0.0064, 0.9125, 0.0729, 0.0055, 0.0013, 0.0006, 0.0001, 0.0007}, // AA
        {0.0008, 0.0232, 0.9121, 0.0546, 0.0058, 0.0025, 0.0001, 0.0009}, // A
        {0.0002, 0.0030, 0.0541, 0.8750, 0.0565, 0.0077, 0.0011, 0.0024}, // BBB
        {0.0001, 0.0006, 0.0064, 0.0721, 0.8290, 0.0779, 0.0049, 0.0090}, // BB
        {0.0000, 0.0004, 0.0019, 0.0053, 0.0741, 0.8386, 0.0246, 0.0551}, // B
        {0.0000, 0.0000, 0.0022, 0.0033, 0.0114, 0.0452, 0.7430, 0.1949}  // CCC
    };
    
    @Override
    public double calculateProbabilityOfDefault(double timeHorizon, Map<String, Object> parameters) {
        // Extract parameters
        String initialRating = (String) parameters.getOrDefault("initialRating", "BBB");
        
        // Get index for initial rating
        int ratingIndex = getRatingIndex(initialRating);
        
        // For one-year horizon, return direct probability from transition matrix
        if (Math.abs(timeHorizon - 1.0) < 0.01) {
            return transitionMatrix[ratingIndex][7]; // Default column
        }
        
        // For multi-year horizon, compute power of transition matrix
        double[][] transitionPower = matrixPower(transitionMatrix, (int) Math.round(timeHorizon));
        return transitionPower[ratingIndex][7];
    }
    
    /**
     * Compute power of a matrix (simple implementation)
     */
    private double[][] matrixPower(double[][] matrix, int power) {
        if (power <= 0) {
            throw new IllegalArgumentException("Power must be positive");
        }
        
        if (power == 1) {
            return matrix;
        }
        
        int n = matrix.length;
        double[][] result = new double[n][n];
        
        // Initialize with identity matrix
        for (int i = 0; i < n; i++) {
            result[i][i] = 1.0;
        }
        
        // Compute matrix power using binary exponentiation
        double[][] temp = deepCopy(matrix);
        
        while (power > 0) {
            if (power % 2 == 1) {
                result = matrixMultiply(result, temp);
            }
            temp = matrixMultiply(temp, temp);
            power /= 2;
        }
        
        return result;
    }
    
    /**
     * Multiply two matrices
     */
    private double[][] matrixMultiply(double[][] a, double[][] b) {
        int n = a.length;
        double[][] result = new double[n][n];
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++) {
                    sum += a[i][k] * b[k][j];
                }
                result[i][j] = sum;
            }
        }
        
        return result;
    }
    
    /**
     * Deep copy a matrix
     */
    private double[][] deepCopy(double[][] matrix) {
        int n = matrix.length;
        double[][] copy = new double[n][n];
        
        for (int i = 0; i < n; i++) {
            System.arraycopy(matrix[i], 0, copy[i], 0, n);
        }
        
        return copy;
    }
    
    /**
     * Get the index for a rating
     */
    private int getRatingIndex(String rating) {
        switch (rating.toUpperCase()) {
            case "AAA": return 0;
            case "AA": return 1;
            case "A": return 2;
            case "BBB": return 3;
            case "BB": return 4;
            case "B": return 5;
            case "CCC": return 6;
            default:
                throw new IllegalArgumentException("Unknown rating: " + rating);
        }
    }
    
    /**
     * Run Monte Carlo simulation for portfolio credit risk
     * 
     * @param exposures Map of obligor IDs to exposure amounts
     * @param ratings Map of obligor IDs to credit ratings
     * @param recoveryRates Map of obligor IDs to recovery rates
     * @param correlationMatrix Correlation matrix between obligors
     * @param numSimulations Number of Monte Carlo simulations
     * @return Map of statistics about portfolio loss distribution
     */
    public Map<String, Double> simulatePortfolioLosses(
            Map<String, Double> exposures,
            Map<String, String> ratings,
            Map<String, Double> recoveryRates,
            Map<String, Map<String, Double>> correlationMatrix,
            int numSimulations) {
        
        String[] obligorIds = exposures.keySet().toArray(new String[0]);
        int numObligors = obligorIds.length;
        
        // Extract default probabilities for each obligor
        double[] defaultProbabilities = new double[numObligors];
        for (int i = 0; i < numObligors; i++) {
            String rating = ratings.get(obligorIds[i]);
            int ratingIndex = getRatingIndex(rating);
            defaultProbabilities[i] = transitionMatrix[ratingIndex][7];
        }
        
        // Generate correlated uniform random variables
        double[][] correlatedUniform = generateCorrelatedUniform(
            defaultProbabilities, correlationMatrix, obligorIds, numSimulations
        );
        
        // Run simulations
        double[] portfolioLosses = new double[numSimulations];
        
        for (int sim = 0; sim < numSimulations; sim++) {
            double totalLoss = 0.0;
            
            for (int i = 0; i < numObligors; i++) {
                boolean defaulted = correlatedUniform[sim][i] <= defaultProbabilities[i];
                
                if (defaulted) {
                    double exposure = exposures.get(obligorIds[i]);
                    double recoveryRate = recoveryRates.get(obligorIds[i]);
                    double loss = exposure * (1.0 - recoveryRate);
                    totalLoss += loss;
                }
            }
            
            portfolioLosses[sim] = totalLoss;
        }
        
        // Calculate statistics
        Map<String, Double> results = new HashMap<>();
        
        // Sort losses for percentile calculations
        double[] sortedLosses = portfolioLosses.clone();
        java.util.Arrays.sort(sortedLosses);
        
        // Expected loss
        double expectedLoss = 0.0;
        for (double loss : portfolioLosses) {
            expectedLoss += loss;
        }
        expectedLoss /= numSimulations;
        results.put("expectedLoss", expectedLoss);
        
        // Value at Risk (VaR)
        int var95Index = (int) Math.ceil(0.95 * numSimulations) - 1;
        int var99Index = (int) Math.ceil(0.99 * numSimulations) - 1;
        results.put("var95", sortedLosses[var95Index]);
        results.put("var99", sortedLosses[var99Index]);
        
        // Conditional VaR (CVaR) / Expected Shortfall
        double cvar95 = 0.0;
        for (int i = var95Index; i < numSimulations; i++) {
            cvar95 += sortedLosses[i];
        }
        cvar95 /= (numSimulations - var95Index);
        results.put("cvar95", cvar95);
        
        // Standard deviation
        double variance = 0.0;
        for (double loss : portfolioLosses) {
            variance += (loss - expectedLoss) * (loss - expectedLoss);
        }
        variance /= numSimulations;
        results.put("standardDeviation", Math.sqrt(variance));
        
        return results;
    }
    
    /**
     * Generate correlated uniform random variables
     */
    private double[][] generateCorrelatedUniform(
            double[] defaultProbabilities,
            Map<String, Map<String, Double>> correlationMatrix,
            String[] obligorIds,
            int numSimulations) {
        
        int numObligors = obligorIds.length;
        Random random = new Random();
        
        // Generate normal random variables
        double[][] normalRVs = new double[numSimulations][numObligors];
        for (int sim = 0; sim < numSimulations; sim++) {
            for (int i = 0; i < numObligors; i++) {
                normalRVs[sim][i] = random.nextGaussian();
            }
        }
        
        // Correlate the normal random variables
        double[][] correlatedNormal = new double[numSimulations][numObligors];
        
        for (int sim = 0; sim < numSimulations; sim++) {
            for (int i = 0; i < numObligors; i++) {
                double sum = 0.0;
                
                // Add contributions from other obligors based on correlation
                for (int j = 0; j < numObligors; j++) {
                    if (i == j) {
                        // Self correlation is 1.0
                        sum += normalRVs[sim][j];
                    } else {
                        // Use correlation matrix
                        double correlation = correlationMatrix.get(obligorIds[i]).get(obligorIds[j]);
                        sum += correlation * normalRVs[sim][j];
                    }
                }
                
                correlatedNormal[sim][i] = sum / Math.sqrt(numObligors);
            }
        }
        
        // Transform to uniform(0,1) using normal CDF
        double[][] correlatedUniform = new double[numSimulations][numObligors];
        
        for (int sim = 0; sim < numSimulations; sim++) {
            for (int i = 0; i < numObligors; i++) {
                // Convert to standard normal CDF using our own implementation
                correlatedUniform[sim][i] = standardNormalCDF(correlatedNormal[sim][i]);
            }
        }
        
        return correlatedUniform;
    }
} 