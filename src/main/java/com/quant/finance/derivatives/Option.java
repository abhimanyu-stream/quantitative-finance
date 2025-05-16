package com.quant.finance.derivatives;

import com.quant.finance.instruments.Instrument;
import com.quant.finance.utils.DateUtils;

import java.time.LocalDate;
import java.util.HashMap;
import java.util.Map;

/**
 * Class representing financial options (calls and puts)
 */
public class Option extends Instrument {
    
    public enum OptionType {
        CALL, PUT
    }
    
    public enum ExerciseStyle {
        EUROPEAN, AMERICAN, BERMUDAN
    }
    
    private final OptionType optionType;
    private final ExerciseStyle exerciseStyle;
    private final double strikePrice;
    private final LocalDate expiryDate;
    private final String underlyingId;
    
    public Option(String id, String name, LocalDate issueDate, OptionType optionType, 
                  ExerciseStyle exerciseStyle, double strikePrice, LocalDate expiryDate,
                  String underlyingId) {
        super(id, name, issueDate);
        this.optionType = optionType;
        this.exerciseStyle = exerciseStyle;
        this.strikePrice = strikePrice;
        this.expiryDate = expiryDate;
        this.underlyingId = underlyingId;
    }
    
    @Override
    public double presentValue(LocalDate valuationDate, Map<String, Object> marketData) {
        if (valuationDate.isAfter(expiryDate)) {
            return 0.0; // Expired option has zero value
        }
        
        // Extract required market data
        double spotPrice = (double) marketData.getOrDefault("spotPrice", 0.0);
        double riskFreeRate = (double) marketData.getOrDefault("riskFreeRate", 0.02);
        double volatility = (double) marketData.getOrDefault("volatility", 0.2);
        
        // Calculate time to expiry in years
        double timeToExpiry = DateUtils.yearFraction(valuationDate, expiryDate);
        
        // Use Black-Scholes for European options
        if (exerciseStyle == ExerciseStyle.EUROPEAN) {
            return blackScholesPrice(optionType, spotPrice, strikePrice, timeToExpiry, riskFreeRate, volatility);
        } else {
            // For American options, use binomial model or other approximation
            // This is a simplified version - a real implementation would use a more sophisticated model
            double europeanPrice = blackScholesPrice(optionType, spotPrice, strikePrice, timeToExpiry, riskFreeRate, volatility);
            
            // Apply a premium for American options
            // In reality, we would use a proper American option pricing model
            if (optionType == OptionType.PUT) {
                // American puts are worth more than European due to early exercise
                return europeanPrice * 1.05;
            } else {
                // For non-dividend paying stocks, American calls = European calls
                return europeanPrice;
            }
        }
    }
    
    @Override
    public Map<String, Double> calculateRisks(LocalDate valuationDate, Map<String, Object> marketData) {
        Map<String, Double> risks = new HashMap<>();
        
        // Extract required market data
        double spotPrice = (double) marketData.getOrDefault("spotPrice", 0.0);
        double riskFreeRate = (double) marketData.getOrDefault("riskFreeRate", 0.02);
        double volatility = (double) marketData.getOrDefault("volatility", 0.2);
        
        // Time to expiry in years
        double timeToExpiry = DateUtils.yearFraction(valuationDate, expiryDate);
        
        // Calculate option Greeks
        if (exerciseStyle == ExerciseStyle.EUROPEAN) {
            // Delta - sensitivity to underlying price
            risks.put("delta", calculateDelta(spotPrice, riskFreeRate, volatility, timeToExpiry));
            
            // Gamma - rate of change of delta
            risks.put("gamma", calculateGamma(spotPrice, riskFreeRate, volatility, timeToExpiry));
            
            // Theta - sensitivity to time decay
            risks.put("theta", calculateTheta(spotPrice, riskFreeRate, volatility, timeToExpiry));
            
            // Vega - sensitivity to volatility
            risks.put("vega", calculateVega(spotPrice, riskFreeRate, volatility, timeToExpiry));
            
            // Rho - sensitivity to interest rate
            risks.put("rho", calculateRho(spotPrice, riskFreeRate, volatility, timeToExpiry));
        } else {
            // For American options, greeks would be calculated using finite difference methods
            // or other numerical approaches. This is simplified.
            risks.put("delta", 0.0);
            risks.put("gamma", 0.0);
            risks.put("theta", 0.0);
            risks.put("vega", 0.0);
            risks.put("rho", 0.0);
        }
        
        return risks;
    }
    
    // Black-Scholes pricing formula
    private double blackScholesPrice(OptionType type, double spot, double strike, double time, 
                                    double rate, double volatility) {
        if (time <= 0) {
            return Math.max(0, type == OptionType.CALL ? spot - strike : strike - spot);
        }
        
        double d1 = (Math.log(spot / strike) + (rate + 0.5 * volatility * volatility) * time) 
                    / (volatility * Math.sqrt(time));
        double d2 = d1 - volatility * Math.sqrt(time);
        
        if (type == OptionType.CALL) {
            return spot * cumulativeNormalDistribution(d1) 
                   - strike * Math.exp(-rate * time) * cumulativeNormalDistribution(d2);
        } else {
            return strike * Math.exp(-rate * time) * cumulativeNormalDistribution(-d2) 
                   - spot * cumulativeNormalDistribution(-d1);
        }
    }
    
    // Standard normal cumulative distribution function approximation
    private double cumulativeNormalDistribution(double x) {
        // Approximation from Abramowitz and Stegun
        double b1 = 0.31938153;
        double b2 = -0.356563782;
        double b3 = 1.781477937;
        double b4 = -1.821255978;
        double b5 = 1.330274429;
        double p = 0.2316419;
        double c = 0.39894228;
        
        if (x >= 0.0) {
            double t = 1.0 / (1.0 + p * x);
            return 1.0 - c * Math.exp(-x * x / 2.0) * t * (t * (t * (t * (t * b5 + b4) + b3) + b2) + b1);
        } else {
            double t = 1.0 / (1.0 - p * x);
            return c * Math.exp(-x * x / 2.0) * t * (t * (t * (t * (t * b5 + b4) + b3) + b2) + b1);
        }
    }
    
    // Calculate Delta - first derivative of price with respect to underlying price
    private double calculateDelta(double spotPrice, double riskFreeRate, double volatility, double timeToExpiry) {
        if (timeToExpiry <= 0) {
            return 0.0;
        }
        
        double d1 = (Math.log(spotPrice / strikePrice) + (riskFreeRate + 0.5 * volatility * volatility) * timeToExpiry) 
                    / (volatility * Math.sqrt(timeToExpiry));
                    
        if (optionType == OptionType.CALL) {
            return cumulativeNormalDistribution(d1);
        } else {
            return cumulativeNormalDistribution(d1) - 1.0;
        }
    }
    
    // Calculate Gamma - second derivative of price with respect to underlying price
    private double calculateGamma(double spotPrice, double riskFreeRate, double volatility, double timeToExpiry) {
        if (timeToExpiry <= 0) {
            return 0.0;
        }
        
        double d1 = (Math.log(spotPrice / strikePrice) + (riskFreeRate + 0.5 * volatility * volatility) * timeToExpiry) 
                    / (volatility * Math.sqrt(timeToExpiry));
        
        // Standard normal probability density function
        double pdf = Math.exp(-0.5 * d1 * d1) / Math.sqrt(2.0 * Math.PI);
        
        return pdf / (spotPrice * volatility * Math.sqrt(timeToExpiry));
    }
    
    // Calculate Theta - derivative of price with respect to time to expiry
    private double calculateTheta(double spotPrice, double riskFreeRate, double volatility, double timeToExpiry) {
        if (timeToExpiry <= 0) {
            return 0.0;
        }
        
        double d1 = (Math.log(spotPrice / strikePrice) + (riskFreeRate + 0.5 * volatility * volatility) * timeToExpiry) 
                    / (volatility * Math.sqrt(timeToExpiry));
        double d2 = d1 - volatility * Math.sqrt(timeToExpiry);
        
        // Standard normal probability density function
        double pdf = Math.exp(-0.5 * d1 * d1) / Math.sqrt(2.0 * Math.PI);
        
        if (optionType == OptionType.CALL) {
            double term1 = -spotPrice * pdf * volatility / (2.0 * Math.sqrt(timeToExpiry));
            double term2 = -riskFreeRate * strikePrice * Math.exp(-riskFreeRate * timeToExpiry) * cumulativeNormalDistribution(d2);
            return term1 + term2;
        } else {
            double term1 = -spotPrice * pdf * volatility / (2.0 * Math.sqrt(timeToExpiry));
            double term2 = riskFreeRate * strikePrice * Math.exp(-riskFreeRate * timeToExpiry) * cumulativeNormalDistribution(-d2);
            return term1 + term2;
        }
    }
    
    // Calculate Vega - derivative of price with respect to volatility
    private double calculateVega(double spotPrice, double riskFreeRate, double volatility, double timeToExpiry) {
        if (timeToExpiry <= 0) {
            return 0.0;
        }
        
        double d1 = (Math.log(spotPrice / strikePrice) + (riskFreeRate + 0.5 * volatility * volatility) * timeToExpiry) 
                    / (volatility * Math.sqrt(timeToExpiry));
        
        // Standard normal probability density function
        double pdf = Math.exp(-0.5 * d1 * d1) / Math.sqrt(2.0 * Math.PI);
        
        return spotPrice * Math.sqrt(timeToExpiry) * pdf * 0.01; // Scale to 1% change in volatility
    }
    
    // Calculate Rho - derivative of price with respect to interest rate
    private double calculateRho(double spotPrice, double riskFreeRate, double volatility, double timeToExpiry) {
        if (timeToExpiry <= 0) {
            return 0.0;
        }
        
        double d1 = (Math.log(spotPrice / strikePrice) + (riskFreeRate + 0.5 * volatility * volatility) * timeToExpiry) 
                    / (volatility * Math.sqrt(timeToExpiry));
        double d2 = d1 - volatility * Math.sqrt(timeToExpiry);
        
        if (optionType == OptionType.CALL) {
            return strikePrice * timeToExpiry * Math.exp(-riskFreeRate * timeToExpiry) 
                   * cumulativeNormalDistribution(d2) * 0.01; // Scale to 1% change in rate
        } else {
            return -strikePrice * timeToExpiry * Math.exp(-riskFreeRate * timeToExpiry) 
                    * cumulativeNormalDistribution(-d2) * 0.01; // Scale to 1% change in rate
        }
    }
    
    // Getters
    public OptionType getOptionType() {
        return optionType;
    }
    
    public ExerciseStyle getExerciseStyle() {
        return exerciseStyle;
    }
    
    public double getStrikePrice() {
        return strikePrice;
    }
    
    public LocalDate getExpiryDate() {
        return expiryDate;
    }
    
    public String getUnderlyingId() {
        return underlyingId;
    }
} 