package com.quant.finance.models;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.time.LocalDate;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 * Class representing an interest rate yield curve
 */
public class YieldCurve {
    private final LocalDate valuationDate;
    private final TreeMap<Double, Double> tenorRates; // key: tenor in years, value: rate
    private final PolynomialSplineFunction interpolator;
    
    /**
     * Constructor for a yield curve with a set of tenor-rate pairs
     * 
     * @param valuationDate The date for which this yield curve is valid
     * @param tenorRates Map of tenors (in years) to interest rates (as decimals)
     */
    public YieldCurve(LocalDate valuationDate, Map<Double, Double> tenorRates) {
        this.valuationDate = valuationDate;
        this.tenorRates = new TreeMap<>(tenorRates);
        
        // Convert tenors and rates to arrays for interpolation
        double[] tenors = new double[this.tenorRates.size()];
        double[] rates = new double[this.tenorRates.size()];
        
        int i = 0;
        for (Map.Entry<Double, Double> entry : this.tenorRates.entrySet()) {
            tenors[i] = entry.getKey();
            rates[i] = entry.getValue();
            i++;
        }
        
        // Create interpolator
        LinearInterpolator linearInterpolator = new LinearInterpolator();
        this.interpolator = linearInterpolator.interpolate(tenors, rates);
    }
    
    /**
     * Get the zero rate for a specific tenor
     * 
     * @param tenor The tenor in years
     * @return The zero rate as a decimal
     */
    public double getRate(double tenor) {
        if (tenor <= 0) {
            throw new IllegalArgumentException("Tenor must be positive");
        }
        
        // Check if the tenor is within the interpolation range
        double minTenor = tenorRates.firstKey();
        double maxTenor = tenorRates.lastKey();
        
        if (tenor < minTenor) {
            // Use the first rate for tenors below the minimum
            return tenorRates.get(minTenor);
        } else if (tenor > maxTenor) {
            // Use the last rate for tenors above the maximum
            return tenorRates.get(maxTenor);
        } else {
            // Interpolate for tenors within the range
            return interpolator.value(tenor);
        }
    }
    
    /**
     * Get the discount factor for a specific tenor
     * 
     * @param tenor The tenor in years
     * @return The discount factor
     */
    public double getDiscountFactor(double tenor) {
        double rate = getRate(tenor);
        return Math.exp(-rate * tenor);
    }
    
    /**
     * Calculate the forward rate between two tenors
     * 
     * @param startTenor The starting tenor in years
     * @param endTenor The ending tenor in years
     * @return The forward rate as a decimal
     */
    public double getForwardRate(double startTenor, double endTenor) {
        if (startTenor >= endTenor) {
            throw new IllegalArgumentException("End tenor must be greater than start tenor");
        }
        
        double startRate = getRate(startTenor);
        double endRate = getRate(endTenor);
        
        // Calculate the forward rate
        double forwardRate = (endRate * endTenor - startRate * startTenor) / (endTenor - startTenor);
        
        return forwardRate;
    }
    
    /**
     * Create a par yield curve from a zero rate curve
     * 
     * @return A new YieldCurve object representing par yields
     */
    public YieldCurve toParYieldCurve() {
        Map<Double, Double> parRates = new HashMap<>();
        
        for (Double tenor : tenorRates.keySet()) {
            double sumDiscountFactors = 0.0;
            
            // Calculate the sum of discount factors up to the tenor
            for (int i = 1; i <= tenor.intValue(); i++) {
                sumDiscountFactors += getDiscountFactor(i);
            }
            
            // Calculate the par rate
            double finalDiscountFactor = getDiscountFactor(tenor);
            double parRate = (1 - finalDiscountFactor) / sumDiscountFactors;
            
            parRates.put(tenor, parRate);
        }
        
        return new YieldCurve(valuationDate, parRates);
    }
    
    /**
     * Build a yield curve from Treasury securities (simplified)
     * 
     * @param valuationDate The valuation date
     * @param treasuryYields Map of tenors to treasury yields
     * @return A new YieldCurve object
     */
    public static YieldCurve buildFromTreasuryYields(LocalDate valuationDate, Map<Double, Double> treasuryYields) {
        // In a real system, this would involve more sophisticated bootstrapping
        // and potentially curve fitting. This is a simplified version.
        return new YieldCurve(valuationDate, treasuryYields);
    }
    
    /**
     * Build a yield curve using the Nelson-Siegel model (simplified)
     * 
     * @param valuationDate The valuation date
     * @param beta0 Long-term level parameter
     * @param beta1 Short-term slope parameter
     * @param beta2 Medium-term curvature parameter
     * @param tau Time decay parameter
     * @param tenors Array of tenors for which to calculate rates
     * @return A new YieldCurve object
     */
    public static YieldCurve buildNelsonSiegelCurve(LocalDate valuationDate, double beta0, double beta1, 
                                                   double beta2, double tau, double[] tenors) {
        Map<Double, Double> rates = new HashMap<>();
        
        for (double tenor : tenors) {
            if (tenor <= 0) continue;
            
            double expTerm = Math.exp(-tenor / tau);
            double term1 = beta0;
            double term2 = beta1 * ((1 - expTerm) / (tenor / tau));
            double term3 = beta2 * ((1 - expTerm) / (tenor / tau) - expTerm);
            
            double rate = term1 + term2 + term3;
            rates.put(tenor, rate);
        }
        
        return new YieldCurve(valuationDate, rates);
    }
    
    /**
     * Calculate the price of a zero-coupon bond
     * 
     * @param maturityTenor The time to maturity in years
     * @param faceValue The face value of the bond
     * @return The price of the zero-coupon bond
     */
    public double priceZeroCouponBond(double maturityTenor, double faceValue) {
        double discountFactor = getDiscountFactor(maturityTenor);
        return faceValue * discountFactor;
    }
    
    // Getters
    public LocalDate getValuationDate() {
        return valuationDate;
    }
    
    public Map<Double, Double> getTenorRates() {
        return new HashMap<>(tenorRates);
    }
} 