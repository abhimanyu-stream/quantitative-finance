package com.quant.finance.models;

import com.quant.finance.utils.DateUtils;
import com.quant.finance.utils.FinancialMathUtils;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Abstract base class for yield curve models
 */
public abstract class YieldCurveModels {
    
    protected LocalDate valuationDate;
    protected TreeMap<Double, Double> tenorRates; // Map of tenor (in years) to zero rates
    
    /**
     * Constructor
     * 
     * @param valuationDate Valuation date
     * @param tenorRates Map of tenors (in years) to zero rates
     */
    public YieldCurveModels(LocalDate valuationDate, Map<Double, Double> tenorRates) {
        this.valuationDate = valuationDate;
        this.tenorRates = new TreeMap<>(tenorRates);
    }
    
    /**
     * Get the zero rate for a given tenor
     * 
     * @param tenor Tenor in years
     * @return Zero rate (as a decimal)
     */
    public abstract double getRate(double tenor);
    
    /**
     * Get the discount factor for a given tenor
     * 
     * @param tenor Tenor in years
     * @return Discount factor
     */
    public double getDiscountFactor(double tenor) {
        double rate = getRate(tenor);
        return FinancialMathUtils.calculateDiscountFactor(rate, tenor);
    }
    
    /**
     * Get the forward rate between two tenors
     * 
     * @param startTenor Start tenor in years
     * @param endTenor End tenor in years
     * @return Forward rate (as a decimal)
     */
    public double getForwardRate(double startTenor, double endTenor) {
        if (startTenor >= endTenor) {
            throw new IllegalArgumentException("Start tenor must be less than end tenor");
        }
        
        double startRate = getRate(startTenor);
        double endRate = getRate(endTenor);
        
        // Calculate forward rate from discount factors
        double startDF = Math.exp(-startRate * startTenor);
        double endDF = Math.exp(-endRate * endTenor);
        
        double forwardPeriod = endTenor - startTenor;
        return -Math.log(endDF / startDF) / forwardPeriod;
    }
    
    /**
     * Price a zero-coupon bond with a given maturity and face value
     * 
     * @param maturityInYears Maturity in years
     * @param faceValue Face value of the bond
     * @return Bond price
     */
    public double priceZeroCouponBond(double maturityInYears, double faceValue) {
        return faceValue * getDiscountFactor(maturityInYears);
    }
    
    /**
     * Create a yield curve model using the specified interpolation method
     * 
     * @param valuationDate Valuation date
     * @param tenorRates Map of tenors (in years) to zero rates
     * @param interpolationMethod Interpolation method to use
     * @return YieldCurveModels instance
     */
    public static YieldCurveModels createModel(LocalDate valuationDate, Map<Double, Double> tenorRates, 
                                              String interpolationMethod) {
        switch (interpolationMethod.toLowerCase()) {
            case "linear":
                return new LinearInterpolationYieldCurve(valuationDate, tenorRates);
            case "loglinear":
                return new LogLinearInterpolationYieldCurve(valuationDate, tenorRates);
            case "cubicspline":
                return new CubicSplineInterpolationYieldCurve(valuationDate, tenorRates);
            case "nelsonsiegel":
                return new NelsonSiegelYieldCurve(valuationDate, tenorRates);
            default:
                throw new IllegalArgumentException("Unknown interpolation method: " + interpolationMethod);
        }
    }
}

/**
 * Linear interpolation yield curve model
 */
class LinearInterpolationYieldCurve extends YieldCurveModels {
    
    public LinearInterpolationYieldCurve(LocalDate valuationDate, Map<Double, Double> tenorRates) {
        super(valuationDate, tenorRates);
    }
    
    @Override
    public double getRate(double tenor) {
        if (tenorRates.containsKey(tenor)) {
            return tenorRates.get(tenor);
        }
        
        // Find surrounding tenors
        Double lowerTenor = tenorRates.floorKey(tenor);
        Double upperTenor = tenorRates.ceilingKey(tenor);
        
        // Handle edge cases
        if (lowerTenor == null) {
            return tenorRates.get(upperTenor);
        }
        
        if (upperTenor == null) {
            return tenorRates.get(lowerTenor);
        }
        
        // Perform linear interpolation
        double lowerRate = tenorRates.get(lowerTenor);
        double upperRate = tenorRates.get(upperTenor);
        
        return FinancialMathUtils.linearInterpolate(tenor, lowerTenor, lowerRate, upperTenor, upperRate);
    }
}

/**
 * Log-linear interpolation yield curve model
 * Interpolates linearly in the log of discount factors
 */
class LogLinearInterpolationYieldCurve extends YieldCurveModels {
    
    public LogLinearInterpolationYieldCurve(LocalDate valuationDate, Map<Double, Double> tenorRates) {
        super(valuationDate, tenorRates);
    }
    
    @Override
    public double getRate(double tenor) {
        if (tenorRates.containsKey(tenor)) {
            return tenorRates.get(tenor);
        }
        
        // Find surrounding tenors
        Double lowerTenor = tenorRates.floorKey(tenor);
        Double upperTenor = tenorRates.ceilingKey(tenor);
        
        // Handle edge cases
        if (lowerTenor == null) {
            return tenorRates.get(upperTenor);
        }
        
        if (upperTenor == null) {
            return tenorRates.get(lowerTenor);
        }
        
        // Convert rates to discount factors
        double lowerRate = tenorRates.get(lowerTenor);
        double upperRate = tenorRates.get(upperTenor);
        
        double lowerDF = Math.exp(-lowerRate * lowerTenor);
        double upperDF = Math.exp(-upperRate * upperTenor);
        
        // Interpolate in log space
        double logLowerDF = Math.log(lowerDF);
        double logUpperDF = Math.log(upperDF);
        
        double interpolatedLogDF = FinancialMathUtils.linearInterpolate(
            tenor, lowerTenor, logLowerDF, upperTenor, logUpperDF
        );
        
        // Convert back to discount factor and then to rate
        double interpolatedDF = Math.exp(interpolatedLogDF);
        return -Math.log(interpolatedDF) / tenor;
    }
}

/**
 * Cubic spline interpolation yield curve model
 */
class CubicSplineInterpolationYieldCurve extends YieldCurveModels {
    private double[] tenors;
    private double[] rates;
    
    public CubicSplineInterpolationYieldCurve(LocalDate valuationDate, Map<Double, Double> tenorRates) {
        super(valuationDate, tenorRates);
        initializeArrays();
    }
    
    private void initializeArrays() {
        int size = tenorRates.size();
        tenors = new double[size];
        rates = new double[size];
        
        int i = 0;
        for (Map.Entry<Double, Double> entry : tenorRates.entrySet()) {
            tenors[i] = entry.getKey();
            rates[i] = entry.getValue();
            i++;
        }
    }
    
    @Override
    public double getRate(double tenor) {
        if (tenorRates.containsKey(tenor)) {
            return tenorRates.get(tenor);
        }
        
        // Use cubic spline interpolation
        return FinancialMathUtils.cubicSplineInterpolate(tenor, tenors, rates);
    }
}

/**
 * Nelson-Siegel yield curve model
 * Fits a parametric model to the yield curve
 */
class NelsonSiegelYieldCurve extends YieldCurveModels {
    private double beta0; // Long-term interest rate level
    private double beta1; // Short-term component
    private double beta2; // Medium-term component
    private double tau;   // Time decay parameter
    
    public NelsonSiegelYieldCurve(LocalDate valuationDate, Map<Double, Double> tenorRates) {
        super(valuationDate, tenorRates);
        calibrateParameters();
    }
    
    /**
     * Calibrate the Nelson-Siegel parameters to match the input rates
     */
    private void calibrateParameters() {
        // For simplicity, we'll use fixed parameters rather than calibration
        // In a real implementation, this would involve an optimization algorithm
        beta0 = 0.03;    // 3% long-term rate
        beta1 = -0.02;   // Negative slope (upward sloping curve)
        beta2 = 0.01;    // Hump/trough
        tau = 2.0;       // Time decay parameter
        
        // Adjust beta0 to match the long-term rate
        Double maxTenor = tenorRates.lastKey();
        if (maxTenor >= 10.0) {
            beta0 = tenorRates.get(maxTenor);
        }
    }
    
    @Override
    public double getRate(double tenor) {
        if (tenor <= 0) {
            throw new IllegalArgumentException("Tenor must be positive");
        }
        
        if (tenorRates.containsKey(tenor)) {
            return tenorRates.get(tenor);
        }
        
        // Nelson-Siegel formula
        double term1 = beta0;
        double term2 = beta1 * (1 - Math.exp(-tenor / tau)) / (tenor / tau);
        double term3 = beta2 * ((1 - Math.exp(-tenor / tau)) / (tenor / tau) - Math.exp(-tenor / tau));
        
        return term1 + term2 + term3;
    }
    
    /**
     * Get the parameters of the Nelson-Siegel model
     * 
     * @return Map of parameter names to values
     */
    public Map<String, Double> getParameters() {
        Map<String, Double> params = new HashMap<>();
        params.put("beta0", beta0);
        params.put("beta1", beta1);
        params.put("beta2", beta2);
        params.put("tau", tau);
        return params;
    }
} 