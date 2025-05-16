package com.quant.finance.utils;

import java.time.LocalDate;
import java.util.Map;

/**
 * Utility class providing basic financial mathematics calculations
 * that can be used across different financial instruments and models.
 */
public class FinancialMathUtils {
    
    /**
     * Calculate discount factor from a continuously compounded rate
     * 
     * @param rate Interest rate (as a decimal)
     * @param timeInYears Time in years
     * @return Discount factor
     */
    public static double calculateDiscountFactor(double rate, double timeInYears) {
        return Math.exp(-rate * timeInYears);
    }
    
    /**
     * Calculate discount factor from a continuously compounded rate
     * 
     * @param rate Interest rate (as a decimal)
     * @param startDate Start date
     * @param endDate End date
     * @return Discount factor
     */
    public static double calculateDiscountFactor(double rate, LocalDate startDate, LocalDate endDate) {
        double timeInYears = DateUtils.yearFraction(startDate, endDate);
        return calculateDiscountFactor(rate, timeInYears);
    }
    
    /**
     * Calculate the present value of a future cash flow
     * 
     * @param cashFlow Future cash flow amount
     * @param rate Discount rate (as a decimal)
     * @param timeInYears Time to cash flow in years
     * @return Present value
     */
    public static double presentValue(double cashFlow, double rate, double timeInYears) {
        return cashFlow * calculateDiscountFactor(rate, timeInYears);
    }
    
    /**
     * Calculate the future value of a present amount
     * 
     * @param presentAmount Present amount
     * @param rate Interest rate (as a decimal)
     * @param timeInYears Time in years
     * @return Future value
     */
    public static double futureValue(double presentAmount, double rate, double timeInYears) {
        return presentAmount * Math.exp(rate * timeInYears);
    }
    
    /**
     * Calculate the internal rate of return for a series of cash flows
     * 
     * @param cashFlows Array of cash flows (negative for outflows, positive for inflows)
     * @param times Array of times (in years) corresponding to each cash flow
     * @return Internal rate of return (as a decimal)
     */
    public static double calculateIRR(double[] cashFlows, double[] times) {
        // Newton-Raphson method for IRR calculation
        double guess = 0.1; // Initial guess
        double tolerance = 1e-10;
        int maxIterations = 100;
        
        for (int i = 0; i < maxIterations; i++) {
            double f = 0.0;
            double df = 0.0;
            
            for (int j = 0; j < cashFlows.length; j++) {
                double discountFactor = Math.exp(-guess * times[j]);
                f += cashFlows[j] * discountFactor;
                df -= cashFlows[j] * times[j] * discountFactor;
            }
            
            double delta = f / df;
            guess -= delta;
            
            if (Math.abs(delta) < tolerance) {
                return guess;
            }
        }
        
        throw new RuntimeException("IRR calculation did not converge");
    }
    
    /**
     * Calculate the net present value (NPV) of a series of cash flows
     * 
     * @param cashFlows Array of cash flows
     * @param times Array of times (in years) corresponding to each cash flow
     * @param rate Discount rate (as a decimal)
     * @return Net present value
     */
    public static double calculateNPV(double[] cashFlows, double[] times, double rate) {
        double npv = 0.0;
        
        for (int i = 0; i < cashFlows.length; i++) {
            npv += presentValue(cashFlows[i], rate, times[i]);
        }
        
        return npv;
    }
    
    /**
     * Calculate yield to maturity for a bond given its price and cash flows
     * 
     * @param price Bond price
     * @param cashFlows Array of cash flows (coupon payments and principal)
     * @param times Array of times (in years) to each cash flow
     * @return Yield to maturity (as a decimal)
     */
    public static double calculateYieldToMaturity(double price, double[] cashFlows, double[] times) {
        // Create a new array with initial price as negative cash flow at time 0
        double[] adjustedCashFlows = new double[cashFlows.length + 1];
        double[] adjustedTimes = new double[times.length + 1];
        
        adjustedCashFlows[0] = -price;
        adjustedTimes[0] = 0.0;
        
        System.arraycopy(cashFlows, 0, adjustedCashFlows, 1, cashFlows.length);
        System.arraycopy(times, 0, adjustedTimes, 1, times.length);
        
        return calculateIRR(adjustedCashFlows, adjustedTimes);
    }
    
    /**
     * Calculate the modified duration of a bond
     * 
     * @param cashFlows Array of cash flows
     * @param times Array of times (in years) to each cash flow
     * @param yield Yield (as a decimal)
     * @return Modified duration
     */
    public static double calculateModifiedDuration(double[] cashFlows, double[] times, double yield) {
        double price = calculateNPV(cashFlows, times, yield);
        double macaulayDuration = calculateMacaulayDuration(cashFlows, times, yield);
        
        return macaulayDuration / (1 + yield);
    }
    
    /**
     * Calculate the Macaulay duration of a bond
     * 
     * @param cashFlows Array of cash flows
     * @param times Array of times (in years) to each cash flow
     * @param yield Yield (as a decimal)
     * @return Macaulay duration
     */
    public static double calculateMacaulayDuration(double[] cashFlows, double[] times, double yield) {
        double price = calculateNPV(cashFlows, times, yield);
        double weightedSum = 0.0;
        
        for (int i = 0; i < cashFlows.length; i++) {
            double pv = presentValue(cashFlows[i], yield, times[i]);
            weightedSum += times[i] * pv;
        }
        
        return weightedSum / price;
    }
    
    /**
     * Calculate the convexity of a bond
     * 
     * @param cashFlows Array of cash flows
     * @param times Array of times (in years) to each cash flow
     * @param yield Yield (as a decimal)
     * @return Convexity
     */
    public static double calculateConvexity(double[] cashFlows, double[] times, double yield) {
        double price = calculateNPV(cashFlows, times, yield);
        double weightedSum = 0.0;
        
        for (int i = 0; i < cashFlows.length; i++) {
            double pv = presentValue(cashFlows[i], yield, times[i]);
            weightedSum += times[i] * (times[i] + 1) * pv;
        }
        
        return weightedSum / (price * Math.pow(1 + yield, 2));
    }
    
    /**
     * Calculate the price change for a given yield change, using duration and convexity
     * 
     * @param price Current price
     * @param duration Modified duration
     * @param convexity Convexity
     * @param yieldChange Change in yield (in decimal)
     * @return Price change as a percentage
     */
    public static double calculatePriceChangeWithConvexity(double price, double duration, 
                                                          double convexity, double yieldChange) {
        double firstOrderTerm = -duration * yieldChange;
        double secondOrderTerm = 0.5 * convexity * yieldChange * yieldChange;
        
        return firstOrderTerm + secondOrderTerm;
    }
    
    /**
     * Linear interpolation between two points
     * 
     * @param x Point to interpolate at
     * @param x0 First x coordinate
     * @param y0 First y coordinate
     * @param x1 Second x coordinate
     * @param y1 Second y coordinate
     * @return Interpolated value
     */
    public static double linearInterpolate(double x, double x0, double y0, double x1, double y1) {
        if (Math.abs(x1 - x0) < 1e-10) {
            return (y0 + y1) / 2.0;
        }
        
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    }
    
    /**
     * Cubic spline interpolation
     * This is a simplified implementation for demonstration purposes
     * 
     * @param x Value to interpolate at
     * @param xValues Array of x coordinates
     * @param yValues Array of y coordinates
     * @return Interpolated value
     */
    public static double cubicSplineInterpolate(double x, double[] xValues, double[] yValues) {
        // Find the interval
        int i = 0;
        while (i < xValues.length - 1 && xValues[i + 1] < x) {
            i++;
        }
        
        // If outside the range, use linear extrapolation
        if (i >= xValues.length - 1) {
            return linearInterpolate(x, xValues[xValues.length - 2], yValues[yValues.length - 2],
                                    xValues[xValues.length - 1], yValues[yValues.length - 1]);
        }
        
        if (i < 0) {
            return linearInterpolate(x, xValues[0], yValues[0], xValues[1], yValues[1]);
        }
        
        // For simplicity, we'll just use linear interpolation in this implementation
        // A full cubic spline would require setting up and solving a system of equations
        return linearInterpolate(x, xValues[i], yValues[i], xValues[i + 1], yValues[i + 1]);
    }
} 