package com.quant.finance.instruments;

import com.quant.finance.models.YieldCurveModels;
import com.quant.finance.utils.DateUtils;
import com.quant.finance.utils.FinancialMathUtils;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class for fixed income analytics calculations
 */
public class FixedIncomeAnalytics {
    
    /**
     * Calculate bond cash flows
     * 
     * @param bond Bond object
     * @return List of cash flow dates and amounts
     */
    public static List<Map.Entry<LocalDate, Double>> calculateBondCashFlows(Bond bond) {
        List<Map.Entry<LocalDate, Double>> cashFlows = new ArrayList<>();
        
        // Determine payment frequency in months
        int paymentsPerYear = bond.getPaymentsPerYear();
        int paymentFrequencyMonths = 12 / paymentsPerYear;
        
        LocalDate currentDate = bond.getIssueDate();
        LocalDate maturityDate = bond.getMaturityDate();
        
        // Get notional and coupon rate from bond properties
        double notionalAmount = 1000.0; // Default value
        double couponRate = 0.0;        // Default value
        
        if (bond.getProperty("notionalAmount") != null) {
            notionalAmount = (double) bond.getProperty("notionalAmount");
        }
        
        if (bond.getProperty("couponRate") != null) {
            couponRate = (double) bond.getProperty("couponRate");
        }
        
        // Add coupon payments
        while (currentDate.isBefore(maturityDate)) {
            currentDate = currentDate.plusMonths(paymentFrequencyMonths);
            
            if (!currentDate.isAfter(maturityDate)) {
                double couponPayment = notionalAmount * couponRate / paymentsPerYear;
                cashFlows.add(Map.entry(currentDate, couponPayment));
            }
        }
        
        // Add principal repayment at maturity
        double lastPayment = cashFlows.get(cashFlows.size() - 1).getValue() + notionalAmount;
        cashFlows.set(cashFlows.size() - 1, Map.entry(maturityDate, lastPayment));
        
        return cashFlows;
    }
    
    /**
     * Calculate bond price using a yield curve
     * 
     * @param bond Bond object
     * @param valuationDate Valuation date
     * @param yieldCurve Yield curve for discounting
     * @return Bond price
     */
    public static double calculateBondPriceWithYieldCurve(Bond bond, LocalDate valuationDate, YieldCurveModels yieldCurve) {
        List<Map.Entry<LocalDate, Double>> cashFlows = calculateBondCashFlows(bond);
        double price = 0.0;
        
        for (Map.Entry<LocalDate, Double> cashFlow : cashFlows) {
            LocalDate cashFlowDate = cashFlow.getKey();
            double amount = cashFlow.getValue();
            
            if (cashFlowDate.isAfter(valuationDate)) {
                double timeToPayment = DateUtils.yearFraction(valuationDate, cashFlowDate);
                double discountFactor = yieldCurve.getDiscountFactor(timeToPayment);
                price += amount * discountFactor;
            }
        }
        
        return price;
    }
    
    /**
     * Calculate bond duration measures
     * 
     * @param bond Bond object
     * @param valuationDate Valuation date
     * @param yield Yield to maturity
     * @return Map containing Macaulay, modified, and effective duration
     */
    public static Map<String, Double> calculateBondDurations(Bond bond, LocalDate valuationDate, double yield) {
        // Get cash flows
        List<Map.Entry<LocalDate, Double>> cashFlows = calculateBondCashFlows(bond);
        
        // Convert to arrays for financial math calculations
        double[] amounts = new double[cashFlows.size()];
        double[] times = new double[cashFlows.size()];
        
        for (int i = 0; i < cashFlows.size(); i++) {
            Map.Entry<LocalDate, Double> cashFlow = cashFlows.get(i);
            LocalDate cashFlowDate = cashFlow.getKey();
            double amount = cashFlow.getValue();
            
            if (cashFlowDate.isAfter(valuationDate)) {
                times[i] = DateUtils.yearFraction(valuationDate, cashFlowDate);
                amounts[i] = amount;
            }
        }
        
        // Calculate price
        double price = FinancialMathUtils.calculateNPV(amounts, times, yield);
        
        // Calculate durations
        double macaulayDuration = FinancialMathUtils.calculateMacaulayDuration(amounts, times, yield);
        double modifiedDuration = FinancialMathUtils.calculateModifiedDuration(amounts, times, yield);
        
        // Calculate effective duration using price sensitivity
        double yieldShift = 0.0001; // 1 basis point
        double priceUp = FinancialMathUtils.calculateNPV(amounts, times, yield - yieldShift);
        double priceDown = FinancialMathUtils.calculateNPV(amounts, times, yield + yieldShift);
        double effectiveDuration = (priceUp - priceDown) / (2 * price * yieldShift);
        
        // Return results
        Map<String, Double> durations = new HashMap<>();
        durations.put("macaulayDuration", macaulayDuration);
        durations.put("modifiedDuration", modifiedDuration);
        durations.put("effectiveDuration", effectiveDuration);
        
        return durations;
    }
    
    /**
     * Calculate bond convexity
     * 
     * @param bond Bond object
     * @param valuationDate Valuation date
     * @param yield Yield to maturity
     * @return Convexity value
     */
    public static double calculateBondConvexity(Bond bond, LocalDate valuationDate, double yield) {
        // Get cash flows
        List<Map.Entry<LocalDate, Double>> cashFlows = calculateBondCashFlows(bond);
        
        // Convert to arrays for financial math calculations
        double[] amounts = new double[cashFlows.size()];
        double[] times = new double[cashFlows.size()];
        
        for (int i = 0; i < cashFlows.size(); i++) {
            Map.Entry<LocalDate, Double> cashFlow = cashFlows.get(i);
            LocalDate cashFlowDate = cashFlow.getKey();
            double amount = cashFlow.getValue();
            
            if (cashFlowDate.isAfter(valuationDate)) {
                times[i] = DateUtils.yearFraction(valuationDate, cashFlowDate);
                amounts[i] = amount;
            }
        }
        
        return FinancialMathUtils.calculateConvexity(amounts, times, yield);
    }
    
    /**
     * Calculate bond price sensitivity to yield changes
     * 
     * @param bond Bond object
     * @param valuationDate Valuation date
     * @param yield Current yield
     * @param yieldChange Yield change in decimal
     * @return Expected percentage price change
     */
    public static double calculateBondPriceSensitivity(Bond bond, LocalDate valuationDate, 
                                                     double yield, double yieldChange) {
        // Calculate current price
        double currentPrice = bond.presentValue(valuationDate, Map.of("yieldRate", yield));
        
        // Calculate duration and convexity
        Map<String, Double> durations = calculateBondDurations(bond, valuationDate, yield);
        double modifiedDuration = durations.get("modifiedDuration");
        double convexity = calculateBondConvexity(bond, valuationDate, yield);
        
        // Calculate price change with convexity adjustment
        return FinancialMathUtils.calculatePriceChangeWithConvexity(
            currentPrice, modifiedDuration, convexity, yieldChange
        );
    }
    
    /**
     * Calculate key rate durations for a bond
     * 
     * @param bond Bond object
     * @param valuationDate Valuation date
     * @param yieldCurve Yield curve
     * @param keyRateTenors Tenors at which to calculate key rate durations
     * @return Map of tenors to key rate durations
     */
    public static Map<Double, Double> calculateKeyRateDurations(Bond bond, LocalDate valuationDate, 
                                                              YieldCurveModels yieldCurve, 
                                                              double[] keyRateTenors) {
        // Calculate base price
        double basePrice = calculateBondPriceWithYieldCurve(bond, valuationDate, yieldCurve);
        
        // Calculate key rate durations
        Map<Double, Double> keyRateDurations = new HashMap<>();
        double shift = 0.0001; // 1 basis point
        
        for (double tenor : keyRateTenors) {
            // Create a shocked yield curve
            LocalDate tenorDate = valuationDate.plusDays((int)(tenor * 365));
            double originalRate = yieldCurve.getRate(tenor);
            
            // Create new yield curve with shocked rate at the key tenor
            Map<Double, Double> shockedRates = new HashMap<>();
            for (double t : keyRateTenors) {
                if (Math.abs(t - tenor) < 0.0001) {
                    shockedRates.put(t, originalRate + shift);
                } else {
                    shockedRates.put(t, yieldCurve.getRate(t));
                }
            }
            
            YieldCurveModels shockedCurve = YieldCurveModels.createModel(
                valuationDate, shockedRates, "linear"
            );
            
            // Calculate price with shocked curve
            double shockedPrice = calculateBondPriceWithYieldCurve(bond, valuationDate, shockedCurve);
            
            // Calculate key rate duration
            double keyRateDuration = -(shockedPrice - basePrice) / (basePrice * shift);
            keyRateDurations.put(tenor, keyRateDuration);
        }
        
        return keyRateDurations;
    }
    
    /**
     * Calculate zero-coupon yield curve using bootstrap method
     * 
     * @param valuationDate Valuation date
     * @param bonds List of bonds for bootstrapping
     * @return Map of tenors to zero rates
     */
    public static Map<Double, Double> bootstrapYieldCurve(LocalDate valuationDate, List<Bond> bonds) {
        Map<Double, Double> zeroRates = new HashMap<>();
        
        // Sort bonds by maturity
        bonds.sort((b1, b2) -> b1.getMaturityDate().compareTo(b2.getMaturityDate()));
        
        for (Bond bond : bonds) {
            double tenor = DateUtils.yearFraction(valuationDate, bond.getMaturityDate());
            
            // For short tenors, use simple bootstrapping
            if (tenor <= 1.0) {
                // For short maturity, approximate zero rate from bond yield
                double bondYield = bond.getYield(valuationDate, bond.presentValue(valuationDate, Map.of()));
                zeroRates.put(tenor, bondYield);
            } else {
                // For longer maturities, we need to consider existing zero rates
                List<Map.Entry<LocalDate, Double>> cashFlows = calculateBondCashFlows(bond);
                double bondPrice = bond.presentValue(valuationDate, Map.of());
                
                // Discount all cash flows except the final one using existing zero rates
                double presentValueExcludingFinal = 0.0;
                double finalPayment = 0.0;
                LocalDate finalDate = null;
                
                for (Map.Entry<LocalDate, Double> cashFlow : cashFlows) {
                    LocalDate cashFlowDate = cashFlow.getKey();
                    double amount = cashFlow.getValue();
                    
                    if (cashFlowDate.isAfter(valuationDate)) {
                        double timeToPayment = DateUtils.yearFraction(valuationDate, cashFlowDate);
                        
                        if (cashFlowDate.equals(bond.getMaturityDate())) {
                            finalPayment = amount;
                            finalDate = cashFlowDate;
                        } else {
                            // Interpolate zero rate for this tenor
                            double zeroRate = interpolateZeroRate(timeToPayment, zeroRates);
                            double discountFactor = Math.exp(-zeroRate * timeToPayment);
                            presentValueExcludingFinal += amount * discountFactor;
                        }
                    }
                }
                
                // Solve for zero rate at the final payment date
                double timeToFinal = DateUtils.yearFraction(valuationDate, finalDate);
                double discountFactorFinal = (bondPrice - presentValueExcludingFinal) / finalPayment;
                double zeroRate = -Math.log(discountFactorFinal) / timeToFinal;
                
                zeroRates.put(timeToFinal, zeroRate);
            }
        }
        
        return zeroRates;
    }
    
    /**
     * Interpolate zero rate for a given tenor
     * 
     * @param tenor Tenor in years
     * @param zeroRates Map of tenors to zero rates
     * @return Interpolated zero rate
     */
    private static double interpolateZeroRate(double tenor, Map<Double, Double> zeroRates) {
        // Find surrounding tenors
        Double lowerTenor = null;
        Double upperTenor = null;
        
        for (double t : zeroRates.keySet()) {
            if (t <= tenor && (lowerTenor == null || t > lowerTenor)) {
                lowerTenor = t;
            }
            if (t >= tenor && (upperTenor == null || t < upperTenor)) {
                upperTenor = t;
            }
        }
        
        if (lowerTenor == null) {
            return zeroRates.get(upperTenor);
        }
        
        if (upperTenor == null) {
            return zeroRates.get(lowerTenor);
        }
        
        if (Math.abs(lowerTenor - tenor) < 0.00001) {
            return zeroRates.get(lowerTenor);
        }
        
        if (Math.abs(upperTenor - tenor) < 0.00001) {
            return zeroRates.get(upperTenor);
        }
        
        // Linear interpolation
        double lowerRate = zeroRates.get(lowerTenor);
        double upperRate = zeroRates.get(upperTenor);
        
        return FinancialMathUtils.linearInterpolate(tenor, lowerTenor, lowerRate, upperTenor, upperRate);
    }
    
    /**
     * Calculate the par yield for a given maturity
     * 
     * @param valuationDate Valuation date
     * @param maturityDate Maturity date
     * @param paymentsPerYear Coupon payments per year
     * @param yieldCurve Yield curve for discounting
     * @return Par yield (coupon rate that makes bond price equal to par)
     */
    public static double calculateParYield(LocalDate valuationDate, LocalDate maturityDate, 
                                         int paymentsPerYear, YieldCurveModels yieldCurve) {
        // Calculate discount factors for all payment dates
        List<Double> discountFactors = new ArrayList<>();
        
        int paymentFrequencyMonths = 12 / paymentsPerYear;
        LocalDate currentDate = valuationDate;
        
        while (currentDate.isBefore(maturityDate)) {
            currentDate = currentDate.plusMonths(paymentFrequencyMonths);
            
            if (!currentDate.isAfter(maturityDate)) {
                double timeToPayment = DateUtils.yearFraction(valuationDate, currentDate);
                double discountFactor = yieldCurve.getDiscountFactor(timeToPayment);
                discountFactors.add(discountFactor);
            }
        }
        
        // Calculate sum of discount factors
        double sumDiscountFactors = 0.0;
        for (double df : discountFactors) {
            sumDiscountFactors += df;
        }
        
        // Par yield formula: (1 - Final DF) / Sum(DF) * Payments Per Year
        double finalDF = discountFactors.get(discountFactors.size() - 1);
        return (1.0 - finalDF) / (sumDiscountFactors / paymentsPerYear);
    }
    
    /**
     * Calculate the forward rate between two dates
     * 
     * @param valuationDate Current valuation date
     * @param startDate Start date
     * @param endDate End date
     * @param yieldCurve Yield curve
     * @return Forward rate
     */
    public static double calculateForwardRate(LocalDate valuationDate, LocalDate startDate, 
                                            LocalDate endDate, YieldCurveModels yieldCurve) {
        double startTenor = DateUtils.yearFraction(valuationDate, startDate);
        double endTenor = DateUtils.yearFraction(valuationDate, endDate);
        
        return yieldCurve.getForwardRate(startTenor, endTenor);
    }
} 