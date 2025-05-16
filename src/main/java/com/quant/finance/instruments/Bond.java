package com.quant.finance.instruments;

import com.quant.finance.utils.DateUtils;

import java.time.LocalDate;
import java.time.temporal.ChronoUnit;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Represents a fixed-income bond instrument
 */
public class Bond extends Instrument {
    private final double principal;
    private final double couponRate;
    private final int paymentsPerYear;
    private final LocalDate maturityDate;
    private final List<LocalDate> cashflowDates;
    
    public Bond(String id, String name, LocalDate issueDate, double principal, double couponRate, 
                int paymentsPerYear, LocalDate maturityDate) {
        super(id, name, issueDate);
        this.principal = principal;
        this.couponRate = couponRate;
        this.paymentsPerYear = paymentsPerYear;
        this.maturityDate = maturityDate;
        this.cashflowDates = generateCashflowDates();
    }
    
    private List<LocalDate> generateCashflowDates() {
        List<LocalDate> dates = new ArrayList<>();
        LocalDate currentDate = issueDate;
        
        while (currentDate.isBefore(maturityDate)) {
            // Add months based on payment frequency
            int monthsToAdd = 12 / paymentsPerYear;
            currentDate = currentDate.plusMonths(monthsToAdd);
            
            if (!currentDate.isAfter(maturityDate)) {
                dates.add(currentDate);
            }
        }
        
        return dates;
    }
    
    /**
     * Calculate present value using discounted cash flow method
     */
    @Override
    public double presentValue(LocalDate valuationDate, Map<String, Object> marketData) {
        // Extract yield curve from market data
        double yieldRate = (double) marketData.getOrDefault("yieldRate", 0.05);
        
        // Calculate present value
        double pv = 0.0;
        
        // Add coupon payments
        double couponAmount = principal * couponRate / paymentsPerYear;
        
        for (LocalDate cashflowDate : cashflowDates) {
            if (cashflowDate.isAfter(valuationDate)) {
                double timeInYears = DateUtils.yearFraction(valuationDate, cashflowDate);
                double discountFactor = 1.0 / Math.pow(1.0 + yieldRate, timeInYears);
                pv += couponAmount * discountFactor;
            }
        }
        
        // Add principal repayment at maturity
        if (maturityDate.isAfter(valuationDate)) {
            double timeToMaturity = DateUtils.yearFraction(valuationDate, maturityDate);
            double discountFactor = 1.0 / Math.pow(1.0 + yieldRate, timeToMaturity);
            pv += principal * discountFactor;
        }
        
        return pv;
    }
    
    /**
     * Calculate key bond risk metrics (duration, convexity, etc.)
     */
    @Override
    public Map<String, Double> calculateRisks(LocalDate valuationDate, Map<String, Object> marketData) {
        Map<String, Double> risks = new HashMap<>();
        double yieldRate = (double) marketData.getOrDefault("yieldRate", 0.05);
        
        // Calculate Macaulay Duration
        double duration = calculateMacaulayDuration(valuationDate, yieldRate);
        risks.put("macaulayDuration", duration);
        
        // Calculate Modified Duration
        double modifiedDuration = duration / (1 + yieldRate / paymentsPerYear);
        risks.put("modifiedDuration", modifiedDuration);
        
        // Calculate approximate price change for 1% yield change
        double priceChange = -modifiedDuration * 0.01;
        risks.put("priceChangeFor1PercentYieldChange", priceChange);
        
        // Calculate convexity (simplified)
        double convexity = calculateConvexity(valuationDate, yieldRate);
        risks.put("convexity", convexity);
        
        return risks;
    }
    
    private double calculateMacaulayDuration(LocalDate valuationDate, double yieldRate) {
        double pv = presentValue(valuationDate, Map.of("yieldRate", yieldRate));
        double weightedTime = 0.0;
        double couponAmount = principal * couponRate / paymentsPerYear;
        
        for (LocalDate cashflowDate : cashflowDates) {
            if (cashflowDate.isAfter(valuationDate)) {
                double timeInYears = DateUtils.yearFraction(valuationDate, cashflowDate);
                double discountFactor = 1.0 / Math.pow(1.0 + yieldRate, timeInYears);
                weightedTime += timeInYears * couponAmount * discountFactor;
            }
        }
        
        // Principal repayment
        double timeToMaturity = DateUtils.yearFraction(valuationDate, maturityDate);
        double discountFactor = 1.0 / Math.pow(1.0 + yieldRate, timeToMaturity);
        weightedTime += timeToMaturity * principal * discountFactor;
        
        return weightedTime / pv;
    }
    
    private double calculateConvexity(LocalDate valuationDate, double yieldRate) {
        double pv = presentValue(valuationDate, Map.of("yieldRate", yieldRate));
        double weightedTimeSquared = 0.0;
        double couponAmount = principal * couponRate / paymentsPerYear;
        
        for (LocalDate cashflowDate : cashflowDates) {
            if (cashflowDate.isAfter(valuationDate)) {
                double timeInYears = DateUtils.yearFraction(valuationDate, cashflowDate);
                double discountFactor = 1.0 / Math.pow(1.0 + yieldRate, timeInYears);
                weightedTimeSquared += timeInYears * (timeInYears + 1) * couponAmount * discountFactor;
            }
        }
        
        // Principal repayment
        double timeToMaturity = DateUtils.yearFraction(valuationDate, maturityDate);
        double discountFactor = 1.0 / Math.pow(1.0 + yieldRate, timeToMaturity);
        weightedTimeSquared += timeToMaturity * (timeToMaturity + 1) * principal * discountFactor;
        
        return weightedTimeSquared / (pv * Math.pow(1 + yieldRate, 2));
    }
    
    // Getters
    public double getPrincipal() {
        return principal;
    }
    
    public double getCouponRate() {
        return couponRate;
    }
    
    public int getPaymentsPerYear() {
        return paymentsPerYear;
    }
    
    public LocalDate getMaturityDate() {
        return maturityDate;
    }
    
    public List<LocalDate> getCashflowDates() {
        return new ArrayList<>(cashflowDates);
    }
    
    public double getYield(LocalDate valuationDate, double marketPrice) {
        // Simple yield calculation using Newton-Raphson method
        double guess = couponRate; // Initial guess = coupon rate
        double epsilon = 0.0000001;
        int maxIterations = 100;
        
        for (int i = 0; i < maxIterations; i++) {
            double price = presentValue(valuationDate, Map.of("yieldRate", guess));
            double diff = price - marketPrice;
            
            if (Math.abs(diff) < epsilon) {
                return guess;
            }
            
            // Approximate derivative using small change in yield
            double delta = 0.0001;
            double priceUp = presentValue(valuationDate, Map.of("yieldRate", guess + delta));
            double derivative = (priceUp - price) / delta;
            
            // Newton-Raphson update
            guess = guess - diff / derivative;
        }
        
        return guess; // Return best guess after max iterations
    }
} 