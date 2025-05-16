package com.quant.finance.derivatives;

import com.quant.finance.instruments.Instrument;
import com.quant.finance.utils.DateUtils;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class representing a Credit Default Swap (CDS) contract
 */
public class CreditDefaultSwap extends Instrument {
    private final double notionalAmount;
    private final double spreadBps;  // CDS spread in basis points
    private final int paymentsPerYear;
    private final LocalDate effectiveDate;
    private final LocalDate maturityDate;
    private final String referenceEntity;
    private final List<LocalDate> paymentDates;
    private final double recoveryRate;
    
    /**
     * Constructor for a Credit Default Swap
     * 
     * @param id Unique identifier
     * @param name Descriptive name
     * @param issueDate Trade date
     * @param notionalAmount Contract notional amount
     * @param spreadBps CDS spread in basis points
     * @param paymentsPerYear Number of premium payments per year
     * @param effectiveDate Start date of protection
     * @param maturityDate End date of protection
     * @param referenceEntity Reference entity name
     * @param recoveryRate Expected recovery rate in case of default (0.0 to 1.0)
     */
    public CreditDefaultSwap(String id, String name, LocalDate issueDate, double notionalAmount,
                            double spreadBps, int paymentsPerYear, LocalDate effectiveDate,
                            LocalDate maturityDate, String referenceEntity, double recoveryRate) {
        super(id, name, issueDate);
        this.notionalAmount = notionalAmount;
        this.spreadBps = spreadBps;
        this.paymentsPerYear = paymentsPerYear;
        this.effectiveDate = effectiveDate;
        this.maturityDate = maturityDate;
        this.referenceEntity = referenceEntity;
        this.recoveryRate = recoveryRate;
        this.paymentDates = generatePaymentDates();
    }
    
    private List<LocalDate> generatePaymentDates() {
        List<LocalDate> dates = new ArrayList<>();
        LocalDate currentDate = effectiveDate;
        
        // Calculate payment frequency in months
        int monthsPerPayment = 12 / paymentsPerYear;
        
        while (currentDate.isBefore(maturityDate)) {
            currentDate = currentDate.plusMonths(monthsPerPayment);
            
            if (!currentDate.isAfter(maturityDate)) {
                dates.add(currentDate);
            }
        }
        
        return dates;
    }
    
    @Override
    public double presentValue(LocalDate valuationDate, Map<String, Object> marketData) {
        // Extract market data
        double riskFreeRate = (double) marketData.getOrDefault("riskFreeRate", 0.02);
        double impliedHazardRate = (double) marketData.getOrDefault("hazardRate", convertSpreadToHazardRate());
        
        // Calculate premium leg and protection leg values
        double premiumLegPV = calculatePremiumLegPV(valuationDate, riskFreeRate, impliedHazardRate);
        double protectionLegPV = calculateProtectionLegPV(valuationDate, riskFreeRate, impliedHazardRate);
        
        // Return value from perspective of protection buyer (pays premium, receives protection)
        return protectionLegPV - premiumLegPV;
    }
    
    /**
     * Convert CDS spread to hazard rate (simplified model)
     * Assuming flat hazard rate and risk-free rate
     */
    private double convertSpreadToHazardRate() {
        // Convert spread from bps to decimal
        double spreadDecimal = spreadBps / 10000.0;
        
        // Simple approximation: hazard rate ≈ spread / (1 - recovery rate)
        return spreadDecimal / (1 - recoveryRate);
    }
    
    /**
     * Calculate the present value of the premium leg (payments from protection buyer)
     */
    private double calculatePremiumLegPV(LocalDate valuationDate, double riskFreeRate, double hazardRate) {
        double premiumLegPV = 0.0;
        double spreadDecimal = spreadBps / 10000.0;
        
        for (int i = 0; i < paymentDates.size(); i++) {
            LocalDate paymentDate = paymentDates.get(i);
            
            if (paymentDate.isAfter(valuationDate)) {
                // Calculate time to payment
                double timeToPayment = DateUtils.yearFraction(valuationDate, paymentDate);
                
                // Determine start of period
                LocalDate periodStart = (i == 0) ? effectiveDate : paymentDates.get(i - 1);
                if (periodStart.isBefore(valuationDate)) {
                    periodStart = valuationDate;
                }
                
                double periodLength = DateUtils.yearFraction(periodStart, paymentDate);
                
                // Calculate survival probability to payment date
                double survivalProb = Math.exp(-hazardRate * timeToPayment);
                
                // Calculate discount factor
                double discountFactor = Math.exp(-riskFreeRate * timeToPayment);
                
                // Premium payment (assuming no accrual on default)
                double premium = notionalAmount * spreadDecimal * periodLength;
                
                // Premium PV = Premium * Discount Factor * Survival Probability
                premiumLegPV += premium * discountFactor * survivalProb;
            }
        }
        
        return premiumLegPV;
    }
    
    /**
     * Calculate the premium leg PV with a specific spread value
     * Used for implied spread calculations
     */
    private double calculatePremiumLegPV(LocalDate valuationDate, double riskFreeRate, 
                                        double hazardRate, double customSpreadBps) {
        double premiumLegPV = 0.0;
        double spreadDecimal = customSpreadBps / 10000.0;
        
        for (int i = 0; i < paymentDates.size(); i++) {
            LocalDate paymentDate = paymentDates.get(i);
            
            if (paymentDate.isAfter(valuationDate)) {
                // Calculate time to payment
                double timeToPayment = DateUtils.yearFraction(valuationDate, paymentDate);
                
                // Determine start of period
                LocalDate periodStart = (i == 0) ? effectiveDate : paymentDates.get(i - 1);
                if (periodStart.isBefore(valuationDate)) {
                    periodStart = valuationDate;
                }
                
                double periodLength = DateUtils.yearFraction(periodStart, paymentDate);
                
                // Calculate survival probability to payment date
                double survivalProb = Math.exp(-hazardRate * timeToPayment);
                
                // Calculate discount factor
                double discountFactor = Math.exp(-riskFreeRate * timeToPayment);
                
                // Premium payment (assuming no accrual on default)
                double premium = notionalAmount * spreadDecimal * periodLength;
                
                // Premium PV = Premium * Discount Factor * Survival Probability
                premiumLegPV += premium * discountFactor * survivalProb;
            }
        }
        
        return premiumLegPV;
    }
    
    /**
     * Calculate the present value of the protection leg (payment in case of default)
     */
    private double calculateProtectionLegPV(LocalDate valuationDate, double riskFreeRate, double hazardRate) {
        // Time to maturity in years
        double timeToMaturity = DateUtils.yearFraction(valuationDate, maturityDate);
        
        // Protection amount (loss given default)
        double lossGivenDefault = notionalAmount * (1 - recoveryRate);
        
        // For simplicity, we'll use a continuous-time model to value the protection leg
        // Integrate protection payment * discount factor * default probability density over time
        
        // Discretize the time period into monthly steps for numerical integration
        int numSteps = (int) Math.ceil(timeToMaturity * 12);
        double stepSize = timeToMaturity / numSteps;
        double protectionLegPV = 0.0;
        
        for (int i = 1; i <= numSteps; i++) {
            double time = i * stepSize;
            
            // Discount factor at time t
            double discountFactor = Math.exp(-riskFreeRate * time);
            
            // Survival probability to time t-dt
            double survivalProbBefore = Math.exp(-hazardRate * (time - stepSize));
            
            // Survival probability to time t
            double survivalProbAfter = Math.exp(-hazardRate * time);
            
            // Default probability during this small time interval
            double defaultProb = survivalProbBefore - survivalProbAfter;
            
            // Contribution to protection leg value
            protectionLegPV += lossGivenDefault * discountFactor * defaultProb;
        }
        
        return protectionLegPV;
    }
    
    @Override
    public Map<String, Double> calculateRisks(LocalDate valuationDate, Map<String, Object> marketData) {
        Map<String, Double> risks = new HashMap<>();
        
        // Extract market data
        double riskFreeRate = (double) marketData.getOrDefault("riskFreeRate", 0.02);
        double hazardRate = (double) marketData.getOrDefault("hazardRate", convertSpreadToHazardRate());
        
        // Base present value
        double basePV = presentValue(valuationDate, marketData);
        
        // Calculate credit spread sensitivity (CS01 - impact of 1bp increase in spread)
        // Create shifted market data
        Map<String, Object> shiftedMarketData = new HashMap<>(marketData);
        double shiftedHazardRate = hazardRate * (1 + 0.0001 / (spreadBps / 10000.0));
        shiftedMarketData.put("hazardRate", shiftedHazardRate);
        
        double shiftedPV = presentValue(valuationDate, shiftedMarketData);
        double cs01 = basePV - shiftedPV;
        risks.put("cs01", cs01);
        
        // Calculate interest rate sensitivity (IR01 - impact of 1bp increase in rates)
        shiftedMarketData = new HashMap<>(marketData);
        shiftedMarketData.put("riskFreeRate", riskFreeRate + 0.0001);
        
        shiftedPV = presentValue(valuationDate, shiftedMarketData);
        double ir01 = basePV - shiftedPV;
        risks.put("ir01", ir01);
        
        // Calculate risky duration - weighted average time to payment
        double premiumLegPV = calculatePremiumLegPV(valuationDate, riskFreeRate, hazardRate);
        double weightedTime = 0.0;
        
        for (int i = 0; i < paymentDates.size(); i++) {
            LocalDate paymentDate = paymentDates.get(i);
            
            if (paymentDate.isAfter(valuationDate)) {
                double timeToPayment = DateUtils.yearFraction(valuationDate, paymentDate);
                
                // Determine start of period
                LocalDate periodStart = (i == 0) ? effectiveDate : paymentDates.get(i - 1);
                if (periodStart.isBefore(valuationDate)) {
                    periodStart = valuationDate;
                }
                
                double periodLength = DateUtils.yearFraction(periodStart, paymentDate);
                double spreadDecimal = spreadBps / 10000.0;
                double premium = notionalAmount * spreadDecimal * periodLength;
                double survivalProb = Math.exp(-hazardRate * timeToPayment);
                double discountFactor = Math.exp(-riskFreeRate * timeToPayment);
                double paymentPV = premium * discountFactor * survivalProb;
                
                weightedTime += timeToPayment * paymentPV;
            }
        }
        
        double riskyDuration = weightedTime / premiumLegPV;
        risks.put("riskyDuration", riskyDuration);
        
        // Recovery rate sensitivity - impact of 1% change in recovery rate
        shiftedMarketData = new HashMap<>(marketData);
        CreditDefaultSwap shiftedCDS = new CreditDefaultSwap(
            id, name, issueDate, notionalAmount, spreadBps, paymentsPerYear,
            effectiveDate, maturityDate, referenceEntity, recoveryRate + 0.01
        );
        
        shiftedPV = shiftedCDS.presentValue(valuationDate, shiftedMarketData);
        double recoveryRateSensitivity = (basePV - shiftedPV) / 0.01;
        risks.put("recoveryRateSensitivity", recoveryRateSensitivity);
        
        return risks;
    }
    
    /**
     * Calculate the implied CDS spread from market data
     * 
     * @param valuationDate Valuation date
     * @param marketData Market data including hazard rate
     * @return Implied CDS spread in basis points
     */
    public double calculateImpliedSpread(LocalDate valuationDate, Map<String, Object> marketData) {
        double riskFreeRate = (double) marketData.getOrDefault("riskFreeRate", 0.02);
        double hazardRate = (double) marketData.get("hazardRate");
        
        if (hazardRate <= 0) {
            throw new IllegalArgumentException("Hazard rate must be positive");
        }
        
        // Calculate protection leg PV for a unit spread
        double protectionLegPV = calculateProtectionLegPV(valuationDate, riskFreeRate, hazardRate);
        
        // Calculate premium leg PV for a 1bp spread without modifying the field
        double premiumLegPVper1bp = calculatePremiumLegPV(valuationDate, riskFreeRate, hazardRate, 1.0);
        
        // Implied spread = (Protection Leg PV) / (Premium Leg PV per bp)
        return protectionLegPV / premiumLegPVper1bp;
    }
    
    // Getters
    public double getNotionalAmount() {
        return notionalAmount;
    }
    
    public double getSpreadBps() {
        return spreadBps;
    }
    
    public int getPaymentsPerYear() {
        return paymentsPerYear;
    }
    
    public LocalDate getEffectiveDate() {
        return effectiveDate;
    }
    
    public LocalDate getMaturityDate() {
        return maturityDate;
    }
    
    public String getReferenceEntity() {
        return referenceEntity;
    }
    
    public List<LocalDate> getPaymentDates() {
        return new ArrayList<>(paymentDates);
    }
    
    public double getRecoveryRate() {
        return recoveryRate;
    }
} 