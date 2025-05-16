package com.quant.finance.instruments;

import com.quant.finance.models.YieldCurve;
import com.quant.finance.utils.DateUtils;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class representing an Interest Rate Swap
 */
public class InterestRateSwap extends Instrument {
    private final double notionalAmount;
    private final double fixedRate;
    private final String floatingRateIndex; // e.g., "LIBOR", "SOFR"
    private final double floatingRateSpread;
    private final int paymentsPerYear;
    private final LocalDate effectiveDate;
    private final LocalDate maturityDate;
    private final List<LocalDate> paymentDates;
    
    public InterestRateSwap(String id, String name, LocalDate issueDate, 
                           double notionalAmount, double fixedRate, String floatingRateIndex,
                           double floatingRateSpread, int paymentsPerYear, 
                           LocalDate effectiveDate, LocalDate maturityDate) {
        super(id, name, issueDate);
        this.notionalAmount = notionalAmount;
        this.fixedRate = fixedRate;
        this.floatingRateIndex = floatingRateIndex;
        this.floatingRateSpread = floatingRateSpread;
        this.paymentsPerYear = paymentsPerYear;
        this.effectiveDate = effectiveDate;
        this.maturityDate = maturityDate;
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
        // Extract yield curves from market data
        YieldCurve discountCurve = (YieldCurve) marketData.get("discountCurve");
        YieldCurve forwardCurve = (YieldCurve) marketData.get("forwardCurve");
        
        if (discountCurve == null || forwardCurve == null) {
            throw new IllegalArgumentException("Market data must contain both discount and forward curves");
        }
        
        // Calculate fixed and floating legs
        double fixedLegPV = calculateFixedLegPV(valuationDate, discountCurve);
        double floatingLegPV = calculateFloatingLegPV(valuationDate, discountCurve, forwardCurve);
        
        // Swap value from the perspective of the fixed rate payer
        return floatingLegPV - fixedLegPV;
    }
    
    private double calculateFixedLegPV(LocalDate valuationDate, YieldCurve discountCurve) {
        double fixedLegPV = 0.0;
        double fixedPayment = notionalAmount * fixedRate / paymentsPerYear;
        
        for (LocalDate paymentDate : paymentDates) {
            if (paymentDate.isAfter(valuationDate)) {
                double timeToPayment = DateUtils.yearFraction(valuationDate, paymentDate);
                double discountFactor = discountCurve.getDiscountFactor(timeToPayment);
                fixedLegPV += fixedPayment * discountFactor;
            }
        }
        
        return fixedLegPV;
    }
    
    private double calculateFloatingLegPV(LocalDate valuationDate, YieldCurve discountCurve, YieldCurve forwardCurve) {
        double floatingLegPV = 0.0;
        
        // For the floating leg, we need to project the forward rates for each period
        for (int i = 0; i < paymentDates.size(); i++) {
            LocalDate paymentDate = paymentDates.get(i);
            
            if (paymentDate.isAfter(valuationDate)) {
                // Determine start date for this payment period
                LocalDate startDate = (i == 0) ? effectiveDate : paymentDates.get(i - 1);
                if (startDate.isBefore(valuationDate)) {
                    startDate = valuationDate; // For the current period, use valuation date
                }
                
                double periodLength = DateUtils.yearFraction(startDate, paymentDate);
                double timeToPayment = DateUtils.yearFraction(valuationDate, paymentDate);
                
                // Calculate forward rate for this period
                double forwardRate;
                if (startDate.equals(valuationDate)) {
                    // If startDate is valuationDate, use direct lookup
                    forwardRate = forwardCurve.getRate(periodLength);
                } else {
                    // Otherwise, calculate proper forward rate between startDate and paymentDate
                    double timeToStart = DateUtils.yearFraction(valuationDate, startDate);
                    forwardRate = forwardCurve.getForwardRate(timeToStart, timeToPayment);
                }
                
                // Add spread to forward rate
                double projectedRate = forwardRate + floatingRateSpread;
                
                // Calculate projected payment
                double floatingPayment = notionalAmount * projectedRate * periodLength;
                
                // Discount the payment
                double discountFactor = discountCurve.getDiscountFactor(timeToPayment);
                floatingLegPV += floatingPayment * discountFactor;
            }
        }
        
        return floatingLegPV;
    }
    
    @Override
    public Map<String, Double> calculateRisks(LocalDate valuationDate, Map<String, Object> marketData) {
        Map<String, Double> risks = new HashMap<>();
        
        // Extract yield curves from market data
        YieldCurve discountCurve = (YieldCurve) marketData.get("discountCurve");
        YieldCurve forwardCurve = (YieldCurve) marketData.get("forwardCurve");
        
        // Calculate DV01 (Dollar Value of 01) - sensitivity to 1bp parallel shift in rates
        double basePV = presentValue(valuationDate, marketData);
        
        // Create shifted discount curve (up 1bp)
        Map<Double, Double> shiftedRates = new HashMap<>(discountCurve.getTenorRates());
        shiftedRates.replaceAll((tenor, rate) -> rate + 0.0001); // 1bp = 0.0001
        YieldCurve shiftedDiscountCurve = new YieldCurve(valuationDate, shiftedRates);
        
        // Replace discount curve with shifted curve
        Map<String, Object> shiftedMarketData = new HashMap<>(marketData);
        shiftedMarketData.put("discountCurve", shiftedDiscountCurve);
        
        double shiftedPV = presentValue(valuationDate, shiftedMarketData);
        double dv01 = basePV - shiftedPV;
        risks.put("dv01", dv01);
        
        // Calculate effective duration
        risks.put("effectiveDuration", dv01 / basePV * 10000);
        
        // Calculate floating rate sensitivity
        // Create shifted forward curve (up 1bp)
        Map<Double, Double> shiftedForwardRates = new HashMap<>(forwardCurve.getTenorRates());
        shiftedForwardRates.replaceAll((tenor, rate) -> rate + 0.0001);
        YieldCurve shiftedForwardCurve = new YieldCurve(valuationDate, shiftedForwardRates);
        
        // Replace forward curve with shifted curve
        shiftedMarketData = new HashMap<>(marketData);
        shiftedMarketData.put("forwardCurve", shiftedForwardCurve);
        
        shiftedPV = presentValue(valuationDate, shiftedMarketData);
        double floatingRateSensitivity = basePV - shiftedPV;
        risks.put("floatingRateSensitivity", floatingRateSensitivity);
        
        return risks;
    }
    
    /**
     * Calculate the par swap rate - the fixed rate at which the swap's present value is zero
     */
    public double calculateParSwapRate(LocalDate valuationDate, Map<String, Object> marketData) {
        // Extract yield curves from market data
        YieldCurve discountCurve = (YieldCurve) marketData.get("discountCurve");
        YieldCurve forwardCurve = (YieldCurve) marketData.get("forwardCurve");
        
        if (discountCurve == null || forwardCurve == null) {
            throw new IllegalArgumentException("Market data must contain both discount and forward curves");
        }
        
        // Calculate floating leg PV with the current fixed rate
        double floatingLegPV = calculateFloatingLegPV(valuationDate, discountCurve, forwardCurve);
        
        // Calculate PV of a single basis point on the fixed leg
        double annuityFactor = 0.0;
        for (LocalDate paymentDate : paymentDates) {
            if (paymentDate.isAfter(valuationDate)) {
                double timeToPayment = DateUtils.yearFraction(valuationDate, paymentDate);
                double discountFactor = discountCurve.getDiscountFactor(timeToPayment);
                annuityFactor += discountFactor / paymentsPerYear;
            }
        }
        
        // Par swap rate = Floating Leg PV / (Notional * Annuity Factor)
        return floatingLegPV / (notionalAmount * annuityFactor);
    }
    
    // Getters
    public double getNotionalAmount() {
        return notionalAmount;
    }
    
    public double getFixedRate() {
        return fixedRate;
    }
    
    public String getFloatingRateIndex() {
        return floatingRateIndex;
    }
    
    public double getFloatingRateSpread() {
        return floatingRateSpread;
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
    
    public List<LocalDate> getPaymentDates() {
        return new ArrayList<>(paymentDates);
    }
} 