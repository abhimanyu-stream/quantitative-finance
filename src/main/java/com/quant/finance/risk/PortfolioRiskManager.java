package com.quant.finance.risk;

import com.quant.finance.instruments.Instrument;
import com.quant.finance.statistics.PortfolioStatistics;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

/**
 * Class for portfolio risk management, aggregating risks across instruments
 * and performing portfolio-level risk calculations
 */
public class PortfolioRiskManager {
    private final List<Instrument> instruments;
    private final Map<String, Double> positions; // instrument ID to position size
    private final Map<String, Map<String, Double>> instrumentRisks; // instrument ID to risk metrics
    
    /**
     * Constructor for a portfolio risk manager
     */
    public PortfolioRiskManager() {
        this.instruments = new ArrayList<>();
        this.positions = new HashMap<>();
        this.instrumentRisks = new HashMap<>();
    }
    
    /**
     * Add an instrument to the portfolio
     * 
     * @param instrument Instrument to add
     * @param positionSize Size of the position
     */
    public void addPosition(Instrument instrument, double positionSize) {
        instruments.add(instrument);
        positions.put(instrument.getId(), positionSize);
    }
    
    /**
     * Calculate present value of the entire portfolio
     * 
     * @param valuationDate Valuation date
     * @param marketData Market data for valuation
     * @return Total portfolio value
     */
    public double calculatePortfolioValue(LocalDate valuationDate, Map<String, Object> marketData) {
        double totalValue = 0.0;
        
        for (Instrument instrument : instruments) {
            double instrumentValue = instrument.presentValue(valuationDate, marketData);
            double positionSize = positions.get(instrument.getId());
            totalValue += instrumentValue * positionSize;
        }
        
        return totalValue;
    }
    
    /**
     * Calculate all risk metrics for each instrument in the portfolio
     * 
     * @param valuationDate Valuation date
     * @param marketData Market data for valuation
     */
    public void calculateInstrumentRisks(LocalDate valuationDate, Map<String, Object> marketData) {
        for (Instrument instrument : instruments) {
            Map<String, Double> risks = instrument.calculateRisks(valuationDate, marketData);
            instrumentRisks.put(instrument.getId(), risks);
        }
    }
    
    /**
     * Calculate portfolio-level interest rate sensitivity (DV01)
     * 
     * @return Total DV01 for the portfolio
     */
    public double calculatePortfolioDV01() {
        double totalDV01 = 0.0;
        
        for (Instrument instrument : instruments) {
            String instrumentId = instrument.getId();
            double positionSize = positions.get(instrumentId);
            
            Map<String, Double> risks = instrumentRisks.get(instrumentId);
            if (risks != null && risks.containsKey("dv01")) {
                totalDV01 += risks.get("dv01") * positionSize;
            }
        }
        
        return totalDV01;
    }
    
    /**
     * Calculate Value at Risk (VaR) for the portfolio using historical simulation
     * 
     * @param historicalReturns Map of instrument IDs to their historical return series
     * @param confidenceLevel Confidence level (e.g., 0.95 for 95%)
     * @param timeHorizon Time horizon in days
     * @return Portfolio VaR
     */
    public double calculateHistoricalVaR(Map<String, double[]> historicalReturns, 
                                        double confidenceLevel, int timeHorizon) {
        // Get instrument IDs in consistent order
        List<String> instrumentIds = instruments.stream()
                                   .map(Instrument::getId)
                                   .collect(Collectors.toList());
        
        // Create portfolio return series
        List<Double> portfolioReturns = new ArrayList<>();
        
        // For each historical date
        int numDays = historicalReturns.get(instrumentIds.get(0)).length;
        
        for (int day = 0; day < numDays; day++) {
            double portfolioReturn = 0.0;
            
            // Add weighted return from each instrument
            for (String id : instrumentIds) {
                double positionSize = positions.get(id);
                double weight = positionSize / calculateTotalPositionSize();
                double instrumentReturn = historicalReturns.get(id)[day];
                portfolioReturn += weight * instrumentReturn;
            }
            
            portfolioReturns.add(portfolioReturn);
        }
        
        // Convert to array for statistics calculations
        double[] portfolioReturnArray = portfolioReturns.stream()
                                      .mapToDouble(Double::doubleValue)
                                      .toArray();
        
        // Calculate portfolio volatility
        DescriptiveStatistics stats = new DescriptiveStatistics(portfolioReturnArray);
        double portfolioVolatility = stats.getStandardDeviation();
        
        // Scale volatility to time horizon
        double scaledVolatility = portfolioVolatility * Math.sqrt(timeHorizon);
        
        // Calculate current portfolio value
        double portfolioValue = calculateTotalPositionSize();
        
        // Calculate VaR
        double zScore = getZScore(confidenceLevel);
        double var = portfolioValue * scaledVolatility * zScore;
        
        return var;
    }
    
    /**
     * Calculate portfolio Value at Risk using parametric method (variance-covariance)
     * 
     * @param volatilities Map of instrument IDs to their volatilities
     * @param correlationMatrix Correlation matrix between instruments
     * @param confidenceLevel Confidence level (e.g., 0.95 for 95%)
     * @param timeHorizon Time horizon in days
     * @return Portfolio VaR
     */
    public double calculateParametricVaR(Map<String, Double> volatilities,
                                        Map<String, Map<String, Double>> correlationMatrix,
                                        double confidenceLevel, int timeHorizon) {
        // Get instrument IDs in consistent order
        List<String> instrumentIds = instruments.stream()
                                   .map(Instrument::getId)
                                   .collect(Collectors.toList());
        
        // Calculate portfolio weights
        double[] weights = new double[instrumentIds.size()];
        double totalSize = calculateTotalPositionSize();
        
        for (int i = 0; i < instrumentIds.size(); i++) {
            weights[i] = positions.get(instrumentIds.get(i)) / totalSize;
        }
        
        // Create volatility vector
        double[] vols = new double[instrumentIds.size()];
        for (int i = 0; i < instrumentIds.size(); i++) {
            vols[i] = volatilities.get(instrumentIds.get(i));
        }
        
        // Create correlation matrix
        double[][] corrMatrix = new double[instrumentIds.size()][instrumentIds.size()];
        for (int i = 0; i < instrumentIds.size(); i++) {
            for (int j = 0; j < instrumentIds.size(); j++) {
                String id1 = instrumentIds.get(i);
                String id2 = instrumentIds.get(j);
                
                if (i == j) {
                    corrMatrix[i][j] = 1.0;
                } else {
                    corrMatrix[i][j] = correlationMatrix.get(id1).get(id2);
                }
            }
        }
        
        // Convert to covariance matrix
        double[][] covMatrix = new double[instrumentIds.size()][instrumentIds.size()];
        for (int i = 0; i < instrumentIds.size(); i++) {
            for (int j = 0; j < instrumentIds.size(); j++) {
                covMatrix[i][j] = corrMatrix[i][j] * vols[i] * vols[j];
            }
        }
        
        // Calculate portfolio variance
        RealMatrix weightMatrix = new Array2DRowRealMatrix(weights);
        RealMatrix covarianceMatrix = new Array2DRowRealMatrix(covMatrix);
        
        double portfolioVariance = weightMatrix.transpose()
                                  .multiply(covarianceMatrix)
                                  .multiply(weightMatrix)
                                  .getEntry(0, 0);
        
        double portfolioVolatility = Math.sqrt(portfolioVariance);
        
        // Scale to time horizon
        double scaledVolatility = portfolioVolatility * Math.sqrt(timeHorizon);
        
        // Calculate VaR
        double zScore = getZScore(confidenceLevel);
        double portfolioValue = calculateTotalPositionSize();
        double var = portfolioValue * scaledVolatility * zScore;
        
        return var;
    }
    
    /**
     * Calculate Conditional Value at Risk (Expected Shortfall)
     * 
     * @param historicalReturns Map of instrument IDs to their historical return series
     * @param confidenceLevel Confidence level (e.g., 0.95 for 95%)
     * @return Conditional Value at Risk
     */
    public double calculateCVaR(Map<String, double[]> historicalReturns, double confidenceLevel) {
        // Create weighted portfolio returns
        List<String> instrumentIds = instruments.stream()
                                   .map(Instrument::getId)
                                   .collect(Collectors.toList());
        
        int numDays = historicalReturns.get(instrumentIds.get(0)).length;
        double[] portfolioReturns = new double[numDays];
        
        for (int day = 0; day < numDays; day++) {
            double portfolioReturn = 0.0;
            
            for (String id : instrumentIds) {
                double positionSize = positions.get(id);
                double weight = positionSize / calculateTotalPositionSize();
                double instrumentReturn = historicalReturns.get(id)[day];
                portfolioReturn += weight * instrumentReturn;
            }
            
            portfolioReturns[day] = portfolioReturn;
        }
        
        // Calculate CVaR from portfolio returns
        double portfolioValue = calculateTotalPositionSize();
        return PortfolioStatistics.calculateHistoricalCVaR(portfolioReturns, confidenceLevel, portfolioValue);
    }
    
    /**
     * Calculate stress test impact based on defined market shocks
     * 
     * @param valuationDate Current valuation date
     * @param baseMarketData Base market data
     * @param stressScenarios Map of stress scenario name to market data adjustments
     * @return Map of scenario name to P&L impact
     */
    public Map<String, Double> calculateStressTests(LocalDate valuationDate, 
                                                  Map<String, Object> baseMarketData,
                                                  Map<String, Map<String, Object>> stressScenarios) {
        // Calculate base portfolio value
        double baseValue = calculatePortfolioValue(valuationDate, baseMarketData);
        
        // Results map
        Map<String, Double> stressResults = new HashMap<>();
        
        // For each stress scenario
        for (Map.Entry<String, Map<String, Object>> scenario : stressScenarios.entrySet()) {
            // Create stressed market data
            Map<String, Object> stressedMarketData = new HashMap<>(baseMarketData);
            
            // Apply stress adjustments
            for (Map.Entry<String, Object> adjustment : scenario.getValue().entrySet()) {
                String key = adjustment.getKey();
                Object value = adjustment.getValue();
                
                if (baseMarketData.containsKey(key)) {
                    // If the key already exists, apply the adjustment
                    if (baseMarketData.get(key) instanceof Double && value instanceof Double) {
                        // For numeric values, apply percentage or absolute changes
                        double baseVal = (double) baseMarketData.get(key);
                        double adjustmentValue = (double) value;
                        
                        // Determine if it's a percentage or absolute change
                        if (key.endsWith("_pct")) {
                            // Percentage change
                            String actualKey = key.substring(0, key.length() - 4);
                            stressedMarketData.put(actualKey, baseVal * (1 + adjustmentValue));
                        } else {
                            // Absolute change
                            stressedMarketData.put(key, baseVal + adjustmentValue);
                        }
                    } else {
                        // For non-numeric values, just replace
                        stressedMarketData.put(key, value);
                    }
                } else {
                    // If the key doesn't exist, just add it
                    stressedMarketData.put(key, value);
                }
            }
            
            // Calculate stressed portfolio value
            double stressedValue = calculatePortfolioValue(valuationDate, stressedMarketData);
            
            // Calculate P&L impact
            double plImpact = stressedValue - baseValue;
            stressResults.put(scenario.getKey(), plImpact);
        }
        
        return stressResults;
    }
    
    /**
     * Calculate the marginal contribution to risk for each instrument
     * 
     * @param riskMeasure Risk measure to analyze (e.g., "dv01", "var")
     * @return Map of instrument ID to marginal risk contribution
     */
    public Map<String, Double> calculateRiskContribution(String riskMeasure) {
        Map<String, Double> contributions = new HashMap<>();
        
        // Total portfolio risk
        double totalRisk = 0.0;
        
        if (riskMeasure.equals("dv01")) {
            totalRisk = calculatePortfolioDV01();
        } else {
            // Use the risk values already calculated for each instrument
            for (Instrument instrument : instruments) {
                String id = instrument.getId();
                double positionSize = positions.get(id);
                
                Map<String, Double> risks = instrumentRisks.get(id);
                if (risks != null && risks.containsKey(riskMeasure)) {
                    totalRisk += risks.get(riskMeasure) * positionSize;
                }
            }
        }
        
        // Calculate contribution for each instrument
        for (Instrument instrument : instruments) {
            String id = instrument.getId();
            double positionSize = positions.get(id);
            
            Map<String, Double> risks = instrumentRisks.get(id);
            if (risks != null && risks.containsKey(riskMeasure)) {
                double instrumentRisk = risks.get(riskMeasure);
                double contribution = (instrumentRisk * positionSize) / totalRisk;
                contributions.put(id, contribution);
            } else {
                contributions.put(id, 0.0);
            }
        }
        
        return contributions;
    }
    
    /**
     * Calculate risk budgeting - optimal position sizes to achieve target risk allocation
     * 
     * @param targetRiskAllocations Map of instrument ID to target risk allocation (as percentage)
     * @param riskMeasure Risk measure to use (e.g., "dv01", "var")
     * @return Map of instrument ID to recommended position size
     */
    public Map<String, Double> calculateRiskBudgeting(Map<String, Double> targetRiskAllocations,
                                                    String riskMeasure) {
        Map<String, Double> recommendedPositions = new HashMap<>();
        double totalTargetAllocation = targetRiskAllocations.values().stream()
                                     .mapToDouble(Double::doubleValue)
                                     .sum();
        
        // If total allocation doesn't sum to 1, normalize it
        if (Math.abs(totalTargetAllocation - 1.0) > 0.000001) {
            for (String id : targetRiskAllocations.keySet()) {
                targetRiskAllocations.put(id, targetRiskAllocations.get(id) / totalTargetAllocation);
            }
        }
        
        // Get current risk sensitivities for each instrument
        Map<String, Double> riskSensitivities = new HashMap<>();
        for (Instrument instrument : instruments) {
            String id = instrument.getId();
            Map<String, Double> risks = instrumentRisks.get(id);
            
            if (risks != null && risks.containsKey(riskMeasure)) {
                riskSensitivities.put(id, risks.get(riskMeasure));
            } else {
                riskSensitivities.put(id, 0.0);
            }
        }
        
        // Calculate position adjustment factors
        double totalExposureNeeded = calculateTotalPositionSize();
        
        for (Instrument instrument : instruments) {
            String id = instrument.getId();
            
            if (riskSensitivities.get(id) > 0 && targetRiskAllocations.containsKey(id)) {
                double targetRisk = targetRiskAllocations.get(id);
                double riskPerUnit = riskSensitivities.get(id);
                
                // Calculate position size that would give the target risk allocation
                double recommendedSize = (totalExposureNeeded * targetRisk) / riskPerUnit;
                recommendedPositions.put(id, recommendedSize);
            } else {
                // Keep current position if no risk sensitivity or no target allocation
                recommendedPositions.put(id, positions.get(id));
            }
        }
        
        return recommendedPositions;
    }
    
    /**
     * Get the Z-score for a given confidence level
     * 
     * @param confidenceLevel Confidence level (e.g., 0.95 for 95%)
     * @return Z-score
     */
    private double getZScore(double confidenceLevel) {
        // Simplified approach - for a more accurate calculation, use a normal distribution function
        if (confidenceLevel == 0.95) {
            return 1.645;
        } else if (confidenceLevel == 0.99) {
            return 2.326;
        } else if (confidenceLevel == 0.975) {
            return 1.96;
        } else if (confidenceLevel == 0.90) {
            return 1.282;
        } else {
            // Default approximation
            return 1.645;
        }
    }
    
    /**
     * Calculate total position size across all instruments
     * 
     * @return Total position size
     */
    private double calculateTotalPositionSize() {
        return positions.values().stream()
                      .mapToDouble(Double::doubleValue)
                      .sum();
    }
} 