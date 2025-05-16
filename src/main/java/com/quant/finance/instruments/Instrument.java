package com.quant.finance.instruments;

import java.time.LocalDate;
import java.util.HashMap;
import java.util.Map;

/**
 * Base class for all financial instruments
 */
public abstract class Instrument {
    protected String id;
    protected String name;
    protected LocalDate issueDate;
    protected Map<String, Object> properties;
    
    public Instrument(String id, String name, LocalDate issueDate) {
        this.id = id;
        this.name = name;
        this.issueDate = issueDate;
        this.properties = new HashMap<>();
    }
    
    /**
     * Calculate the present value of the instrument
     * 
     * @param valuationDate Date on which to perform valuation
     * @param marketData Relevant market data for valuation
     * @return The present value
     */
    public abstract double presentValue(LocalDate valuationDate, Map<String, Object> marketData);
    
    /**
     * Calculate risks associated with this instrument (greeks, etc.)
     * 
     * @param valuationDate Date on which to perform risk calculations
     * @param marketData Relevant market data
     * @return Map of risk metrics
     */
    public abstract Map<String, Double> calculateRisks(LocalDate valuationDate, Map<String, Object> marketData);
    
    // Getters and setters
    public String getId() {
        return id;
    }
    
    public String getName() {
        return name;
    }
    
    public LocalDate getIssueDate() {
        return issueDate;
    }
    
    public void setProperty(String key, Object value) {
        properties.put(key, value);
    }
    
    public Object getProperty(String key) {
        return properties.get(key);
    }
} 