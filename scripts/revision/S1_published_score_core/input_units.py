"""Fixed unit conversion; no file access."""
import numpy as np

def canonical(variable,value,unit):
    try:x=float(value)
    except (TypeError,ValueError):return np.nan
    if not np.isfinite(x) or x<0 or (variable=='creatinine' and x<=0):return np.nan
    u=str(unit).strip().lower().replace('μ','u').replace('µ','u')
    if variable in {'sodium','potassium','chloride','bicarbonate'} and u in {'meq/l','mmol/l'}:return x
    if variable in {'bun','creatinine'} and u=='mg/dl':return x
    if variable=='creatinine' and u in {'umol/l','micromol/l'}:return x/88.4
    if variable=='hemoglobin' and u=='g/dl':return x
    return np.nan
