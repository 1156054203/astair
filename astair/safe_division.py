from __future__ import division

def non_zero_division(x, y, sign):
    """Ensures safe division of x by y, assuming that the result will be zero if y is 0."""
    if y == 0 and x!=0:
        return int(0)
    elif y == 0 and x==0:
        return sign
    else:
        return x / y


def non_zero_division_NA(x, y):
    """Ensures safe division of x by y, assuming that the result will be zero if x is 0, and NA if y is zero."""
    if y == 0:
        return "NA"
    elif x == 0:
        return int(0)
    else:
        return x / y


def safe_rounder(data, precision, multiple):
    """Ensures no attempts at rounding non-numeric data are done."""
    if isinstance(data, str):
        return data
    else:
        if multiple == True:
            return round(data*100, precision)
        else:
            return round(data, precision)
