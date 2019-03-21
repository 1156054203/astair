from __future__ import division

def non_zero_division(x, y):
    """Ensures safe division of x by y, assuming that the result will be zero if y is 0."""
    if y == 0:
        return int(0)
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
