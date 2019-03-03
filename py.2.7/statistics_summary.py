from __future__ import division

import math
from numpy import percentile


def general_statistics_summary(raw_data_list):
    """Gives a statistics summary about the input data."""
    quantile25, median, quantile75 = percentile(raw_data_list, [25, 50, 75])
    mu = round(sum(raw_data_list) / len(raw_data_list), 3)
    sigma = round(math.sqrt((sum([(x - mu) ** 2 for x in raw_data_list]) / (len(raw_data_list) - 1))), 3)
    return mu, round(median, 3), sigma, round(quantile25, 3), round(quantile75, 3), round(min(raw_data_list),3), round(max(raw_data_list),3)

