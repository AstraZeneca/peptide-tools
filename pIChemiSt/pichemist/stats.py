import math


def mean(input_list):
    """Calculates the mean for an input list."""
    return sum(input_list) / len(input_list)


def stddev(input_list):
    """Calculates the standard deviation for an input list."""
    mn = mean(input_list)
    variance = sum([(e-mn)**2 for e in input_list])
    return math.sqrt(variance)


def stderr(input_list):
    """Calculates the standard error for an input list."""
    mn = mean(input_list)
    variance = sum([(e-mn)**2 for e in input_list])
    return math.sqrt(variance) / math.sqrt(len(input_list))
