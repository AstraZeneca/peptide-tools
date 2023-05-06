import math


def mean(lst):
    """calculates mean"""
    return sum(lst) / len(lst)


def stddev(lst):
    """returns the standard deviation of lst"""
    mn = mean(lst)
    variance = sum([(e-mn)**2 for e in lst])
    return math.sqrt(variance)


def stderr(lst):
    """returns the standard error of the mean of lst"""
    mn = mean(lst)
    variance = sum([(e-mn)**2 for e in lst])
    return math.sqrt(variance) / math.sqrt(len(lst))
