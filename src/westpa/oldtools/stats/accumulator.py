import numpy

NAN = float('nan')


class RunningStatsAccumulator:
    def __init__(self, shape, dtype=numpy.float64, count_dtype=numpy.uint, weight_dtype=numpy.float64, mask_value=NAN):
        self.sum = numpy.zeros(shape, dtype)
        self.sqsum = numpy.zeros(shape, dtype)
        self.weight = numpy.zeros(shape, weight_dtype)
        self.count = numpy.zeros(shape, count_dtype)
        self.mask_value = mask_value

    def incorporate(self, index, value, weight):
        self.count[index] += 1
        self.weight[index] += weight
        self.sum[index] += weight * value
        self.sqsum[index] += weight * value * value

    def average(self):
        valid = self.count > 0
        avg = numpy.empty_like(self.sum)
        avg[valid] = self.sum[valid] / self.weight[valid]
        avg[~valid] = self.mask_value
        return avg

    mean = average

    def std(self):
        valid = self.count > 0
        vavg = self.average()[valid]
        std = numpy.empty_like(self.sqsum)
        std[valid] = (self.sqsum[valid] / self.weight[valid] - vavg * vavg) ** 0.5
        std[~valid] = self.mask_value
        return std
