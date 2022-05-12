import sys
import numpy as np

if __name__ == "__main__":
    for line in sys.stdin:
        lower95, upper95 = np.array(line.split()).astype(np.float)

        mean = (upper95 + lower95) / 2
        sd = (upper95 - mean) / 1.959964

        print(*[round(x, 2) for x in np.random.normal(mean, sd, 10)], sep = "\t")
