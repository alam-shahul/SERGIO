import numpy as np
from scipy.stats import t, norm

result = np.load('simple_clean.npy')
stds = result.std(axis=1, ddof=1)
upper = result.mean(axis=1) + t.ppf(0.975, 999) * stds / np.sqrt(1000)
lower = result.mean(axis=1) - t.ppf(0.975, 999) * stds / np.sqrt(1000)

print('AR')
print(lower)
print(upper)
print(result.mean(axis=1))

result = np.load('simple_clean_no_ar.npy')
stds = result.std(axis=1, ddof=1)
upper = result.mean(axis=1) + t.ppf(0.975, 999) * stds / np.sqrt(1000)
lower = result.mean(axis=1) - t.ppf(0.975, 999) * stds / np.sqrt(1000)

print('NO AR')
print(lower)
print(upper)
print(result.mean(axis=1))

result = np.load('simple_clean_pseudo.npy')
stds = result.std(axis=1, ddof=1)
upper = result.mean(axis=1) + t.ppf(0.975, 999) * stds / np.sqrt(1000)
lower = result.mean(axis=1) - t.ppf(0.975, 999) * stds / np.sqrt(1000)

print('PSEUDO AR')
print(lower[[0,2,4,5]])
print(upper[[0,2,4,5]])
print(result.mean(axis=1)[[0,2,4,5]])
