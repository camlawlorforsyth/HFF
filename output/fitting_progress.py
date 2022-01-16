
import numpy as np

'''
# total bins fit after 7 days of running time
astro = np.array([242, 81, 175, 157, 0, 0])
quixote = np.array([0, 0, 0, 42, 267, 192])
total = astro + quixote
days = 7
per_day = int(np.round(np.sum(total)/days))
days_remaining = 15
fitted = np.sum(total) + per_day*days_remaining
full_sample = 3634
'''

'''
# total bins fit after 10 days of running time
astro = np.array([349, 113, 255, 222, 0, 0])
quixote = np.array([0, 0, 0, 47, 407, 268])
total = astro + quixote
days = 10
per_day = int(np.round(np.sum(total)/days))
days_remaining = 12
fitted = np.sum(total) + per_day*days_remaining
full_sample = 3634
'''

'''
# total bins fit after ~13.67 days of running time
astro = np.array([440, 171, 371, 299, 0, 0])
quixote = np.array([0, 0, 0, 47, 571, 363])
total = astro + quixote
print(np.sum(total))
days = 13.67
per_day = int(np.round(np.sum(total)/days))
print(per_day)
days_remaining = 8.5
fitted = int(np.sum(total) + per_day*days_remaining)
print(fitted)
full_sample = 3634
'''

# total bins fit after ~20.5 days of running time
astro = np.array([562, 239, 533, 469, 0, 0])
quixote = np.array([0, 0, 0, 47, 829, 552])
total = astro + quixote
print('Total so far: {}'.format(np.sum(total)))
days = 20.5
per_day = int(np.round(np.sum(total)/days))
print('Number of bins per day: {}'.format(per_day))
days_remaining = 2.75
fitted = int(np.sum(total) + per_day*days_remaining)
print(fitted)
full_sample = 3634

# had to kill m717_ID_1358_bin_7 as it was taking a very long time
# (>~ 7 days which is strange)
