# Python modules
from __future__ import division
import json
import sys

# 3rd party modules
import numpy as np
#import matplotlib.pyplot as pyplot

# This project's modules
import hlsvd


def main(data_filename=None):
    with open(data_filename) as f:
        input_ = json.load(f)

    # Signals (input data) are stored as correlated arrays of real/imaginary
    # pairs. Here I knit them into a single complex list.
    input_['signals'] = [complex(r, i) for r, i in zip(input_['signals_r'],
                                                       input_['signals_i'])]
    del input_['signals_r']
    del input_['signals_i']

    observed = input_["signals"]
    step_size = input_["step_size"]

    # Call hlsvd()
    results = hlsvd.hlsvd(observed, input_["n_singular_values"], step_size)

    nsv_found, singular_values, frequencies, \
        damping_factors, amplitudes, phases = results

    estimated = _create_hlsvd_sum_fid(frequencies, damping_factors, amplitudes,
                                      phases, len(observed), step_size)
    #print("\nnsv_found = \n")
    #print nsv_found,
    #print("\nfrequencies = \n")
    #print frequencies,
    #print("\namplitudes = \n")
    #print amplitudes,
    #print("\n phases = \n")
    #print phases,
    #print("\n damping_factors = \n")
    #print damping_factors,
    data={"damping_factors":damping_factors, "phases":phases, "amplitudes":amplitudes,"frequencies":frequencies, "nsv_found":nsv_found}
    data_filename+='hlsvd'
    with open (data_filename,'w') as outfile:
        json.dump(data,outfile)
    # We're only concerned with the real portion of the data
    observed = np.real(np.array(observed))
    estimated = np.real(estimated)

    rmse = np.sqrt(np.sum((observed - estimated) ** 2) / len(observed))

    nrmse = rmse / max(observed) - min(observed)

    print("RMSE = %.2f, normalized RMSE = %.2f%%" % (rmse, nrmse / 100))

    # The observed values are always noisier, so plotting them first allows
    # the cleaner estimation to be displayed on top.
    #pyplot.plot(observed, color='r')
    #pyplot.plot(estimated, color='b')
    #pyplot.show()


# This function written by Brian Soher for Vespa:
# https://scion.duhs.duke.edu/vespa/
# Modified by me (Philip Semanchuk); any errors are mine.
def _create_hlsvd_sum_fid(frequencies, decays, areas, phases,
                          ndata_points, dwell_time):
    """Construct time domain signal from the estimated parameters."""
    result = np.zeros((len(frequencies), ndata_points), dtype=complex)
    t = np.arange(ndata_points) * dwell_time

    # K is an invariant for the loop below
    K = 1j * 2 * np.pi

    for i, decay in enumerate(decays):
        if decays[i]:
            # We frequently force the exp() function here into very small
            # values that raise the following exception:
            #    FloatingPointError: underflow encountered in exp
            # Handling those exceptions with try/except is not an option
            # because raising the exception interrupts the calculation.
            # Instead, we temporarily disable numpy's error reporting,
            # allow it to silently generate NaNs, and then change those NaNs
            # to 0.0.
            old_settings = np.seterr(all='ignore')
            print areas[i],

            line_result = areas[i] * np.exp((t/decays[i]) +
                          K * (frequencies[i]*t+phases[i]/360.0))

            zeros = np.zeros_like(line_result)

            result[i,:] = np.where(np.isnan(line_result), zeros, line_result)

            np.seterr(**old_settings)
        else:
            result[i,:] = result[i,:] * 0
    #np.savetxt('uzay_test.txt',t.view(float))
    return np.sum(result, axis=0)

if __name__ == "__main__":
    data_filename = sys.argv[1] if (len(sys.argv) > 1) else 'uzay1.json'
    main(data_filename)
