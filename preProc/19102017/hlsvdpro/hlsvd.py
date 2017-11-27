"""This is a Python wrapper for the hlsvd function exposed by the HLSVDPRO library."""
# Python modules
from __future__ import division
import sys
import ctypes
import os


RADIANS_TO_DEGREES = 180 / 3.1415926

# These limits are defined in hlsvdpro.f as ndpmax and kmax.
MAX_DATA_POINTS = 8192
MAX_SINGULAR_VALUES = 50

# _LIBRARY_NAMES maps platforms to the name of the library on said platform.
_LIBRARY_NAMES = {"osx"    : "libhlsvdpro.dylib",
                  "windows": "hlsvdpro.dll",
                  "linux"  : "libhlsvdpro.so",
                  }

def hlsvd(signals, nsv_sought, dwell_time, libhlsvd=None):
    """
    Input --
    signals:
        An iterable (e.g. list) of complex numbers. A numpy array should
        work too. The list must have <= MAX_DATA_POINTS elements.

    nsv_sought:
        The number of singular values sought. The function will return a
        maximum of this many singular values.

        Must be <= MAX_SINGULAR_VALUES.

    dwell_time:
        Dwell time in milliseconds

    libhlsvd:
        None (the default) or a string.
        If None, this function will try to find the library on its own.
        If a string, it should be a file name (fully qualified or not,
        e.g. "hlsvd.dll" or "/home/philip/lib/libhlsvd.dylib") that
        this function can pass to ctypes.CDLL().

    Output is a 6-tuple of --
     0) nsv_found: the number of singular values found (<= nsv_sought)
     1) A list of floats representing the singular values
     2) A list of floats representing the frequencies (in kilohertz)
     3) A list of floats representing the damping factors (in milliseconds?)
     4) A list of floats representing the amplitudes (in arbitrary units)
     5) A list of floats representing the phases (in degrees)

    Each list's length == nsv_found. The five lists are correlated (element N
    of one list is associated with element N in the other four lists) and are
    sorted by singular value with the largest (strongest signal) first.

    The frequencies and damping factors have been adjusted by the dwell time.
    """
    libhlsvd = _load_the_library(libhlsvd)

    # At this point libhlsvd should be a valid reference to the hlsvd library.
    # I test that assertion here as best as I can by generating a throwaway
    # reference to the function I'm going to call.
    f = libhlsvd.hlsvdpw_python_

    # OK, all is well. I create and populate the variables I need as parameters
    # to the function.
    n_data_points = len(signals)
    print("RMSE = %.2f" % (n_data_points))
    InputDoubleArrayType = ctypes.c_double * MAX_DATA_POINTS
    OutputDoubleArrayType = ctypes.c_double * MAX_SINGULAR_VALUES

    # In-line comments refer to the names of the corresponding variables
    # in hlsvdpro.f.

    # Input params
    real_signals = InputDoubleArrayType()
    imaginary_signals = InputDoubleArrayType()
     
    for i, signal in enumerate(signals):
        real_signals[i] = signal.real                   # signal_r
        imaginary_signals[i] = signal.imag              # signal_i

  
    nsv_sought = ctypes.c_long(nsv_sought)              # kuser

    # I copied the forumula for calculating the size of the Hankel matrix
    # from the Fortran code.
    mcol = n_data_points // 2
    lrow = n_data_points - mcol + 1
    mcol = ctypes.c_long(mcol)                          # mcoL
    lrow = ctypes.c_long(lrow)                          # Lrow
   
    n_data_points = ctypes.c_long(n_data_points)        # ndp
   
    # Output params
    frequencies = OutputDoubleArrayType()               # freq
    damping_factors = OutputDoubleArrayType()           # damp
    amplitudes = OutputDoubleArrayType()                # ampl
    phases = OutputDoubleArrayType()                    # fase
    singular_values = OutputDoubleArrayType()           # Lsinval
    nsv_found = ctypes.c_long()                         # kfit

    libhlsvd.hlsvdpw_python_(real_signals,
                             imaginary_signals,
                             ctypes.pointer(n_data_points),
                             ctypes.pointer(lrow),
                             ctypes.pointer(mcol),
                             ctypes.pointer(nsv_sought),
                             ctypes.pointer(nsv_found),
                             ctypes.pointer(singular_values),
                             ctypes.pointer(amplitudes),
                             ctypes.pointer(phases),
                             ctypes.pointer(damping_factors),
                             ctypes.pointer(frequencies),
                             )

    # I tease the returned variables into Python types before passing them
    # back to the caller. (Slicing a ctypes array returns a list.) The Fortran
    # code has already sorted them the way we like (largest singular value
    # first).
    nsv_found       = nsv_found.value
    singular_values = singular_values[:nsv_found]
    frequencies     = frequencies[:nsv_found]
    damping_factors = damping_factors[:nsv_found]
    amplitudes      = amplitudes[:nsv_found]
    phases          = phases[:nsv_found]

    damping_factors = [1 / df for df in damping_factors]
    damping_factors = [df * dwell_time for df in damping_factors]

    frequencies = [frequency / dwell_time for frequency in frequencies]

    phases = [phase * RADIANS_TO_DEGREES for phase in phases]

    return (nsv_found, singular_values, frequencies, damping_factors, amplitudes, phases)

###################   Internal use functions   ################

def _load_the_library(libhlsvd=None):
    platform = sys.platform.lower()

    if "linux" in platform:
        platform = "linux"
    elif "darwin" in platform:
        platform = "osx"
    elif platform.startswith("win"):
        platform = "windows"
    else:
        platform = None

    if not libhlsvd:
        # Construct a name for the library.
        local_path = os.path.dirname(os.path.realpath(__file__))

        local_path = os.path.join(local_path, 'bin')

        libhlsvd = os.path.join(local_path, _LIBRARY_NAMES[platform])

    if platform == "windows":
        # Under Windows, setting the CWD helps Windows to find the HLSVDPRO DLL.
        # ref: http://msdn.microsoft.com/en-us/library/windows/desktop/ms682586(v=vs.85).aspx#search_order_for_desktop_applications
        restore_cwd = os.getcwd()
        libpath = os.path.dirname(libhlsvd)
        os.chdir(libpath)
    else:
        restore_cwd = None

    try:
        libhlsvd = ctypes.CDLL(libhlsvd)
    except OSError as exception_instance:
        raise OSError("Unable to load '%s'. The OS reports: %s" % (libhlsvd,
                                                                   str(exception_instance)))
    finally:
        if restore_cwd:
            os.chdir(restore_cwd)

    return libhlsvd
