"""The Python interface to FINUFFT is divided into two parts: the simple
interface (through the ``nufft*`` functions) and the more advanced plan
interface (through the ``Plan`` class). The former allows the user to perform
an NUFFT in a single call while the latter allows for more efficient reuse of
resources when the same NUFFT is applied several times to different data by
saving FFTW plans, sorting the nonuniform points, and so on.
"""


# start delvewheel patch
def _delvewheel_patch_1_10_0():
    import ctypes
    import os
    import platform
    import sys
    libs_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'finufft.libs'))
    is_conda_cpython = platform.python_implementation() == 'CPython' and (hasattr(ctypes.pythonapi, 'Anaconda_GetVersion') or 'packaged by conda-forge' in sys.version)
    if sys.version_info[:2] >= (3, 8) and not is_conda_cpython or sys.version_info[:2] >= (3, 10):
        if os.path.isdir(libs_dir):
            os.add_dll_directory(libs_dir)
    else:
        load_order_filepath = os.path.join(libs_dir, '.load-order-finufft-2.5.1')
        if os.path.isfile(load_order_filepath):
            import ctypes.wintypes
            with open(os.path.join(libs_dir, '.load-order-finufft-2.5.1')) as file:
                load_order = file.read().split()
            kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)
            kernel32.LoadLibraryExW.restype = ctypes.wintypes.HMODULE
            kernel32.LoadLibraryExW.argtypes = ctypes.wintypes.LPCWSTR, ctypes.wintypes.HANDLE, ctypes.wintypes.DWORD
            for lib in load_order:
                lib_path = os.path.join(os.path.join(libs_dir, lib))
                if os.path.isfile(lib_path) and not kernel32.LoadLibraryExW(lib_path, None, 8):
                    raise OSError('Error loading {}; {}'.format(lib, ctypes.FormatError(ctypes.get_last_error())))


_delvewheel_patch_1_10_0()
del _delvewheel_patch_1_10_0
# end delvewheel patch

# that was the docstring for the package finufft.

__all__ = ["nufft1d1","nufft1d2","nufft1d3","nufft2d1","nufft2d2","nufft2d3","nufft3d1","nufft3d2","nufft3d3","Plan"]
# etc..

# let's just get guru and nufft1d1 working first...
from finufft._interfaces import Plan
from finufft._interfaces import nufft1d1,nufft1d2,nufft1d3
from finufft._interfaces import nufft2d1,nufft2d2,nufft2d3
from finufft._interfaces import nufft3d1,nufft3d2,nufft3d3

__version__ = '2.5.1'