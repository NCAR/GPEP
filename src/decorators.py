"""

Decorators

A decorator is a function that takes another function and extends the behavior of the latter function
without explicitly modifying it.

"""

# Package Import
import functools
import time

def timer(func):
    """Print the runtime of the decorated function"""
    @functools.wraps(func)
    def wrapper_timer(*args,**kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time-start_time
        print(f'Function {func.__name__!r} completed in {run_time:.4f} seconds')
        return value
    return wrapper_timer

