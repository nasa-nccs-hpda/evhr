import concurrent.futures
import functools
import time

import ipywidgets as widgets


def progress_bar(float_progress: widgets.widget_float.FloatProgress, 
                 expected_time: int, 
                 increments=100):

    def _progress_bar(func):

        def timed_progress_bar(future, 
                               float_progress, 
                               expected_time, 
                               increments=100):
            """
            Display progress bar for expected_time seconds.
            Complete early if future completes.
            Wait for future if it doesn't complete in expected_time.
            """
            interval = expected_time / increments
            for i in range(increments - 1):
                if future.done():
                    # finish the progress bar
                    # not sure if there's a cleaner way to do this?
                    float_progress.value = expected_time
                    float_progress.bar_style = 'success'
                    return
                else:
                    time.sleep(interval)
                    float_progress.value = float_progress.value + interval
            # if the future still hasn't completed, wait for it.
            future.result()
            float_progress.value = expected_time
            float_progress.bar_style = 'success'

        @functools.wraps(func)
        def _func(*args, **kwargs):
            with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
                future = pool.submit(func, *args, **kwargs)
                timed_progress_bar(future,
                                   float_progress,
                                   expected_time,
                                   increments)

            return future.result()

        return _func

    return _progress_bar
