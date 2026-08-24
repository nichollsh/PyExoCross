from __future__ import annotations
import os
import time
import platform

# Determine whether the folder exists or not.
def ensure_dir(path):
    """
    Ensure that a directory exists, creating it if necessary.

    If the directory already exists, no action is taken. Otherwise,
    the directory and any necessary parent directories are created.

    Parameters
    ----------
    path : str
        Path to the directory to ensure exists.
    """
    # If the folder exists, save files directory,otherwise, create it.
    if os.path.exists(path):
        pass
    else:
        # Create the folder.
        os.makedirs(path, exist_ok=True)
        
# Report time
class Timer:
    """
    Timer class for measuring CPU and system time intervals.
    
    This class tracks both CPU time (process time) and system time (wall clock time)
    for performance measurement purposes.
    """
    
    def start(self):
        """
        Start the timer.

        Records the current CPU time and system time as the start point.

        Returns
        -------
        Timer
            Returns self for method chaining.
        """
        self.start_CPU = time.process_time()
        self.start_sys = time.time()
        self.start_save_CPU = 0.0
        self.start_save_sys = 0.0
        from pyexocross import core
        if core.timing_active:
            save = core.timing_summary.get('save', {})
            self.start_save_CPU = save.get('cpu', 0.0)
            self.start_save_sys = save.get('system', 0.0)
        return self

    def end(self, *args):
        """
        End the timer and print the elapsed time intervals.

        Calculates and prints both CPU time and system time intervals
        since the timer was started.

        Parameters
        ----------
        *args
            Ignored, kept for backward compatibility.
        """
        self.end_CPU = time.process_time()
        self.end_sys = time.time()
        self.interval_CPU = self.end_CPU - self.start_CPU
        self.interval_sys = self.end_sys - self.start_sys
        if args and isinstance(args[0], str):
            from pyexocross.base.result import add_timing
            add_timing(args[0], self.interval_CPU, self.interval_sys)
        # print('{:25s} : {}'.format('Running time on CPU', self.interval_CPU), 's')
        # print('{:25s} : {}'.format('Running time on system', self.interval_sys), 's')
        
    def cal(self, *args):
        """
        Calculate elapsed time intervals without printing.

        Computes CPU time and system time intervals since the timer was started,
        but does not print them. Useful for programmatic access to timing data.

        Parameters
        ----------
        *args
            Ignored, kept for backward compatibility.

        Returns
        -------
        tuple of float
            A tuple containing (CPU_time_interval, system_time_interval) in seconds.
        """
        self.end_CPU = time.process_time()
        self.end_sys = time.time()
        self.interval_CPU = self.end_CPU - self.start_CPU
        self.interval_sys = self.end_sys - self.start_sys
        return(self.interval_CPU, self.interval_sys)

# Track process memory usage, e.g. to help size Slurm '--mem' requests.
def get_memory_usage_mb():
    """
    Return the process's current resident memory usage in MB.

    Best-effort: returns None if memory usage can't be determined on this
    platform, rather than raising, so callers can treat it as optional.

    Returns
    -------
    float or None
        Current resident set size (RSS) in MB, or peak RSS if current RSS
        is unavailable. None if neither could be read.
    """
    try:
        with open('/proc/self/status') as f:
            for line in f:
                if line.startswith('VmRSS:'):
                    return int(line.split()[1]) / 1024.0
    except OSError:
        pass
    try:
        import resource
        peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # Linux reports ru_maxrss in KB, macOS in bytes.
        return peak_kb / 1024.0 if platform.system() == 'Linux' else peak_kb / (1024.0 ** 2)
    except Exception:
        return None
