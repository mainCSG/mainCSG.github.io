# Purpose

This file contains useful notes and code for developmental purposes.

# Code Profiling

> **NOTE:** Examples are pulled from Dot Array Class Usage.ipynb

Profile full script using `%prun`

Example:
    <pre><code>%prun dots.g_factors(voltage_configs)
    </pre></code> 

Profile a function line by line using `line_profiler`. After installing liner_profiler via `pip install line_profiler` run the following to load `line_profiler` and profile the `g_factors()` method.

Example: 
    <pre><code>%load_ext line_profiler
    %lprun -f dots.g_factors dots.g_factors(voltage_configs)
    </pre></code>

Profile a functions memory usage using `memory_profiler`. After installing liner_profiler via `pip install memory_profiler` run the following to load `memory_profiler` and profile the `g_factors()` method.

Example: <pre><code>%load_ext memory_profiler
%memit dots.g_factors(voltage_configs)</pre></code>