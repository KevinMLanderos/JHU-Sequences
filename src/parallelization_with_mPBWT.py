import numpy as np
from concurrent.futures import ThreadPoolExecutor
#from mPBWT import iterate_over_positions2

######################################mPBWT
def find_long_matches2(xk, t, ak, dk, k, L):
    """
    Modified find_long_matches to report haplotype matches with start and end positions.
    Parameters:
        xk (list): Contains alleles at position k for all haplotypes.
        t (int): Number of possible alleles (e.g., t-allelic sites).
        ak (list): PBWT ordering of indices at position k.
        dk (list): Divergence array at position k.
        k (int): Current position index.
        L (int): Minimum match length to consider.
    Outputs:
        Reports matches between haplotypes, including start and end positions.
    """
    # Initialize variables
    m = [False] * t  # Tracks the presence of alleles
    i0 = 0           # Start index of the current interval
    M = len(dk)      # Total number of haplotypes
    # Iterate through the haplotypes in PBWT order
    for i in range(M):
        # Check if divergence exceeds the threshold
        if dk[i] > k - L:
            # Identify distinct alleles in the current interval
            distinct_alleles = [j for j in range(t) if m[j]]
            if len(distinct_alleles) >= 2:  # At least two distinct alleles
                # Process matches in the interval
                for ia in range(i0, i):
                    dmin = 0
                    for ib in range(ia + 1, i):
                        if dk[ib] > dmin:
                            dmin = dk[ib]
                        if xk[ak[ia]] == xk[ak[ib]]:  # If alleles match
                            print(f"Match: Haplotype {ak[ia]} and Haplotype {ak[ib]}, Start: {dmin}, End: {k}")
                        else:
                            if k - dmin >= L:  # Ensure the match length meets the threshold
                                print(f"Match: Haplotype {ak[ia]} and Haplotype {ak[ib]}, Start: {dmin}, End: {k - 1}")
            # Reset the start index and allele presence tracking
            i0 = i
            m = [False] * t
        # Update allele presence for the current haplotype
        allele = xk[ak[i]]
        m[allele] = True
    # Final check for the remaining interval
    distinct_alleles = [j for j in range(t) if m[j]]
    if len(distinct_alleles) >= 2:  # At least two distinct alleles
        for ia in range(i0, M):
            dmin = 0
            for ib in range(ia + 1, M):
                if dk[ib] > dmin:
                    dmin = dk[ib]
                if xk[ak[ia]] == xk[ak[ib]]:  # If alleles match
                    print(f"Match: Haplotype {ak[ia]} and {ak[ib]}, Start: {dmin}, End: {k}")
                else:
                    if k - dmin >= L:  # Ensure the match length meets the threshold
                        print(f"Match: Haplotype {ak[ia]} and {ak[ib]}, Start: {dmin}, End: {k - 1}")


# Algorithm 2: Update prefix and divergence arrays
def prefix_and_divergence_org(xk, t, ak, dk, k): # does not consider an input s
    """
    Implements the PREFIXANDDIVERGENCE procedure.
    Parameters:
        xk (list): List of t-allelic variant sites at position k.
        t (int): Number of alleles.
        ak (list): PBWT ordering of indices at position k.
        dk (list): Divergence array at position k.
        k (int): Current position index.
    Returns:
        ak_next (list): Updated PBWT ordering for position k+1.
        dk_next (list): Updated divergence array for position k+1.
    """
    # Initialize arrays for storing a[t], d[t], p[t], and u[t]
    a = [[] for _ in range(t)]  # a[t] will store arrays for each allele
    d = [[] for _ in range(t)]  # d[t] will store divergence arrays for each allele
    p = [k + 1] * t  # Initialize p[j] to k + 1 for all j
    u = [0] * t      # Initialize u[j] to 0 for all j
    M = len(ak)  # Number of haplotypes
    # Loop through all haplotypes
    for i in range(M):
        allele = xk[ak[i]]  # Identify the allele for the current haplotype
        # Append the current haplotype index to the appropriate a[allele]
        a[allele].append(ak[i])
        # Update p[j] values for all alleles
        for j in range(t):
            if dk[i] > p[j]:
                p[j] = dk[i]
        # Store the updated p[allele] in d[allele]
        d[allele].append(p[allele])
        # Reset p[allele] for the current allele
        p[allele] = 0
        # Increment the counter for the current allele
        u[allele] += 1
    # Concatenate arrays a[i] and d[i] for all alleles to form ak+1 and dk+1
    ak_next = [item for sublist in a for item in sublist]
    dk_next = [item for sublist in d for item in sublist]
    return ak_next, dk_next

def iterate_over_positions2(s, X, L, t):
    """
    Combines PREFIXANDDIVERGENCE, FINDLONGMATCHES, FINDLONGMATCHESLASTK, 
    and SETMAXIMALMATCHESLASTK to iterate over all positions k.
    Parameters:
        s (list): New sequence to match.
        X (list of lists): Reference panel of haplotypes (list of sequences).
        L (int): Minimum match length.
        t (int): Number of alleles.
    Outputs:
        Matches for each position k.
    """
    # Include new sequence to  X
    X = X + [s] 
    # Initialize PBWT ordering and divergence arrays
    ak = list(range(len(X)))  # Initial PBWT order
    dk = [0] * len(X)         # Initial divergence array
    M = len(X)                # Number of sequences in the reference panel
    N = len(s)                # Length of the sequences
    # Iterate over positions
    for k in range(N):
        print(f"\nPosition {k}")
        # Extract xk for the current position from the reference panel (used to get k+1 order)
        xk_current = [X[m][k] for m in range(len(X))] # NOTE: NEEDS TO BE ADJUSTED TO THE NEW ORDER
        # xk_current = [X[hap][k] for hap in ak]
        print(f"xk_current: {xk_current}")
        if k!=0:
            # Extract xk for the previous position from the reference panel (used to get matches)
            # xk_prev = [X[m][k-1] for m in range(len(X))] # NOTE: NEEDS TO BE ADJUSTED TO THE NEW ORDER
            xk_prev = [X[hap][k-1] for hap in ak]
            print(f"xk_prev: {xk_prev}")
        # 
        if k!=0:
            # Step 1: Find long matches for the current position
            print("Finding matches...")
            find_long_matches2(xk_prev, t, ak, dk, k, L)
        # Step 2: Obtain Prefix and divergence arrays for k+1
        print("Calculating prefix and divergence array for next position...")
        ak, dk = prefix_and_divergence_org(xk_current, t, ak, dk, k)  # NOTE: CHECK
        print(ak, dk)
    return ak, dk

####################################Parallelization 

def split_arrays(X, num_threads):
    """
    Combines split X/s based on number of threads 
    Parameters:
        X (list of lists): Reference panel of haplotypes (list of sequences).
        num_threads (int): number of threads 
    Outputs:
        list of lists 
    """
    n_cols = len(X[0])  # num of columns
    base_width = n_cols // num_threads  # base number of columns per part
    remainder = n_cols % num_threads    # extra columns to distribute
    X_parts = []
    start = 0
    for i in range(num_threads):
        width = base_width + (1 if i < remainder else 0)
        end = start + width
        X_parts.append([row[start:end] for row in X])
        start = end
    return X_parts

def process_part(s_part, X_part, L, t):
    """
    run mPBWT on X in multiple threads 
    Parameters:
        X (list of lists): part of X in one thread 
        s(list): part of s in one thread 
        L (int): Minimum match length.
        t (int): Number of alleles.
    Outputs:
        list of lists 
    """
    ak2, dk2 = iterate_over_positions2(s_part, X_part, L, t)
    return (ak2, dk2)

def reorder_duplicates(duplicates, current_ak, previous_ak):
    """
    Generate new_dk by reordering Xs and counting matches from right to left
    
    Args:
        Xs (list): List of binary lists representing the patterns
        corrected_ak (list): Corrected ordering to apply to Xs
        
    Returns:
        list: new_dk values
    """
    result = current_ak.copy()
    
    for duplicate_values in duplicates.values():
        # Get current positions of duplicate values in current_ak
        current_positions = [result.index(val) for val in duplicate_values]
        
        # Get positions in previous_ak and create ordering
        prev_positions = [previous_ak.index(val) for val in duplicate_values]
        # Create pairs of (value, position in previous_ak)
        value_prev_pos = list(zip(duplicate_values, prev_positions))
        # Sort by position in previous_ak
        value_prev_pos.sort(key=lambda x: x[1])
        # Extract just the values in correct order
        ordered_values = [x[0] for x in value_prev_pos]
        
        # Sort current positions
        current_positions.sort()
        
        # Create mapping of where each value should go
        for i, pos in enumerate(current_positions):
            result[pos] = ordered_values[i]
    
    return result
def generate_dk(Xs, corrected_ak):
    """
    Generate new_dk by reordering Xs and counting matches from right to left
    Args:
        Xs (list): List of binary lists representing the patterns
        corrected_ak (list): Corrected ordering to apply to Xs
    Returns:
        list: new_dk values
    """
    # Get length of binary patterns
    n = len(Xs[0])
    
    # Create mapping from old position to new position
    # Ensure we only map valid positions
    position_mapping = {val: new_pos for new_pos, val in enumerate(corrected_ak)}
    
    # Reorder Xs based on corrected_ak
    reordered_Xs = [None] * len(Xs)
    for old_pos, row in enumerate(Xs):
        # Check if the position exists in our mapping
        if old_pos >= len(corrected_ak):
            continue
        new_pos = position_mapping[old_pos]
        reordered_Xs[new_pos] = row
    
    # Calculate dk values by counting matches from right to left
    new_dk = []
    for i in range(len(reordered_Xs)):
        if i == 0:  # First row has no row to compare with
            new_dk.append(n)
            continue
            
        # Skip None values that might exist due to mapping issues
        if reordered_Xs[i] is None or reordered_Xs[i-1] is None:
            new_dk.append(n)
            continue
            
        # Compare current row with previous row from right to left
        current_row = reordered_Xs[i]
        previous_row = reordered_Xs[i-1]
        
        # Count matches from right until mismatch
        matches = 0
        for j in range(n-1, -1, -1):  # Start from rightmost position
            if current_row[j] == previous_row[j]:
                matches += 1
            else:
                break
                
        # dk value is total length minus matches
        new_dk.append(n - matches)
    
    return new_dk

def generate_ak_dk_parallel(X,s,num_threads,t,L):
    """
    run mPBWT in parallel 
        
    Returns:
        list: ak and dk from each thread 
    """
    X_parts = split_arrays(X, num_threads)
    s_parts = split_arrays([s],num_threads)
    s_parts = [inner_list[0] for inner_list in s_parts]
    results = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [
            executor.submit(process_part, s_part, X_part, L, t)
            for s_part, X_part in zip(s_parts, X_parts)
        ]
        for future in futures:
            results.append(future.result())

    Xs = X+[s]
    Xs_parts = split_arrays(Xs,num_threads)
    Xs_parts = [
        [''.join(str(num) for num in row) for row in part]
        for part in X_parts
    ]

    for idx in range(1,len(results)):
        chunk_str = Xs_parts[idx]
        duplicates = {x: [i for i, v in enumerate(chunk_str) if v == x] 
                for x in set(chunk_str) if chunk_str.count(x) > 1}
        if duplicates:
            current_ak,current_dk=results[idx]
            pre_ak,pre_dk = results[idx-1]
            corrected_ak = reorder_duplicates(duplicates, current_ak, pre_ak)
            corrected_dk = generate_dk(Xs, corrected_ak)
            # rerun current_dk 
            results[idx] = (corrected_ak, corrected_dk)
    return results

#########################################Example
X = [[0, 1, 1, 1, 1, 0, 0, 0, 0, 0],
    [1, 1, 0, 1, 0, 0, 0, 1, 0, 1],
    [0, 1, 0, 1, 0, 0, 0, 1, 0, 1],
    [0, 0, 1, 1, 0, 1, 0, 1, 0, 0],
    [1, 0, 1, 1, 0, 0, 0, 0, 1, 0],
    [1, 0, 1, 1, 0, 0, 0, 0, 1, 1],
    [0, 1, 1, 1, 1, 0, 0, 0, 0, 1],
    [0, 1, 0, 0, 1, 0, 1, 0, 0, 1],
    [0, 1, 0, 0, 0, 0, 0, 0, 1, 0]]
s=[1, 0, 1, 1, 0, 0, 0, 0, 1, 0]
t = 2
L = 3
num_threads = 2
results = generate_ak_dk_parallel(X,s,num_threads,t,L)
results
#############################################Evaluation
import numpy as np
import time
import matplotlib.pyplot as plt
import seaborn as sns
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from typing import Dict, List

# Use your existing mPBWT functions
#from mPBWT import iterate_over_positions2

def generate_test_data(num_sequences, sequence_length, num_alleles, seed=None):
    """Generate test data with multiple alleles"""
    if seed is not None:
        np.random.seed(seed)
        
    # Generate reference panel
    X = np.random.randint(0, num_alleles, size=(num_sequences, sequence_length))
    
    # Generate target sequence as a mix of random and reference sequences
    s = []
    for j in range(sequence_length):
        if np.random.random() < 0.7:  # 70% chance to copy from X
            s.append(X[np.random.randint(0, num_sequences), j])
        else:
            s.append(np.random.randint(0, num_alleles))
            
    return X.tolist(), s

def run_single_test(params, num_runs=3):
    """Run a single test scenario with given parameters"""
    X, s = generate_test_data(
        params['num_sequences'],
        params['sequence_length'],
        params['num_alleles'],
        seed=42
    )
    
    # Run parallel version
    parallel_times = []
    for _ in range(num_runs):
        start_time = time.time()
        results_parallel = generate_ak_dk_parallel(X, s, params['num_threads'], 
                                                 params['num_alleles'], params['L'])
        parallel_times.append(time.time() - start_time)
    
    # Run sequential version
    sequential_times = []
    for _ in range(num_runs):
        start_time = time.time()
        results_sequential = iterate_over_positions2(s, X, params['L'], params['num_alleles'])
        sequential_times.append(time.time() - start_time)
    
    return {
        'parallel_times': parallel_times,
        'sequential_times': sequential_times,
        'avg_parallel': np.mean(parallel_times),
        'avg_sequential': np.mean(sequential_times),
        'avg_speedup': np.mean(sequential_times) / np.mean(parallel_times)
    }

# Define test scenarios
test_scenarios = [
    {
        'name': 'Small Binary',
        'params': {
            'num_sequences': 10,
            'sequence_length': 100,
            'num_alleles': 2,
            'num_threads': 4,
            'L': 3
        }
    },
    {
        'name': 'Medium Binary',
        'params': {
            'num_sequences': 50,
            'sequence_length': 200,
            'num_alleles': 2,
            'num_threads': 4,
            'L': 3
        }
    },
    {
        'name': 'Large Binary',
        'params': {
            'num_sequences': 100,
            'sequence_length': 500,
            'num_alleles': 2,
            'num_threads': 4,
            'L': 3
        }
    },
    {
        'name': 'Small Multi-allelic',
        'params': {
            'num_sequences': 10,
            'sequence_length': 100,
            'num_alleles': 100,
            'num_threads': 4,
            'L': 3
        }
    },
    {
        'name': 'Medium Multi-allelic',
        'params': {
            'num_sequences': 50,
            'sequence_length': 200,
            'num_alleles': 100,
            'num_threads': 4,
            'L': 3
        }
    },
    {
        'name': 'Large Multi-allelic',
        'params': {
            'num_sequences': 100,
            'sequence_length': 500,
            'num_alleles': 100,
            'num_threads': 4,
            'L': 3
        }
    }
]

# Run tests and collect results
all_results = []
for scenario in test_scenarios:
    print(f"\nRunning scenario: {scenario['name']}")
    print(f"Parameters: {scenario['params']}")
    
    results = run_single_test(scenario['params'])
    results['scenario'] = scenario['name']
    results['params'] = scenario['params']
    all_results.append(results)
    
    print(f"Average sequential time: {results['avg_sequential']:.4f} seconds")
    print(f"Average parallel time: {results['avg_parallel']:.4f} seconds")
    print(f"Average speedup: {results['avg_speedup']:.2f}x")

# Create visualization
plt.style.use('seaborn')
fig = plt.figure(figsize=(15, 10))
fig.suptitle('MPBWT Performance Analysis', fontsize=16)

# 1. Execution times comparison
ax1 = plt.subplot(221)
scenario_names = [r['scenario'] for r in all_results]
sequential_times = [r['avg_sequential'] for r in all_results]
parallel_times = [r['avg_parallel'] for r in all_results]

x = np.arange(len(scenario_names))
width = 0.35

ax1.bar(x - width/2, sequential_times, width, label='Sequential')
ax1.bar(x + width/2, parallel_times, width, label='Parallel')
ax1.set_ylabel('Time (seconds)')
ax1.set_title('Average Execution Times')
ax1.set_xticks(x)
ax1.set_xticklabels(scenario_names, rotation=45, ha='right')
ax1.legend()

# 2. Speedup comparison
ax2 = plt.subplot(222)
speedups = [r['avg_speedup'] for r in all_results]
bars = ax2.bar(scenario_names, speedups)
ax2.set_ylabel('Speedup Factor')
ax2.set_title('Parallel Speedup by Scenario')
ax2.set_xticklabels(scenario_names, rotation=45, ha='right')

for bar in bars:
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
            f'{height:.2f}x', ha='center', va='bottom')

# 3. Scaling with sequence length
ax3 = plt.subplot(223)
binary_scenarios = [(r['params']['sequence_length'], r['avg_speedup']) 
                   for r in all_results if r['params']['num_alleles'] == 2]
multi_scenarios = [(r['params']['sequence_length'], r['avg_speedup']) 
                  for r in all_results if r['params']['num_alleles'] == 100]

if binary_scenarios:
    lengths_binary, speedups_binary = zip(*binary_scenarios)
    ax3.plot(lengths_binary, speedups_binary, 'b-o', label='Binary')
if multi_scenarios:
    lengths_multi, speedups_multi = zip(*multi_scenarios)
    ax3.plot(lengths_multi, speedups_multi, 'r-o', label='Multi-allelic')

ax3.set_xlabel('Sequence Length')
ax3.set_ylabel('Speedup Factor')
ax3.set_title('Speedup vs Sequence Length')
ax3.legend()

# 4. Execution time distribution
ax4 = plt.subplot(224)
data_to_plot = []
labels = []
for r in all_results:
    data_to_plot.extend(r['sequential_times'])
    data_to_plot.extend(r['parallel_times'])
    labels.extend(['Sequential'] * len(r['sequential_times']))
    labels.extend(['Parallel'] * len(r['parallel_times']))

df = pd.DataFrame({'Time': data_to_plot, 'Type': labels})
sns.boxplot(x='Type', y='Time', data=df, ax=ax4)
ax4.set_title('Distribution of Execution Times')
ax4.set_ylabel('Time (seconds)')

plt.tight_layout()
plt.show()