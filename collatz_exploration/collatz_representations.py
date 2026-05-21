"""
Collatz Conjecture: Representation-Dependent Structure

Hypothesis: Different representations of the Collatz problem reveal
different structural patterns and converge at different "visibility" rates.

Four representations:
1. Direct: naive iteration
2. Binary: work with bit patterns
3. Modular: analyze constraints mod 2^k
4. Graph: inverse tree structure
"""

import json
from collections import defaultdict, deque
from typing import Dict, List, Tuple, Set
import statistics

class CollatzRepresentation:
    """Base class for Collatz representations."""
    
    def __init__(self, n_max: int = 1000):
        self.n_max = n_max
        self.trajectories = {}
        self.patterns = {}
    
    def collatz_step(self, n: int) -> int:
        if n % 2 == 0:
            return n // 2
        else:
            return 3 * n + 1
    
    def trajectory_length(self, n: int, max_steps: int = 10000) -> int:
        """Count steps to reach 1."""
        seen = set()
        steps = 0
        while n != 1 and steps < max_steps:
            if n in seen:  # cycle detected
                return -1
            seen.add(n)
            n = self.collatz_step(n)
            steps += 1
        return steps if n == 1 else -1


class DirectRepresentation(CollatzRepresentation):
    """Naive iteration: just compute n -> 3n+1 or n/2."""
    
    def analyze(self):
        """Collect trajectory lengths."""
        lengths = []
        for n in range(2, self.n_max):
            length = self.trajectory_length(n)
            if length > 0:
                lengths.append(length)
                self.trajectories[n] = length
        
        return {
            'name': 'direct',
            'count': len(lengths),
            'mean_length': statistics.mean(lengths),
            'max_length': max(lengths),
            'min_length': min(lengths),
            'patterns': 'none - just raw iteration'
        }


class BinaryRepresentation(CollatzRepresentation):
    """Work with binary: track bit patterns, divisibility by 2^k."""
    
    def analyze(self):
        """Analyze bit patterns and divisibility structure."""
        binary_patterns = defaultdict(list)
        power_of_2_divisibility = defaultdict(int)
        
        for n in range(2, self.n_max):
            trajectory = []
            current = n
            steps = 0
            
            while current != 1 and steps < 1000:
                # Count trailing zeros (power of 2 divisibility)
                power = (current & -current).bit_length() - 1
                power_of_2_divisibility[power] += 1
                
                # Binary pattern of lowest 8 bits
                pattern = bin(current & 0xFF)[-8:]
                binary_patterns[pattern].append(current)
                
                current = self.collatz_step(current)
                steps += 1
                trajectory.append(current)
            
            self.trajectories[n] = steps
        
        return {
            'name': 'binary',
            'count': len(self.trajectories),
            'mean_steps': statistics.mean(self.trajectories.values()),
            'power_of_2_distribution': dict(sorted(power_of_2_divisibility.items())[:10]),
            'binary_patterns_found': len(binary_patterns),
            'most_common_patterns': sorted(
                [(k, len(v)) for k, v in binary_patterns.items()],
                key=lambda x: x[1],
                reverse=True
            )[:5]
        }


class ModularRepresentation(CollatzRepresentation):
    """Modular arithmetic: find invariants mod 2^k, mod 3, etc."""
    
    def analyze(self):
        """Look for modular constraints and cycles."""
        mod_cycles = defaultdict(set)
        mod_constraints = defaultdict(lambda: defaultdict(int))
        
        for n in range(2, self.n_max):
            current = n
            trajectory = []
            steps = 0
            
            while current != 1 and steps < 1000:
                # Track trajectory modulo various bases
                for modulus in [2, 4, 8, 3, 6, 9]:
                    mod_val = current % modulus
                    mod_constraints[modulus][mod_val] += 1
                    
                    if steps < 10:  # early steps
                        mod_cycles[f'{modulus}_{mod_val}'].add(current)
                
                current = self.collatz_step(current)
                steps += 1
                trajectory.append(current)
            
            self.trajectories[n] = steps
        
        return {
            'name': 'modular',
            'count': len(self.trajectories),
            'mean_steps': statistics.mean(self.trajectories.values()),
            'modular_invariants': {
                f'mod_{m}': len(v) for m, v in mod_constraints.items()
            },
            'early_cycle_classes': {
                k: len(v) for k, v in mod_cycles.items() if len(v) > 5
            }
        }


class GraphRepresentation(CollatzRepresentation):
    """Graph structure: inverse tree (which n lead to m?)."""
    
    def inverse_collatz(self, n: int) -> List[int]:
        """Find all m such that collatz(m) = n."""
        predecessors = []
        # If n = k/2, then m = 2k
        if True:
            predecessors.append(2 * n)
        # If n = 3m+1, then m = (n-1)/3
        if (n - 1) % 3 == 0 and (n - 1) // 3 > 0:
            m = (n - 1) // 3
            if m % 2 == 1:  # 3m+1 only for odd m
                predecessors.append(m)
        return predecessors
    
    def analyze(self):
        """Build inverse tree and analyze graph structure."""
        graph_depth = {}
        graph_branching = defaultdict(int)
        in_degree = defaultdict(int)
        
        # BFS from 1 backwards
        queue = deque([(1, 0)])  # (node, depth)
        visited = {1}
        
        while queue and len(visited) < self.n_max:
            node, depth = queue.popleft()
            graph_depth[node] = depth
            
            for pred in self.inverse_collatz(node):
                if pred not in visited and pred < 10 * self.n_max:
                    visited.add(pred)
                    queue.append((pred, depth + 1))
                    in_degree[node] += 1
                    graph_branching[depth] += 1
        
        # Forward trajectories for nodes in graph
        for n in visited:
            if n >= 2 and n < self.n_max:
                length = self.trajectory_length(n)
                if length > 0:
                    self.trajectories[n] = length
        
        return {
            'name': 'graph',
            'nodes_in_tree': len(visited),
            'max_depth': max(graph_depth.values()) if graph_depth else 0,
            'nodes_within_range': len([n for n in visited if 2 <= n < self.n_max]),
            'branching_by_depth': dict(sorted(graph_branching.items())[:20]),
            'mean_trajectory_length': statistics.mean(self.trajectories.values()) if self.trajectories else 0
        }


def compare_representations(n_max: int = 500):
    """Run all four representations and compare visibility."""
    
    results = {
        'n_max': n_max,
        'representations': {}
    }
    
    reps = [
        DirectRepresentation(n_max),
        BinaryRepresentation(n_max),
        ModularRepresentation(n_max),
        GraphRepresentation(n_max)
    ]
    
    for rep in reps:
        print(f"Analyzing {rep.analyze()['name']}...")
        analysis = rep.analyze()
        results['representations'][analysis['name']] = analysis
    
    return results


if __name__ == '__main__':
    import sys
    
    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    
    results = compare_representations(n_max)
    
    # Print results
    print("\n" + "="*70)
    print("COLLATZ: Representation-Dependent Structure Analysis")
    print("="*70)
    
    for name, analysis in results['representations'].items():
        print(f"\n--- {name.upper()} REPRESENTATION ---")
        for key, value in analysis.items():
            if key != 'name':
                print(f"  {key}: {value}")
    
    # Save to JSON
    with open('collatz_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n✓ Results saved to collatz_results.json")
