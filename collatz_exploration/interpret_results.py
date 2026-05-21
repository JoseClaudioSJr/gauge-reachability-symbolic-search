"""
Visualization and interpretation of Collatz representation analysis.

Key question: Which representation reveals the most structure?
- Direct: raw iteration efficiency
- Binary: bit-level patterns  
- Modular: algebraic constraints
- Graph: inverse tree structure
"""

import json
import sys

def interpret_results(results_file: str = 'collatz_results.json'):
    """Load and interpret representation comparison."""
    
    try:
        with open(results_file) as f:
            results = json.load(f)
    except FileNotFoundError:
        print(f"Error: {results_file} not found. Run collatz_representations.py first.")
        return
    
    print("\n" + "="*80)
    print("COLLATZ CONJECTURE: What Does Each Representation Reveal?")
    print("="*80)
    
    reps = results['representations']
    
    # Summary table
    print("\n| Representation | Metric | Value |")
    print("|---|---|---|")
    
    for name in ['direct', 'binary', 'modular', 'graph']:
        if name not in reps:
            continue
        rep = reps[name]
        
        if name == 'direct':
            print(f"| {name} | Mean steps to 1 | {rep.get('mean_length', 'N/A')} |")
        elif name == 'binary':
            print(f"| {name} | Binary patterns found | {rep.get('binary_patterns_found', 'N/A')} |")
        elif name == 'modular':
            print(f"| {name} | Modular invariants | {len(rep.get('modular_invariants', {}))} |")
        elif name == 'graph':
            print(f"| {name} | Max tree depth | {rep.get('max_depth', 'N/A')} |")
    
    # Detailed interpretation
    print("\n" + "="*80)
    print("INTERPRETATION: What Structure Does Each Representation Reveal?")
    print("="*80)
    
    if 'direct' in reps:
        print("\n1. DIRECT REPRESENTATION")
        print("-" * 40)
        direct = reps['direct']
        print(f"   Visibility: LOW - only raw trajectory lengths")
        print(f"   Mean steps to convergence: {direct.get('mean_length', 'N/A')}")
        print(f"   Max steps: {direct.get('max_length', 'N/A')}")
        print(f"   Structure revealed: None (black box iteration)")
    
    if 'binary' in reps:
        print("\n2. BINARY REPRESENTATION")
        print("-" * 40)
        binary = reps['binary']
        print(f"   Visibility: MEDIUM - bit patterns and divisibility")
        print(f"   Binary patterns found: {binary.get('binary_patterns_found', 'N/A')}")
        print(f"   Most common bit patterns:")
        for pattern, count in binary.get('most_common_patterns', []):
            print(f"      {pattern}: {count} occurrences")
        print(f"   Structure revealed: Power-of-2 divisibility structure")
    
    if 'modular' in reps:
        print("\n3. MODULAR REPRESENTATION")
        print("-" * 40)
        modular = reps['modular']
        print(f"   Visibility: HIGH - algebraic constraints visible")
        print(f"   Modular invariants checked: {list(modular.get('modular_invariants', {}).keys())}")
        print(f"   Early cycle classes found:")
        for cls, count in list(modular.get('early_cycle_classes', {}).items())[:5]:
            print(f"      {cls}: {count} numbers")
        print(f"   Structure revealed: Algebraic obstructions (mod p)")
    
    if 'graph' in reps:
        print("\n4. GRAPH REPRESENTATION")
        print("-" * 40)
        graph = reps['graph']
        print(f"   Visibility: HIGHEST - inverse tree structure")
        print(f"   Nodes in inverse tree: {graph.get('nodes_in_tree', 'N/A')}")
        print(f"   Maximum depth from 1: {graph.get('max_depth', 'N/A')}")
        print(f"   Branching structure by depth:")
        for depth, count in list(graph.get('branching_by_depth', {}).items())[:10]:
            bars = '█' * min(int(count) // 10 + 1, 50)
            print(f"      Depth {int(depth):2d}: {bars} ({count})")
        print(f"   Structure revealed: Tree topology, reachability from 1")
    
    # Synthesis
    print("\n" + "="*80)
    print("SYNTHESIS: Representation Dependency in Collatz")
    print("="*80)
    print("""
The same Collatz conjecture appears *completely different* in each representation:

- DIRECT:  "Does iteration reach 1?" (binary: yes/no, no insight)
- BINARY:  "What bit patterns emerge?" (structure: divisibility chains)
- MODULAR: "Which mod classes form cycles?" (structure: algebraic obstructions)
- GRAPH:   "What is the tree topology from 1?" (structure: reachability depth)

KEY INSIGHT:
The conjecture is NOT one problem. It's FOUR different problems
depending on how you represent it. Each representation:
  - Makes different aspects visible
  - Makes different algorithmic approaches plausible
  - Suggests different proof strategies

This confirms the principle: **representation determines what becomes visible**.
""")


if __name__ == '__main__':
    results_file = sys.argv[1] if len(sys.argv) > 1 else 'collatz_results.json'
    interpret_results(results_file)
