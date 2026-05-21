"""
Visualization and interpretation of Beale representation analysis.

Key question: Which representation most clearly reveals Beale's obstruction?
- Direct: raw solution enumeration
- Modular: algebraic constraints via residue classes
- GCD: common factor structure  
- Algebraic: genus-based geometric obstruction
"""

import json
import sys

def interpret_results(results_file: str = 'beale_results.json'):
    """Load and interpret representation comparison."""
    
    try:
        with open(results_file) as f:
            results = json.load(f)
    except FileNotFoundError:
        print(f"Error: {results_file} not found. Run beale_representations.py first.")
        return
    
    print("\n" + "="*80)
    print("BEALE CONJECTURE: What Does Each Representation Reveal?")
    print("="*80)
    
    reps = results['representations']
    
    # Summary table
    print("\n| Representation | Key Finding |")
    print("|---|---|")
    
    for name in ['direct', 'modular', 'gcd', 'algebraic']:
        if name not in reps:
            continue
        rep = reps[name]
        
        if name == 'direct':
            coprime = rep.get('coprime_solutions', 0)
            print(f"| {name} | Coprime solutions: {coprime} (Beale violation: {'YES' if coprime > 0 else 'NO'}) |")
        elif name == 'modular':
            forbidden = sum(rep.get('forbidden_classes', {}).values())
            print(f"| {name} | Forbidden residue classes: {forbidden} |")
        elif name == 'gcd':
            coprime_near = rep.get('coprime_near_misses', 0)
            print(f"| {name} | Coprime near-misses: {coprime_near} |")
        elif name == 'algebraic':
            rational_pts = rep.get('rational_points_found', 0)
            print(f"| {name} | Rational points found: {rational_pts} |")
    
    # Detailed interpretation
    print("\n" + "="*80)
    print("INTERPRETATION: What Structure Does Each Representation Reveal?")
    print("="*80)
    
    if 'direct' in reps:
        print("\n1. DIRECT REPRESENTATION (Brute-Force Search)")
        print("-" * 40)
        direct = reps['direct']
        print(f"   Visibility: LOW - binary yes/no for each triple")
        print(f"   Solutions found: {direct.get('solutions_found', 'N/A')}")
        print(f"   Coprime solutions: {direct.get('coprime_solutions', 'N/A')}")
        print(f"   Beale violation: {'FOUND!' if direct.get('coprime_solutions', 0) > 0 else 'None (consistent with Beale)'}")
        print(f"   Near-misses: {direct.get('near_misses', 'N/A')}")
        print(f"   Structure revealed: None (pure enumeration)")
    
    if 'modular' in reps:
        print("\n2. MODULAR REPRESENTATION (Residue Class Analysis)")
        print("-" * 40)
        modular = reps['modular']
        print(f"   Visibility: HIGH - algebraic constraints visible")
        print(f"   Primes checked: {modular.get('primes_checked', [])}")
        
        forbidden = modular.get('forbidden_classes', {})
        if forbidden:
            print(f"   Forbidden residue classes (by prime):")
            for p, count in forbidden.items():
                total = count + modular.get('allowed_classes', {}).get(p, 0)
                pct = 100 * count // total if total > 0 else 0
                print(f"      mod {p}: {count} forbidden out of {total} ({pct}%)")
        
        print(f"   Modular violations in near-misses: {modular.get('modular_violations_in_near_misses', {})}")
        print(f"   Structure revealed: Algebraic obstructions (which combinations impossible mod p)")
    
    if 'gcd' in reps:
        print("\n3. GCD STRUCTURE REPRESENTATION (Common Factor Analysis)")
        print("-" * 40)
        gcd = reps['gcd']
        print(f"   Visibility: MEDIUM-HIGH - common factor patterns")
        print(f"   Total near-misses: {gcd.get('near_misses', 'N/A')}")
        print(f"   Coprime near-misses: {gcd.get('coprime_near_misses', 'N/A')}")
        print(f"   Nontrivial GCD near-misses: {gcd.get('nontrivial_gcd_near_misses', 'N/A')}")
        
        print(f"   GCD distribution: {gcd.get('gcd_distribution', {})}")
        
        coprime_err = gcd.get('coprime_mean_error', None)
        nontrivial_err = gcd.get('nontrivial_mean_error', None)
        
        if coprime_err and nontrivial_err:
            print(f"   Mean error (coprime): {coprime_err:.0f}")
            print(f"   Mean error (nontrivial GCD): {nontrivial_err:.0f}")
            if coprime_err > nontrivial_err:
                print(f"   Evidence: Coprime errors are LARGER (supports Beale!)")
            else:
                print(f"   Evidence: Coprime errors are SMALLER (contradicts Beale)")
        
        print(f"   Structure revealed: Common factor structure in near-misses")
    
    if 'algebraic' in reps:
        print("\n4. ALGEBRAIC GEOMETRY REPRESENTATION (Fermat Curve Analysis)")
        print("-" * 40)
        alg = reps['algebraic']
        print(f"   Visibility: HIGHEST - geometric structure")
        print(f"   Curves analyzed: {alg.get('curves_analyzed', 'N/A')}")
        print(f"   High-genus curves: {alg.get('high_genus_curves', 'N/A')}")
        print(f"   Rational points found: {alg.get('rational_points_found', 'N/A')}")
        print(f"   Faltings prediction: {alg.get('faltings_prediction', 'N/A')}")
        print(f"   Structure revealed: Geometric obstruction (Faltings: high genus ⟹ finitely many points)")
    
    # Synthesis
    print("\n" + "="*80)
    print("SYNTHESIS: Representation Dependency in Beale")
    print("="*80)
    print("""
The same Beale equation A^x + B^y = C^z appears *completely different* in each representation:

- DIRECT:     "Are there coprime solutions?" (binary: yes/no, no insight)
- MODULAR:    "Which residue classes are forbidden?" (structure: algebraic obstructions)
- GCD:        "Do coprime triples stay farther from solutions?" (structure: common factor effect)
- ALGEBRAIC:  "What is the genus of the Fermat curves?" (structure: geometric obstruction)

KEY INSIGHT:
Beale's conjecture is NOT one question. It's FOUR different questions
depending on how you represent it. Each representation:
  - Makes different obstructions visible
  - Suggests different proof strategies
  - Reveals different levels of evidence FOR or AGAINST the conjecture

BEALE'S CLAIM:
"If gcd(A,B,C) = 1, then A^x + B^y ≠ C^z for x,y,z > 2"

In different representations this means:

1. DIRECT:     [Pure enumeration — no insight]

2. MODULAR:    "Coprime residue classes never satisfy x^a + y^b ≡ z^c (mod p)"
               → Suggests modular descent

3. GCD:        "Coprime triples are always farther from solutions than those with common factors"
               → Suggests valuation-based obstruction

4. ALGEBRAIC:  "The Fermat curve x^a + y^b = z^c has genus > 1 for most (a,b,c)"
               → By Faltings: finitely many rational points
               → Suggests Beale is geometric obstruction


REPRESENTATION ENABLES PROOF STRATEGY:
The right representation might contain the seed of the proof!
""")


if __name__ == '__main__':
    results_file = sys.argv[1] if len(sys.argv) > 1 else 'beale_results.json'
    interpret_results(results_file)
