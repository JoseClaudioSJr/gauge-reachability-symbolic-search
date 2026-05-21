"""
Beale Conjecture: Representation-Dependent Obstruction Analysis

Search for A^x + B^y = C^z solutions (with z fixed to 5 for tractability)
and analyze how different representations reveal algebraic obstructions.
"""

import json
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import statistics
from math import gcd
from functools import reduce

class BealeRepresentation:
    """Base class for Beale representations."""
    
    def __init__(self, z: int = 5, a_max: int = 200, b_max: int = 200, 
                 c_max: int = 200, xy_range: Tuple = (3, 5)):
        self.z = z
        self.a_max = a_max
        self.b_max = b_max
        self.c_max = c_max
        self.xy_min, self.xy_max = xy_range
        self.solutions = []
        self.near_misses = []
        self.obstructions = {}
    
    def gcd_triple(self, a, b, c):
        """Compute GCD of triple."""
        return reduce(gcd, [a, b, c])
    
    def check_solution(self, a, b, c, x, y):
        """Check if A^x + B^y = C^z."""
        lhs = a**x + b**y
        rhs = c**self.z
        return lhs == rhs
    
    def find_near_misses(self, tolerance: int = 100):
        """Find (A,B,C,x,y) where A^x + B^y is close to C^z."""
        near = []
        for c in range(2, min(self.c_max, 50)):  # smaller c_max for near misses
            for z_val in [self.z]:
                target = c**z_val
                for x in range(self.xy_min, self.xy_max + 1):
                    for a in range(2, min(self.a_max, 50)):
                        val_a = a**x
                        if val_a > target + tolerance:
                            break
                        for y in range(self.xy_min, self.xy_max + 1):
                            for b in range(2, min(self.b_max, 50)):
                                val_b = b**y
                                lhs = val_a + val_b
                                if abs(lhs - target) <= tolerance and lhs != target:
                                    near.append({
                                        'a': a, 'b': b, 'c': c,
                                        'x': x, 'y': y, 'z': self.z,
                                        'lhs': lhs, 'rhs': target,
                                        'error': lhs - target,
                                        'gcd': self.gcd_triple(a, b, c)
                                    })
        return near


class DirectSearchRepresentation(BealeRepresentation):
    """Brute-force search: iterate all (A,B,C,x,y), check solutions."""
    
    def analyze(self):
        """Search for exact solutions."""
        solutions_found = 0
        coprime_only = 0
        
        for c in range(2, self.c_max):
            target = c**self.z
            if target > 10**10:  # limit search space
                break
            
            for x in range(self.xy_min, self.xy_max + 1):
                for a in range(2, self.a_max):
                    val_a = a**x
                    if val_a > target:
                        break
                    
                    for y in range(self.xy_min, self.xy_max + 1):
                        for b in range(2, self.b_max):
                            val_b = b**y
                            if val_a + val_b == target:
                                g = self.gcd_triple(a, b, c)
                                solutions_found += 1
                                self.solutions.append({
                                    'a': a, 'b': b, 'c': c,
                                    'x': x, 'y': y, 'z': self.z,
                                    'gcd': g
                                })
                                if g == 1:
                                    coprime_only += 1
                            elif val_a + val_b > target:
                                break
        
        self.near_misses = self.find_near_misses(tolerance=1000)
        
        return {
            'name': 'direct',
            'solutions_found': solutions_found,
            'coprime_solutions': coprime_only,
            'near_misses': len(self.near_misses),
            'beale_violation': coprime_only,  # if >0, Beale is false
            'patterns': 'none - exhaustive enumeration'
        }


class ModularObstructionRepresentation(BealeRepresentation):
    """Modular arithmetic: find obstructions mod p."""
    
    def analyze(self):
        """Analyze which residue classes are forbidden."""
        primes = [2, 3, 5, 7, 11, 13]
        forbidden_triples = defaultdict(int)
        allowed_triples = defaultdict(int)
        
        for p in primes:
            for a_mod in range(p):
                for b_mod in range(p):
                    for c_mod in range(p):
                        # Check if there exist x,y such that a^x + b^y ≡ c^z (mod p)
                        possible = False
                        
                        for x in range(self.xy_min, self.xy_max + 1):
                            for y in range(self.xy_min, self.xy_max + 1):
                                a_pow = pow(a_mod, x, p) if a_mod > 0 else 0
                                b_pow = pow(b_mod, y, p) if b_mod > 0 else 0
                                c_pow = pow(c_mod, self.z, p) if c_mod > 0 else 0
                                
                                if (a_pow + b_pow) % p == c_pow:
                                    possible = True
                                    break
                            if possible:
                                break
                        
                        key = f'mod_{p}_({a_mod},{b_mod},{c_mod})'
                        if possible:
                            allowed_triples[p] += 1
                        else:
                            forbidden_triples[p] += 1
        
        # Check near-misses for modular violations
        violations_by_prime = defaultdict(int)
        for miss in self.near_misses:
            for p in primes:
                a_mod = miss['a'] % p
                b_mod = miss['b'] % p
                c_mod = miss['c'] % p
                x, y, z = miss['x'], miss['y'], miss['z']
                
                lhs_mod = (pow(a_mod, x, p) + pow(b_mod, y, p)) % p
                rhs_mod = pow(c_mod, z, p) % p
                
                if lhs_mod != rhs_mod:
                    violations_by_prime[p] += 1
        
        return {
            'name': 'modular',
            'primes_checked': primes,
            'forbidden_classes': dict(forbidden_triples),
            'allowed_classes': dict(allowed_triples),
            'modular_violations_in_near_misses': dict(violations_by_prime),
            'obstructions': 'residue classes that cannot satisfy equation'
        }


class GCDStructureRepresentation(BealeRepresentation):
    """GCD analysis: study common factor structure in near-misses."""
    
    def analyze(self):
        """Analyze GCD patterns in solutions and near-misses."""
        gcd_distribution = defaultdict(int)
        gcd_by_error = defaultdict(list)
        
        self.near_misses = self.find_near_misses(tolerance=10000)
        
        for miss in self.near_misses:
            g = miss['gcd']
            error = miss['error']
            gcd_distribution[g] += 1
            gcd_by_error[g].append(error)
        
        # Analyze if coprime triples are farther from solutions
        coprime_errors = [m['error'] for m in self.near_misses if m['gcd'] == 1]
        nontrivial_errors = [m['error'] for m in self.near_misses if m['gcd'] > 1]
        
        stats = {
            'name': 'gcd',
            'near_misses': len(self.near_misses),
            'gcd_distribution': dict(gcd_distribution),
            'coprime_near_misses': len(coprime_errors),
            'nontrivial_gcd_near_misses': len(nontrivial_errors),
            'coprime_mean_error': statistics.mean(coprime_errors) if coprime_errors else None,
            'nontrivial_mean_error': statistics.mean(nontrivial_errors) if nontrivial_errors else None,
            'beale_evidence': 'if coprime errors >> nontrivial, supports Beale'
        }
        
        return stats


class AlgebraicGeometryRepresentation(BealeRepresentation):
    """Algebraic: analyze Fermat curve structure."""
    
    def analyze(self):
        """Study x^p + y^q = z^r as algebraic curve."""
        # For fixed z, analyze the curve x^a + y^b = z^5 (z fixed)
        # Count rational points, analyze genus, etc.
        
        curve_properties = {}
        
        for a in range(self.xy_min, self.xy_max + 1):
            for b in range(self.xy_min, self.xy_max + 1):
                # Fermat-like curve: x^a + y^b = 1 (projective)
                # Genus of Fermat curve x^a + y^b + z^c = 0 is
                # g = (a-1)(b-1)(c-1)/2 when gcd(a,b,c)=1
                
                from math import gcd as math_gcd
                g_abc = math_gcd(math_gcd(a, b), self.z)
                
                # For Fermat curves with genus > 1, Faltings says finitely many points
                genus = (a - 1) * (b - 1) * (self.z - 1) // 2 if g_abc == 1 else "singular"
                
                key = f'x^{a}_y^{b}_z^{self.z}'
                curve_properties[key] = {
                    'genus': genus,
                    'expected_points': 'finite' if isinstance(genus, int) and genus > 1 else 'unknown'
                }
        
        # Count actual rational points with small coordinates
        rational_points = 0
        for c in range(2, min(self.c_max, 100)):
            for x in range(self.xy_min, self.xy_max + 1):
                for a in range(2, min(self.a_max, 50)):
                    for y in range(self.xy_min, self.xy_max + 1):
                        for b in range(2, min(self.b_max, 50)):
                            if a**x + b**y == c**self.z:
                                rational_points += 1
        
        return {
            'name': 'algebraic',
            'curves_analyzed': len(curve_properties),
            'high_genus_curves': len([v for v in curve_properties.values() 
                                     if isinstance(v['genus'], int) and v['genus'] > 1]),
            'rational_points_found': rational_points,
            'faltings_prediction': 'high genus ⟹ finitely many rational points',
            'obstructions': 'geometric obstruction to rational solutions'
        }


def compare_representations():
    """Run all four representations and compare."""
    
    results = {
        'z': 5,
        'a_max': 200,
        'b_max': 200,
        'c_max': 200,
        'xy_range': [3, 5],
        'representations': {}
    }
    
    reps = [
        DirectSearchRepresentation(),
        ModularObstructionRepresentation(),
        GCDStructureRepresentation(),
        AlgebraicGeometryRepresentation()
    ]
    
    for rep in reps:
        print(f"Analyzing {rep.analyze()['name']}...")
        analysis = rep.analyze()
        results['representations'][analysis['name']] = analysis
    
    return results


if __name__ == '__main__':
    results = compare_representations()
    
    # Print results
    print("\n" + "="*70)
    print("BEALE CONJECTURE: Representation-Dependent Obstruction Analysis")
    print("="*70)
    
    for name, analysis in results['representations'].items():
        print(f"\n--- {name.upper()} REPRESENTATION ---")
        for key, value in analysis.items():
            if key != 'name':
                print(f"  {key}: {value}")
    
    # Save to JSON
    with open('beale_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n✓ Results saved to beale_results.json")
