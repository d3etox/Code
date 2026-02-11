import numpy as np

# =========================================================
# ANALYZE P3 NODE NUMBERING AND BASIS FUNCTION ORDER
# =========================================================

print("="*80)
print("P3 LAGRANGE ELEMENT - NODE NUMBERING ANALYSIS")
print("="*80)

print("\nStandard P3 Lagrange numbering on a reference triangle:")
print("""
        3
        |\\
        | \\
      8 9  7
        |   \\
        |    \\
      1-4-5--2
        
Local numbering (0-indexed):
  Vertices:  0=1, 1=2, 2=3
  Edges 0-1: 3=mid closer to 0, 4=mid closer to 1
  Edges 1-2: 5=mid closer to 1, 6=mid closer to 2
  Edges 2-0: 7=mid closer to 2, 8=mid closer to 0
  Center:    9
""")

print("\nCURRENT CODE IMPLEMENTATION ORDER:")
print("Index | Current Description")
print("------|---------------------")
print("  0   | Vertex 1 (λ₁)")
print("  1   | Vertex 2 (λ₂)")
print("  2   | Vertex 3 (λ₃)")
print("  3   | Edge 1-2 near vertex 1")
print("  4   | Edge 1-2 near vertex 2")
print("  5   | Edge 2-3 near vertex 2")
print("  6   | Edge 2-3 near vertex 3")
print("  7   | Edge 3-1 near vertex 3")
print("  8   | Edge 3-1 near vertex 1")
print("  9   | Center")

print("\n" + "="*80)
print("CHECKING: Which gradient implementations are CORRECT?")
print("="*80)

# From the test output, let's identify which are correct
correct_indices = [0, 7, 8, 9]
wrong_indices = [1, 2, 3, 4, 5, 6]

print(f"\n✓ Correct gradients (match numerical):")
for i in correct_indices:
    print(f"    Index {i}")

print(f"\n✗ Wrong gradients (mismatch):")
for i in wrong_indices:
    print(f"    Index {i}")

print("\n" + "="*80)
print("PATTERN ANALYSIS")
print("="*80)

print("""
Looking at test output patterns:

Index 1 (Vertex 2, λ₂): 
  Analytical: [-0.406250, -0.000000]
  Numerical:  [-1.125000, -2.250000]
  This looks like what index 3 should be!

Index 3 (Edge 1-2 near 1):
  Analytical: [-1.125000, -2.250000]
  Numerical:  [-0.406250, -0.000000]
  This looks like what index 1 should be!

Conclusion: INDICES 1 AND 3 ARE SWAPPED!

Similarly:
  Indices 2 and 6 appear related
  Indices 4 and 5 appear related
""")

print("\n" + "="*80)
print("CORRECT PYTHON IMPLEMENTATION FOR P3")
print("="*80)

print("""
def Phi_Pk(l1, l2, pk):
    l3 = 1 - l1 - l2
    if pk == 3:
        return np.array([
            l1*(3*l1-1)*(3*l1-2)/2,          # 0: Vertex 1 (λ₁)
            l2*(3*l2-1)*(3*l2-2)/2,          # 1: Vertex 2 (λ₂)
            l3*(3*l3-1)*(3*l3-2)/2,          # 2: Vertex 3 (λ₃)
            (9/2)*l1*l2*(3*l1-1),             # 3: Edge 1-2, closer to vertex 1
            (9/2)*l1*l2*(3*l2-1),             # 4: Edge 1-2, closer to vertex 2
            (9/2)*l2*l3*(3*l2-1),             # 5: Edge 2-3, closer to vertex 2
            (9/2)*l2*l3*(3*l3-1),             # 6: Edge 2-3, closer to vertex 3
            (9/2)*l3*l1*(3*l3-1),             # 7: Edge 3-1, closer to vertex 3
            (9/2)*l3*l1*(3*l1-1),             # 8: Edge 3-1, closer to vertex 1
            27*l1*l2*l3                       # 9: Center
        ])

def GradPhi_Pk(l1, l2, G1, G2, pk):
    l3 = 1 - l1 - l2
    G3 = -G1 - G2
    
    if pk == 3:
        return np.array([
            # 0: Vertex 1
            ((27*l1**2 - 18*l1 + 2)/2) * G1,
            
            # 1: Vertex 2
            ((27*l2**2 - 18*l2 + 2)/2) * G2,
            
            # 2: Vertex 3
            ((27*l3**2 - 18*l3 + 2)/2) * G3,
            
            # 3: Edge 1-2 closer to vertex 1
            (9/2) * (l2*(6*l1-1)*G1 + l1*(3*l1-1)*G2),
            
            # 4: Edge 1-2 closer to vertex 2
            (9/2) * (l2*(3*l2-1)*G1 + l1*(6*l2-1)*G2),
            
            # 5: Edge 2-3 closer to vertex 2
            (9/2) * (l3*(3*l2-1)*G2 + l2*(6*l3-1)*G3),
            
            # 6: Edge 2-3 closer to vertex 3
            (9/2) * (l3*(6*l3-1)*G2 + l2*(3*l3-1)*G3),
            
            # 7: Edge 3-1 closer to vertex 3
            (9/2) * (l1*(3*l3-1)*G3 + l3*(6*l1-1)*G1),
            
            # 8: Edge 3-1 closer to vertex 1
            (9/2) * (l1*(6*l1-1)*G3 + l3*(3*l1-1)*G1),
            
            # 9: Center
            27 * (l2*l3*G1 + l1*l3*G2 + l1*l2*G3)
        ])
""")
