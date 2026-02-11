import sympy as sp

# Define symbolic variables
l1, l2, l3 = sp.symbols('lambda_1 lambda_2 lambda_3', real=True)

print("="*80)
print("COMPUTING P2 LAGRANGE BASIS FUNCTION DERIVATIVES")
print("="*80)

# Define the 6 P2 basis functions in terms of l1, l2, l3
phi = {}

# Vertices
phi[0] = l1 * (2*l1 - 1)
phi[1] = l2 * (2*l2 - 1)
phi[2] = l3 * (2*l3 - 1)

# Edges (midpoints)
phi[3] = 4 * l1 * l2  # Edge 1-2
phi[4] = 4 * l2 * l3  # Edge 2-3
phi[5] = 4 * l3 * l1  # Edge 3-1

labels = [
    "Vertex 1 (λ₁)",
    "Vertex 2 (λ₂)",
    "Vertex 3 (λ₃)",
    "Edge 1-2 (midpoint)",
    "Edge 2-3 (midpoint)",
    "Edge 3-1 (midpoint)"
]

for i in range(6):
    print(f"\nBasis function φ₍{i}₎: {labels[i]}")
    print(f"Formula: φ = {phi[i]}")
    
    # Compute derivatives treating l1, l2, l3 as independent
    dphi_dl1 = sp.diff(phi[i], l1)
    dphi_dl2 = sp.diff(phi[i], l2)
    dphi_dl3 = sp.diff(phi[i], l3)
    
    # Simplify
    dphi_dl1 = sp.expand(dphi_dl1)
    dphi_dl2 = sp.expand(dphi_dl2)
    dphi_dl3 = sp.expand(dphi_dl3)
    
    print(f"\n∂φ/∂λ₁ = {dphi_dl1}")
    print(f"∂φ/∂λ₂ = {dphi_dl2}")
    print(f"∂φ/∂λ₃ = {dphi_dl3}")
    
    print("\n" + "-"*80)
    print(f"Gradient in barycentric space:")
    print(f"∇φ = (∂φ/∂λ₁)G₁ + (∂φ/∂λ₂)G₂ + (∂φ/∂λ₃)G₃")
    
    # Substitute G3 = -G1 - G2
    print(f"\nWith G₃ = -G₁ - G₂:")
    coeff_G1 = dphi_dl1 - dphi_dl3
    coeff_G2 = dphi_dl2 - dphi_dl3
    coeff_G1 = sp.expand(coeff_G1)
    coeff_G2 = sp.expand(coeff_G2)
    print(f"∇φ = ({coeff_G1})G₁ + ({coeff_G2})G₂")
    
    print("\n" + "="*80)
