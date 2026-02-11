import numpy as np
from math import pi, sqrt
import sys
sys.path.insert(0, '/Users/charlesfabre/Desktop/Code/EDP_Elements_finis')
from Laplacian_Pk_2D_2025_V_AF import Phi_Pk, GradPhi_Pk, MGrad_Lambda_of_P1, Determinant

# =========================================================
# TEST P2 BASIS FUNCTIONS
# =========================================================

def test_P2_basis():
    """Test P2 basis functions by comparing analytical gradients with numerical gradients."""
    
    # Define a test triangle
    Xi = np.array([0.0, 0.0])
    Xj = np.array([1.0, 0.0])
    Xk = np.array([0.0, 1.0])
    
    # Compute gradient operators
    GradP1 = MGrad_Lambda_of_P1(Xi, Xj, Xk)
    G1 = GradP1[0, :]
    G2 = GradP1[1, :]
    
    print("Triangle vertices:")
    print(f"  Xi = {Xi}, Xj = {Xj}, Xk = {Xk}")
    print(f"\nGradient operators:")
    print(f"  G1 = {G1}")
    print(f"  G2 = {G2}")
    
    # Test several barycentric coordinates
    test_points = [
        (0.5, 0.25),
        (0.7, 0.2),
        (0.3, 0.3),
    ]
    
    epsilon = 1e-5  # For finite differences
    Ndl = 6  # Number of P2 basis functions
    
    errors = []
    
    for l1, l2 in test_points:
        l3 = 1 - l1 - l2
        print(f"\n{'='*70}")
        print(f"Testing at (λ1, λ2, λ3) = ({l1:.3f}, {l2:.3f}, {l3:.3f})")
        print(f"{'='*70}")
        
        # Analytical gradients
        GradPhi_analytical = GradPhi_Pk(l1, l2, G1, G2, pk=2)
        
        # Numerical gradients via finite differences
        GradPhi_numerical = np.zeros((Ndl, 2))
        
        for basis_idx in range(Ndl):
            # Phi at (l1 + eps, l2)
            Phi_plus_l1 = Phi_Pk(l1 + epsilon, l2, pk=2)[basis_idx]
            Phi_minus_l1 = Phi_Pk(l1 - epsilon, l2, pk=2)[basis_idx]
            
            # Phi at (l1, l2 + eps)
            Phi_plus_l2 = Phi_Pk(l1, l2 + epsilon, pk=2)[basis_idx]
            Phi_minus_l2 = Phi_Pk(l1, l2 - epsilon, pk=2)[basis_idx]
            
            # Numerical derivatives in barycentric space
            dPhi_dl1 = (Phi_plus_l1 - Phi_minus_l1) / (2 * epsilon)
            dPhi_dl2 = (Phi_plus_l2 - Phi_minus_l2) / (2 * epsilon)
            
            # Convert to physical gradient using G1 and G2
            GradPhi_numerical[basis_idx, :] = dPhi_dl1 * G1 + dPhi_dl2 * G2
        
        # Compare
        print(f"\n{'idx':<3} {'Gradient (analytical)':<35} {'Gradient (numerical)':<35} {'Error'}")
        print("-" * 110)
        
        max_error = 0
        for i in range(Ndl):
            analytical = GradPhi_analytical[i, :]
            numerical = GradPhi_numerical[i, :]
            error = np.linalg.norm(analytical - numerical)
            max_error = max(max_error, error)
            
            print(f"{i:<3} [{analytical[0]:+.6f}, {analytical[1]:+.6f}]  "
                  f"[{numerical[0]:+.6f}, {numerical[1]:+.6f}]  {error:.2e}")
        
        errors.append(max_error)
        print(f"\nMax error for this point: {max_error:.2e}")
        
        if max_error > 1e-4:
            print("⚠️  LARGE ERROR DETECTED!")
    
    print(f"\n{'='*70}")
    print(f"Summary: Max error across all test points: {max(errors):.2e}")
    print(f"{'='*70}")
    
    if max(errors) < 1e-4:
        print("✓ All tests passed! P2 basis function gradients are correct.")
    else:
        print("✗ Issues found in P2 basis function gradients!")


if __name__ == "__main__":
    test_P2_basis()
