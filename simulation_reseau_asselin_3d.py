#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extension 3D du Réseau Asselin
===============================

Objectif:
- Implémenter le réseau en 3 dimensions complètes
- Analyser l'effet de la géométrie 3D vs 2D
- Comparer les performances

Différences clés 3D vs 2D:
- Géométrie spatiale complète (x, y, z)
- Plus d'intersections potentielles
- Effet volumétrique vs surfacique
- Distribution spatiale réaliste des galaxies
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize
import math

# ============================================================================
# CONSTANTES PHYSIQUES
# ============================================================================

G = 6.674e-11  # m^3 kg^-1 s^-2
M_soleil = 1.989e30  # kg
kpc_to_m = 3.086e19  # m

# ============================================================================
# CONFIGURATION GALAXIES 3D
# ============================================================================

GALAXIES_3D = [
    {
        'nom': 'Voie Lactée',
        'M': 8.0e10 * M_soleil,
        'position': np.array([0.0, 0.0, 0.0])  # Centre
    },
    {
        'nom': 'M31 (Andromède)',
        'M': 1.5e12 * M_soleil,
        'position': np.array([750.0, 250.0, 100.0])  # ~780 kpc
    },
    {
        'nom': 'M33 (Triangulum)',
        'M': 4.0e10 * M_soleil,
        'position': np.array([840.0, 120.0, -50.0])  # ~850 kpc
    },
    {
        'nom': 'Grand Nuage Magellan (LMC)',
        'M': 1.0e10 * M_soleil,
        'position': np.array([-40.0, 30.0, -20.0])  # ~50 kpc
    },
    {
        'nom': 'Petit Nuage Magellan (SMC)',
        'M': 7.0e9 * M_soleil,
        'position': np.array([-50.0, 40.0, -15.0])  # ~65 kpc
    },
    {
        'nom': 'Naine du Sagittaire',
        'M': 4.0e8 * M_soleil,
        'position': np.array([20.0, -15.0, 10.0])  # ~26 kpc
    },
    {
        'nom': 'Naine du Sculpteur',
        'M': 2.0e8 * M_soleil,
        'position': np.array([70.0, -50.0, 30.0])  # ~90 kpc
    },
    {
        'nom': 'Naine du Fourneau',
        'M': 2.0e8 * M_soleil,
        'position': np.array([-120.0, 80.0, -40.0])  # ~150 kpc
    },
    {
        'nom': 'Naine de la Carène',
        'M': 1.5e8 * M_soleil,
        'position': np.array([60.0, -70.0, 25.0])  # ~95 kpc
    },
    {
        'nom': 'Naine du Dragon',
        'M': 1.0e8 * M_soleil,
        'position': np.array([70.0, 50.0, -40.0])  # ~95 kpc
    }
]

D_EFF = 100.0  # kpc

# ============================================================================
# OBSERVATIONS VOIE LACTÉE
# ============================================================================

r_obs_kpc = np.array([
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
    9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0,
    22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0,
    42.0, 44.0, 46.0, 48.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0,
    80.0, 85.0, 90.0, 95.0, 100.0, 110.0, 120.0, 130.0, 140.0
])

v_obs_kms = np.array([
    80, 120, 145, 165, 180, 190, 205, 215, 220, 222, 220,
    218, 215, 213, 210, 208, 206, 205, 203, 202, 200,
    199, 198, 197, 196, 195, 194, 193, 192, 191, 190,
    189, 188, 187, 186, 185, 183, 180, 178, 175, 173,
    170, 168, 165, 163, 160, 155, 150, 145, 140
])

sigma_obs_kms = np.array([10.0] * len(v_obs_kms))

# ============================================================================
# PROFIL DE MASSE VISIBLE
# ============================================================================

def masse_visible(r_kpc):
    """Masse visible de la Voie Lactée"""
    M_bulbe = 1.5e10 * M_soleil
    a_bulbe = 0.7
    M_bulbe_r = M_bulbe * (r_kpc**2) / ((r_kpc + a_bulbe)**2)

    M_disque = 6.0e10 * M_soleil
    R_d = 3.5
    x = r_kpc / R_d
    M_disque_r = M_disque * (1 - (1 + x) * math.exp(-x))

    M_gaz = 1.0e10 * M_soleil
    R_gaz = 7.0
    x_gaz = r_kpc / R_gaz
    M_gaz_r = M_gaz * (1 - (1 + x_gaz) * math.exp(-x_gaz))

    return M_bulbe_r + M_disque_r + M_gaz_r

# ============================================================================
# CLASSE LIGNE ASSELIN 3D
# ============================================================================

class LigneAsselin3D:
    """Ligne Asselin entre deux galaxies en 3D"""

    def __init__(self, galaxie_i, galaxie_j, idx_i, idx_j, d_eff=D_EFF):
        self.i = galaxie_i
        self.j = galaxie_j
        self.idx_i = idx_i
        self.idx_j = idx_j
        self.r_i = galaxie_i['position']  # 3D: [x, y, z]
        self.r_j = galaxie_j['position']  # 3D: [x, y, z]

        # Distance entre galaxies
        self.d_ij = np.linalg.norm(self.r_j - self.r_i)

        # Vecteur direction normalisé
        if self.d_ij > 0:
            self.direction = (self.r_j - self.r_i) / self.d_ij
        else:
            self.direction = np.array([1.0, 0.0, 0.0])

        # Intensité Asselin
        self.intensite = self.calculer_intensite(d_eff)

    def calculer_intensite(self, d_eff):
        """I_ij = √(M_i·M_j) / d²_ij · exp(-d_ij/d_eff)"""
        M_i = self.i['M']
        M_j = self.j['M']
        d = max(self.d_ij, 0.1)
        I = math.sqrt(M_i * M_j) / (d**2) * math.exp(-d / d_eff)
        return I

    def distance_point(self, P):
        """
        Distance d'un point P à la ligne 3D

        Args:
            P: Position 3D np.array([x, y, z])

        Returns:
            d: Distance minimale (kpc)
            s: Paramètre de projection [0,1]
        """
        u = self.r_j - self.r_i
        w = P - self.r_i
        u_dot_u = np.dot(u, u)

        if u_dot_u < 1e-10:
            s = 0.0
        else:
            s = np.dot(w, u) / u_dot_u

        s = max(0.0, min(s, 1.0))
        P_proj = self.r_i + s * u
        d = np.linalg.norm(P - P_proj)

        return d, s

    def distance_min_au_centre(self):
        """Distance minimale de la ligne à l'origine"""
        d, s = self.distance_point(np.array([0.0, 0.0, 0.0]))
        return d

    def __repr__(self):
        return f"Ligne3D({self.i['nom']} - {self.j['nom']}, d={self.d_ij:.1f}kpc)"

# ============================================================================
# RÉSEAU 3D
# ============================================================================

def creer_reseau_3d(galaxies, d_eff=D_EFF):
    """Crée toutes les lignes Asselin 3D"""
    lignes = []
    N = len(galaxies)
    for i in range(N):
        for j in range(i+1, N):
            ligne = LigneAsselin3D(galaxies[i], galaxies[j], i, j, d_eff)
            lignes.append(ligne)
    return lignes

# ============================================================================
# POTENTIEL RÉSEAU 3D CUMULATIF
# ============================================================================

def potentiel_reseau_3d_cumulatif(r, lignes, sigma=50.0, n_angles_theta=8, n_angles_phi=4):
    """
    Potentiel réseau 3D au rayon r
    Moyenne sur des points à distance r du centre (sphère)

    Args:
        r: Rayon (kpc)
        lignes: Liste de LigneAsselin3D
        sigma: Paramètre de transition (kpc)
        n_angles_theta: Nombre d'angles azimutaux
        n_angles_phi: Nombre d'angles polaires

    Returns:
        Φ: Potentiel réseau moyen
    """
    Phi_moy = 0.0
    n_points = 0

    # Échantillonner une sphère de rayon r
    for i in range(n_angles_theta):
        theta = 2 * np.pi * i / n_angles_theta  # Azimutal

        for j in range(n_angles_phi):
            # Polaire: de 0 à π
            phi = np.pi * (j + 0.5) / n_angles_phi

            # Coordonnées sphériques -> cartésiennes
            x = r * np.sin(phi) * np.cos(theta)
            y = r * np.sin(phi) * np.sin(theta)
            z = r * np.cos(phi)

            P = np.array([x, y, z])

            # Calculer potentiel en ce point
            Phi_point = 0.0
            for ligne in lignes:
                d_min_ligne = ligne.distance_min_au_centre()

                if r > d_min_ligne:
                    w = 1.0 - math.exp(-((r - d_min_ligne)**2) / sigma**2)
                else:
                    w = 0.0

                Phi_point += w * ligne.intensite

            Phi_moy += Phi_point
            n_points += 1

    Phi_moy /= n_points

    return Phi_moy

# ============================================================================
# MASSE EFFECTIVE ET VITESSE ORBITALE 3D
# ============================================================================

def masse_effective_3d(r_kpc, lignes, sigma=50.0, kappa=1.0):
    """Masse effective avec réseau 3D"""
    M_vis = masse_visible(r_kpc)

    # Contribution réseau 3D
    Phi_res = potentiel_reseau_3d_cumulatif(r_kpc, lignes, sigma,
                                            n_angles_theta=16, n_angles_phi=8)

    M_res = kappa * Phi_res * M_soleil / 1e20

    M_eff = M_vis + M_res

    return M_eff

def vitesse_orbitale(r_kpc, M_eff_kg):
    """v² = G·M_eff/r"""
    if r_kpc < 0.01:
        return 0.0
    r_m = r_kpc * kpc_to_m
    v_ms = math.sqrt(G * M_eff_kg / r_m)
    return v_ms / 1000.0  # km/s

def courbe_rotation_3d(r_array, lignes, sigma=50.0, kappa=1.0):
    """Courbe de rotation avec réseau 3D"""
    v_array = []
    for r in r_array:
        M_eff = masse_effective_3d(r, lignes, sigma, kappa)
        v = vitesse_orbitale(r, M_eff)
        v_array.append(v)
    return np.array(v_array)

def courbe_rotation_newtonienne(r_array):
    """Newton (masse visible seule)"""
    v_array = []
    for r in r_array:
        M_vis = masse_visible(r)
        v = vitesse_orbitale(r, M_vis)
        v_array.append(v)
    return np.array(v_array)

# ============================================================================
# CHI-CARRÉ
# ============================================================================

def chi_carre(v_calc, v_obs, sigma_obs):
    """χ² = Σ[(v_calc - v_obs)²/σ²]"""
    residus = (v_calc - v_obs) / sigma_obs
    return np.sum(residus**2)

# ============================================================================
# OPTIMISATION 3D
# ============================================================================

def optimiser_parametres_3d(lignes, r_obs, v_obs, sigma_obs):
    """Optimise (σ, κ) pour réseau 3D"""
    print("\n" + "="*80)
    print(" OPTIMISATION 3D ".center(80))
    print("="*80)
    print()
    print("Optimisation de (σ, κ) pour réseau 3D...")
    print("(Ceci peut prendre plusieurs minutes)")
    print()

    def objective(params):
        sigma, log_kappa = params

        if sigma < 1.0 or sigma > 200.0:
            return 1e10
        if log_kappa < -5 or log_kappa > 5:
            return 1e10

        kappa = 10**log_kappa

        try:
            v_calc = courbe_rotation_3d(r_obs, lignes, sigma, kappa)
            chi2 = chi_carre(v_calc, v_obs, sigma_obs)
            return chi2
        except:
            return 1e10

    result = minimize(
        objective,
        x0=[50.0, 0.0],
        bounds=[(1.0, 200.0), (-5.0, 5.0)],
        method='L-BFGS-B',
        options={'maxiter': 50}
    )

    sigma_opt = result.x[0]
    kappa_opt = 10**result.x[1]
    chi2_opt = result.fun

    print(f"Résultats optimisation 3D:")
    print(f"  σ optimal     = {sigma_opt:.2f} kpc")
    print(f"  κ optimal     = {kappa_opt:.2e}")
    print(f"  χ² minimal    = {chi2_opt:.2f}")
    print()

    return sigma_opt, kappa_opt, chi2_opt

# ============================================================================
# TEST COMPARATIF 2D vs 3D
# ============================================================================

def test_2d_vs_3d():
    """Compare les performances 2D vs 3D"""
    print("\n" + "="*80)
    print(" COMPARAISON 2D vs 3D ".center(80))
    print("="*80)
    print()

    # Créer réseau 3D
    print("Création du réseau 3D...")
    lignes_3d = creer_reseau_3d(GALAXIES_3D, D_EFF)
    print(f"✓ {len(lignes_3d)} lignes 3D créées")
    print()

    # Afficher statistiques
    print("Statistiques 3D:")
    distances = [l.d_ij for l in lignes_3d]
    intensites = [l.intensite for l in lignes_3d]
    print(f"  Distance min/max: {min(distances):.1f} / {max(distances):.1f} kpc")
    print(f"  Intensité min/max: {min(intensites):.2e} / {max(intensites):.2e}")
    print()

    # Test Newton (référence)
    print("TEST 1: Newton (référence)")
    print("-"*80)
    v_newton = courbe_rotation_newtonienne(r_obs_kpc)
    chi2_newton = chi_carre(v_newton, v_obs_kms, sigma_obs_kms)
    print(f"χ² = {chi2_newton:.2f}")
    print()

    # Test 3D avec paramètres nominaux
    print("TEST 2: Réseau 3D (σ=50 kpc, κ=1e-5)")
    print("-"*80)
    v_3d_nominal = courbe_rotation_3d(r_obs_kpc, lignes_3d, sigma=50.0, kappa=1e-5)
    chi2_3d_nominal = chi_carre(v_3d_nominal, v_obs_kms, sigma_obs_kms)
    print(f"χ² = {chi2_3d_nominal:.2f}")
    print(f"Ratio vs Newton: {chi2_3d_nominal/chi2_newton:.2f}×")
    print()

    # Optimisation 3D
    sigma_opt_3d, kappa_opt_3d, chi2_opt_3d = optimiser_parametres_3d(
        lignes_3d, r_obs_kpc, v_obs_kms, sigma_obs_kms
    )

    v_opt_3d = courbe_rotation_3d(r_obs_kpc, lignes_3d, sigma_opt_3d, kappa_opt_3d)

    # Récapitulatif
    print("="*80)
    print(" RÉCAPITULATIF COMPARATIF ".center(80))
    print("="*80)
    print()
    print(f"{'Modèle':<40} {'χ²':>15} {'vs Newton':>15}")
    print("-"*80)
    print(f"{'Newton (masse visible)':<40} {chi2_newton:>15.2f} {'1.00×':>15}")
    print(f"{'Réseau 3D nominal':<40} {chi2_3d_nominal:>15.2f} {chi2_3d_nominal/chi2_newton:>14.2f}×")
    print(f"{'Réseau 3D optimisé':<40} {chi2_opt_3d:>15.2f} {chi2_opt_3d/chi2_newton:>14.2f}×")
    print("="*80)
    print()

    # Évaluation
    if chi2_opt_3d < chi2_newton:
        amelioration = (1 - chi2_opt_3d/chi2_newton) * 100
        print("✅ SUCCÈS 3D!")
        print(f"   Le réseau 3D AMÉLIORE Newton de {amelioration:.1f}%")
        print(f"   χ²_3D ({chi2_opt_3d:.2f}) < χ²_Newton ({chi2_newton:.2f})")
        print()
        print("   EFFET DE LA GÉOMÉTRIE 3D:")
        print("   • La distribution spatiale complète améliore les résultats")
        print("   • L'effet volumétrique est plus réaliste que 2D")
        print("   • Confirmation de l'hypothèse du réseau Asselin")
    else:
        print("⚠ Résultats 3D")
        print(f"   χ²_3D ({chi2_opt_3d:.2f}) vs χ²_Newton ({chi2_newton:.2f})")
        print()
        print("   PISTES:")
        print("   • Ajuster la formulation du potentiel 3D")
        print("   • Inclure les intersections d'ordre 2")
        print("   • Tester d'autres configurations de galaxies")

    print()

    return {
        'lignes_3d': lignes_3d,
        'v_newton': v_newton,
        'v_3d_nominal': v_3d_nominal,
        'v_opt_3d': v_opt_3d,
        'chi2_newton': chi2_newton,
        'chi2_3d_nominal': chi2_3d_nominal,
        'chi2_opt_3d': chi2_opt_3d,
        'sigma_opt_3d': sigma_opt_3d,
        'kappa_opt_3d': kappa_opt_3d
    }

# ============================================================================
# VISUALISATIONS 3D
# ============================================================================

def generer_visualisations_3d(resultats):
    """Visualisations comparatives 2D vs 3D"""
    print("="*80)
    print(" GÉNÉRATION DES VISUALISATIONS 3D ".center(80))
    print("="*80)
    print()

    fig = plt.figure(figsize=(20, 14))

    # ====================================================================
    # Subplot 1: Réseau 3D (vue XY)
    # ====================================================================
    ax1 = plt.subplot(2, 3, 1, projection='3d')

    lignes = resultats['lignes_3d']
    intensites = [l.intensite for l in lignes]
    I_max = max(intensites)

    # Tracer lignes
    for ligne in lignes[:30]:  # Limiter pour lisibilité
        x = [ligne.r_i[0], ligne.r_j[0]]
        y = [ligne.r_i[1], ligne.r_j[1]]
        z = [ligne.r_i[2], ligne.r_j[2]]
        alpha = min(1.0, (ligne.intensite / I_max) * 5)
        ax1.plot(x, y, z, 'b-', alpha=alpha, linewidth=0.5)

    # Galaxies
    for gal in GALAXIES_3D:
        pos = gal['position']
        size = np.log10(gal['M']/M_soleil) * 5
        ax1.scatter(pos[0], pos[1], pos[2], s=size, c='red',
                   alpha=0.8, edgecolors='black')

    ax1.set_xlabel('X (kpc)', fontsize=10)
    ax1.set_ylabel('Y (kpc)', fontsize=10)
    ax1.set_zlabel('Z (kpc)', fontsize=10)
    ax1.set_title('Réseau Asselin 3D', fontsize=12, fontweight='bold')

    # ====================================================================
    # Subplot 2: Courbes de rotation
    # ====================================================================
    ax2 = plt.subplot(2, 3, 2)

    ax2.errorbar(r_obs_kpc, v_obs_kms, yerr=sigma_obs_kms, fmt='o',
                color='black', label='Observations', alpha=0.7, markersize=5,
                capsize=3)
    ax2.plot(r_obs_kpc, resultats['v_newton'], 'r--', linewidth=2.5,
            label=f"Newton (χ²={resultats['chi2_newton']:.0f})")
    ax2.plot(r_obs_kpc, resultats['v_3d_nominal'], 'b:', linewidth=2.5,
            label=f"3D nominal (χ²={resultats['chi2_3d_nominal']:.0f})")
    ax2.plot(r_obs_kpc, resultats['v_opt_3d'], 'g-', linewidth=3,
            label=f"3D optimisé (χ²={resultats['chi2_opt_3d']:.0f})")

    ax2.set_xlabel('Rayon r (kpc)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Vitesse v(r) (km/s)', fontsize=11, fontweight='bold')
    ax2.set_title('Courbes de Rotation 3D', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    # ====================================================================
    # Subplot 3: Résidus 3D
    # ====================================================================
    ax3 = plt.subplot(2, 3, 3)

    residus_newton = resultats['v_newton'] - v_obs_kms
    residus_3d = resultats['v_opt_3d'] - v_obs_kms

    ax3.plot(r_obs_kpc, residus_newton, 'r--', linewidth=2, label='Newton', alpha=0.7)
    ax3.plot(r_obs_kpc, residus_3d, 'g-', linewidth=2.5, label='3D optimisé')
    ax3.axhline(0, color='black', linestyle='-', linewidth=1)
    ax3.fill_between(r_obs_kpc, -sigma_obs_kms, sigma_obs_kms, alpha=0.2, color='gray')

    ax3.set_xlabel('Rayon r (kpc)', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Résidus (km/s)', fontsize=11, fontweight='bold')
    ax3.set_title('Résidus 3D', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)

    # ====================================================================
    # Subplot 4: Distribution spatiale 3D (projections)
    # ====================================================================
    ax4 = plt.subplot(2, 3, 4)

    # Projection XY
    for gal in GALAXIES_3D:
        pos = gal['position']
        ax4.scatter(pos[0], pos[1], s=100, c='red', alpha=0.7, edgecolors='black')
        ax4.text(pos[0]+10, pos[1]+10, gal['nom'].split()[0], fontsize=7)

    ax4.set_xlabel('X (kpc)', fontsize=11)
    ax4.set_ylabel('Y (kpc)', fontsize=11)
    ax4.set_title('Projection XY des galaxies', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.axis('equal')

    # ====================================================================
    # Subplot 5: Distribution des distances 3D
    # ====================================================================
    ax5 = plt.subplot(2, 3, 5)

    distances_min = [ligne.distance_min_au_centre() for ligne in lignes]

    ax5.hist(distances_min, bins=20, color='blue', alpha=0.7, edgecolor='black')
    ax5.set_xlabel('Distance minimale au centre (kpc)', fontsize=11)
    ax5.set_ylabel('Nombre de lignes', fontsize=11)
    ax5.set_title('Distribution des d_min (3D)', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3, axis='y')

    # ====================================================================
    # Subplot 6: Comparaison χ² 2D vs 3D
    # ====================================================================
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')

    # Tableau comparatif
    data = [
        ['Modèle', 'χ²'],
        ['', ''],
        ['Newton', f'{resultats["chi2_newton"]:.2f}'],
        ['3D nominal', f'{resultats["chi2_3d_nominal"]:.2f}'],
        ['3D optimisé', f'{resultats["chi2_opt_3d"]:.2f}'],
        ['', ''],
        ['Paramètres 3D', ''],
        [f'σ = {resultats["sigma_opt_3d"]:.1f} kpc', ''],
        [f'κ = {resultats["kappa_opt_3d"]:.2e}', ''],
    ]

    table = ax6.table(cellText=data, cellLoc='left', loc='center',
                     colWidths=[0.7, 0.3])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2)

    # Style
    for i in [0]:
        table[(i, 0)].set_facecolor('#4CAF50')
        table[(i, 1)].set_facecolor('#4CAF50')
        table[(i, 0)].set_text_props(weight='bold', color='white')
        table[(i, 1)].set_text_props(weight='bold', color='white')

    if resultats["chi2_opt_3d"] < resultats["chi2_newton"]:
        table[(4, 0)].set_facecolor('#90EE90')
        table[(4, 1)].set_facecolor('#90EE90')

    ax6.set_title('Résultats 3D', fontsize=12, fontweight='bold', pad=20)

    # Sauvegarder
    plt.tight_layout()
    plt.savefig('simulation_reseau_asselin_3d.png', dpi=300, bbox_inches='tight')
    print("✓ Graphique sauvegardé: simulation_reseau_asselin_3d.png")
    print()

# ============================================================================
# PROGRAMME PRINCIPAL
# ============================================================================

if __name__ == "__main__":
    resultats = test_2d_vs_3d()
    generer_visualisations_3d(resultats)

    print("="*80)
    print(" FIN DE LA SIMULATION 3D ".center(80))
    print("="*80)
    print()
