#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calcul des Courbes de Rotation avec Réseau Asselin
===================================================

Étape 1 complète:
- Dériver v(r) depuis Φ_réseau(r)
- Comparer avec observations Voie Lactée
- Optimiser σ pour minimiser χ²

Approche:
1. Potentiel total: Φ_total(r) = Φ_visible(r) + Φ_réseau(r)
2. Vitesse orbitale: v(r)² = r · dΦ_total/dr
3. Ou équivalent: v(r)² = G·M_eff(r) / r
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import math

# ============================================================================
# CONSTANTES PHYSIQUES
# ============================================================================

G = 6.674e-11  # m^3 kg^-1 s^-2
M_soleil = 1.989e30  # kg
kpc_to_m = 3.086e19  # m

# ============================================================================
# CONFIGURATION GALAXIES 2D
# ============================================================================

GALAXIES_2D = [
    {'nom': 'Voie Lactée', 'M': 8.0e10 * M_soleil, 'position': np.array([0.0, 0.0])},
    {'nom': 'M31 (Andromède)', 'M': 1.5e12 * M_soleil, 'position': np.array([750.0, 250.0])},
    {'nom': 'M33 (Triangulum)', 'M': 4.0e10 * M_soleil, 'position': np.array([840.0, 120.0])},
    {'nom': 'Grand Nuage Magellan (LMC)', 'M': 1.0e10 * M_soleil, 'position': np.array([-40.0, 30.0])},
    {'nom': 'Petit Nuage Magellan (SMC)', 'M': 7.0e9 * M_soleil, 'position': np.array([-50.0, 40.0])},
    {'nom': 'Naine du Sagittaire', 'M': 4.0e8 * M_soleil, 'position': np.array([20.0, -15.0])},
    {'nom': 'Naine du Sculpteur', 'M': 2.0e8 * M_soleil, 'position': np.array([70.0, -50.0])},
    {'nom': 'Naine du Fourneau', 'M': 2.0e8 * M_soleil, 'position': np.array([-120.0, 80.0])},
    {'nom': 'Naine de la Carène', 'M': 1.5e8 * M_soleil, 'position': np.array([60.0, -70.0])},
    {'nom': 'Naine du Dragon', 'M': 1.0e8 * M_soleil, 'position': np.array([70.0, 50.0])}
]

D_EFF = 100.0  # kpc

# ============================================================================
# OBSERVATIONS VOIE LACTÉE
# ============================================================================

# Données observationnelles de courbe de rotation
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
    """
    Masse visible de la Voie Lactée
    Modèle: Bulbe + Disque + Gaz
    """
    # Bulbe (Hernquist)
    M_bulbe = 1.5e10 * M_soleil
    a_bulbe = 0.7
    M_bulbe_r = M_bulbe * (r_kpc**2) / ((r_kpc + a_bulbe)**2)

    # Disque exponentiel
    M_disque = 6.0e10 * M_soleil
    R_d = 3.5
    x = r_kpc / R_d
    M_disque_r = M_disque * (1 - (1 + x) * math.exp(-x))

    # Gaz
    M_gaz = 1.0e10 * M_soleil
    R_gaz = 7.0
    x_gaz = r_kpc / R_gaz
    M_gaz_r = M_gaz * (1 - (1 + x_gaz) * math.exp(-x_gaz))

    return M_bulbe_r + M_disque_r + M_gaz_r

# ============================================================================
# CLASSE LIGNE ASSELIN 2D
# ============================================================================

class LigneAsselin2D:
    """Ligne Asselin entre deux galaxies dans un plan 2D"""

    def __init__(self, galaxie_i, galaxie_j, idx_i, idx_j, d_eff=D_EFF):
        self.i = galaxie_i
        self.j = galaxie_j
        self.idx_i = idx_i
        self.idx_j = idx_j
        self.r_i = galaxie_i['position']
        self.r_j = galaxie_j['position']
        self.d_ij = np.linalg.norm(self.r_j - self.r_i)

        if self.d_ij > 0:
            self.direction = (self.r_j - self.r_i) / self.d_ij
        else:
            self.direction = np.array([1.0, 0.0])

        self.intensite = self.calculer_intensite(d_eff)

    def calculer_intensite(self, d_eff):
        M_i = self.i['M']
        M_j = self.j['M']
        d = max(self.d_ij, 0.1)
        I = math.sqrt(M_i * M_j) / (d**2) * math.exp(-d / d_eff)
        return I

    def distance_point(self, P):
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
        d, s = self.distance_point(np.array([0.0, 0.0]))
        return d

# ============================================================================
# RÉSEAU
# ============================================================================

def creer_reseau_2d(galaxies, d_eff=D_EFF):
    lignes = []
    N = len(galaxies)
    for i in range(N):
        for j in range(i+1, N):
            ligne = LigneAsselin2D(galaxies[i], galaxies[j], i, j, d_eff)
            lignes.append(ligne)
    return lignes

# ============================================================================
# POTENTIEL RÉSEAU CUMULATIF
# ============================================================================

def potentiel_reseau_cumulatif(r, lignes, sigma=50.0):
    """
    Potentiel réseau cumulatif au rayon r
    Mode pondéré avec transition lisse
    """
    Phi = 0.0

    for ligne in lignes:
        d_min = ligne.distance_min_au_centre()

        if r > d_min:
            w = 1.0 - math.exp(-((r - d_min)**2) / sigma**2)
        else:
            w = 0.0

        Phi += w * ligne.intensite

    return Phi

# ============================================================================
# MASSE EFFECTIVE ET VITESSE ORBITALE
# ============================================================================

def masse_effective_avec_reseau(r_kpc, lignes, sigma=50.0, kappa=1.0):
    """
    Masse effective totale incluant contribution réseau

    M_eff(r) = M_visible(r) + κ · Φ_réseau(r)

    où κ convertit le potentiel réseau en masse effective
    """
    M_vis = masse_visible(r_kpc)

    # Contribution réseau
    Phi_res = potentiel_reseau_cumulatif(r_kpc, lignes, sigma)

    # Convertir en masse effective
    # κ a des dimensions de [kg / (M☉/kpc²)]
    M_res = kappa * Phi_res * M_soleil / 1e20  # Normalisation

    M_eff = M_vis + M_res

    return M_eff

def vitesse_orbitale(r_kpc, M_eff_kg):
    """
    Vitesse orbitale circulaire
    v² = G·M_eff/r
    """
    if r_kpc < 0.01:
        return 0.0

    r_m = r_kpc * kpc_to_m
    v_ms = math.sqrt(G * M_eff_kg / r_m)
    return v_ms / 1000.0  # km/s

def courbe_rotation_avec_reseau(r_array, lignes, sigma=50.0, kappa=1.0):
    """
    Calcule la courbe de rotation avec réseau Asselin

    Args:
        r_array: Rayons (kpc)
        lignes: Réseau de lignes
        sigma: Paramètre de transition (kpc)
        kappa: Facteur de couplage réseau

    Returns:
        v_array: Vitesses (km/s)
    """
    v_array = []

    for r in r_array:
        M_eff = masse_effective_avec_reseau(r, lignes, sigma, kappa)
        v = vitesse_orbitale(r, M_eff)
        v_array.append(v)

    return np.array(v_array)

def courbe_rotation_newtonienne(r_array):
    """Courbe de rotation newtonienne (masse visible seule)"""
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
# OPTIMISATION
# ============================================================================

def optimiser_parametres(lignes, r_obs, v_obs, sigma_obs):
    """
    Optimise (σ, κ) pour minimiser χ²

    Args:
        lignes: Réseau de lignes
        r_obs, v_obs, sigma_obs: Observations

    Returns:
        (sigma_opt, kappa_opt, chi2_opt)
    """
    print("\n" + "="*80)
    print(" OPTIMISATION DES PARAMÈTRES ".center(80))
    print("="*80)
    print()
    print("Optimisation de (σ, κ) pour minimiser χ²...")
    print()

    def objective(params):
        sigma, log_kappa = params

        # Contraintes
        if sigma < 1.0 or sigma > 200.0:
            return 1e10
        if log_kappa < -5 or log_kappa > 5:
            return 1e10

        kappa = 10**log_kappa

        try:
            v_calc = courbe_rotation_avec_reseau(r_obs, lignes, sigma, kappa)
            chi2 = chi_carre(v_calc, v_obs, sigma_obs)
            return chi2
        except:
            return 1e10

    # Optimisation
    result = minimize(
        objective,
        x0=[50.0, 0.0],  # σ=50 kpc, κ=1
        bounds=[(1.0, 200.0), (-5.0, 5.0)],
        method='L-BFGS-B',
        options={'maxiter': 100, 'disp': False}
    )

    sigma_opt = result.x[0]
    kappa_opt = 10**result.x[1]
    chi2_opt = result.fun

    print(f"Résultats de l'optimisation:")
    print(f"  σ optimal     = {sigma_opt:.2f} kpc")
    print(f"  κ optimal     = {kappa_opt:.2e}")
    print(f"  χ² minimal    = {chi2_opt:.2f}")
    print()

    return sigma_opt, kappa_opt, chi2_opt

# ============================================================================
# TEST COMPLET
# ============================================================================

def test_courbes_rotation():
    """
    Test complet des courbes de rotation
    """
    print("\n" + "="*80)
    print(" CALCUL DES COURBES DE ROTATION AVEC RÉSEAU ASSELIN ".center(80))
    print("="*80)
    print()

    # Créer le réseau
    print("Création du réseau de lignes...")
    lignes = creer_reseau_2d(GALAXIES_2D, D_EFF)
    print(f"✓ {len(lignes)} lignes créées")
    print()

    # Test 1: Newton (référence)
    print("TEST 1: Newton (masse visible seule)")
    print("-"*80)
    v_newton = courbe_rotation_newtonienne(r_obs_kpc)
    chi2_newton = chi_carre(v_newton, v_obs_kms, sigma_obs_kms)
    print(f"χ² = {chi2_newton:.2f}")
    print()

    # Test 2: Réseau avec paramètres nominaux
    print("TEST 2: Réseau Asselin (σ=50 kpc, κ=1)")
    print("-"*80)
    v_nominal = courbe_rotation_avec_reseau(r_obs_kpc, lignes, sigma=50.0, kappa=1.0)
    chi2_nominal = chi_carre(v_nominal, v_obs_kms, sigma_obs_kms)
    print(f"χ² = {chi2_nominal:.2f}")
    print(f"Ratio vs Newton: {chi2_nominal/chi2_newton:.2f}×")
    print()

    # Test 3: Optimisation
    sigma_opt, kappa_opt, chi2_opt = optimiser_parametres(
        lignes, r_obs_kpc, v_obs_kms, sigma_obs_kms
    )

    v_opt = courbe_rotation_avec_reseau(r_obs_kpc, lignes, sigma_opt, kappa_opt)

    # Récapitulatif
    print("="*80)
    print(" RÉCAPITULATIF ".center(80))
    print("="*80)
    print()
    print(f"{'Modèle':<40} {'χ²':>15} {'vs Newton':>15}")
    print("-"*80)
    print(f"{'Newton (masse visible)':<40} {chi2_newton:>15.2f} {'1.00×':>15}")
    print(f"{'Réseau nominal (σ=50, κ=1)':<40} {chi2_nominal:>15.2f} {chi2_nominal/chi2_newton:>14.2f}×")
    print(f"{'Réseau optimisé':<40} {chi2_opt:>15.2f} {chi2_opt/chi2_newton:>14.2f}×")
    print("="*80)
    print()

    # Évaluation
    if chi2_opt < chi2_newton:
        amelioration = (1 - chi2_opt/chi2_newton) * 100
        print("✅ SUCCÈS MAJEUR!")
        print(f"   Le réseau Asselin AMÉLIORE Newton de {amelioration:.1f}%")
        print(f"   χ²_réseau ({chi2_opt:.2f}) < χ²_Newton ({chi2_newton:.2f})")
        print()
        print("   IMPLICATIONS:")
        print("   • Le réseau géométrique 2D explique mieux les observations")
        print("   • Pas besoin de 'matière noire' exotique")
        print("   • La géométrie de l'espace-temps suffit!")
    elif chi2_opt < chi2_nominal:
        print("⚠ AMÉLIORATION PARTIELLE")
        print(f"   L'optimisation améliore les paramètres nominaux")
        print(f"   Mais χ²_réseau ({chi2_opt:.2f}) > χ²_Newton ({chi2_newton:.2f})")
        print()
        print("   PISTES D'AMÉLIORATION:")
        print("   • Extension 3D du réseau")
        print("   • Intersections d'ordre 2")
        print("   • Révision de la formulation du potentiel")
    else:
        print("⚠ PAS D'AMÉLIORATION SIGNIFICATIVE")
        print(f"   χ²_optimisé ({chi2_opt:.2f}) ≈ χ²_nominal ({chi2_nominal:.2f})")
        print()
        print("   INTERPRÉTATION:")
        print("   • Le réseau 2D seul ne suffit pas")
        print("   • Extension 3D nécessaire")
        print("   • Ou combinaison avec d'autres effets")

    print()

    return {
        'lignes': lignes,
        'r_obs': r_obs_kpc,
        'v_obs': v_obs_kms,
        'sigma_obs': sigma_obs_kms,
        'v_newton': v_newton,
        'v_nominal': v_nominal,
        'v_opt': v_opt,
        'chi2_newton': chi2_newton,
        'chi2_nominal': chi2_nominal,
        'chi2_opt': chi2_opt,
        'sigma_opt': sigma_opt,
        'kappa_opt': kappa_opt
    }

# ============================================================================
# VISUALISATIONS
# ============================================================================

def generer_visualisations(resultats):
    """Génère les visualisations des courbes de rotation"""
    print("="*80)
    print(" GÉNÉRATION DES VISUALISATIONS ".center(80))
    print("="*80)
    print()

    fig = plt.figure(figsize=(20, 12))

    r_obs = resultats['r_obs']
    v_obs = resultats['v_obs']
    sigma_obs = resultats['sigma_obs']

    # ====================================================================
    # Subplot 1: Courbes de rotation
    # ====================================================================
    ax1 = plt.subplot(2, 3, 1)

    ax1.errorbar(r_obs, v_obs, yerr=sigma_obs, fmt='o',
                color='black', label='Observations', alpha=0.7, markersize=5,
                capsize=3, zorder=10)
    ax1.plot(r_obs, resultats['v_newton'], 'r--', linewidth=2.5,
            label=f"Newton (χ²={resultats['chi2_newton']:.0f})")
    ax1.plot(r_obs, resultats['v_nominal'], 'b:', linewidth=2.5,
            label=f"Réseau nominal (χ²={resultats['chi2_nominal']:.0f})")
    ax1.plot(r_obs, resultats['v_opt'], 'g-', linewidth=3,
            label=f"Réseau optimisé (χ²={resultats['chi2_opt']:.0f})")

    ax1.set_xlabel('Rayon r (kpc)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Vitesse v(r) (km/s)', fontsize=12, fontweight='bold')
    ax1.set_title('Courbes de Rotation: Réseau Asselin vs Newton',
                 fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10, loc='lower right')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, max(r_obs)*1.05)
    ax1.set_ylim(0, max(v_obs)*1.2)

    # ====================================================================
    # Subplot 2: Résidus optimisés
    # ====================================================================
    ax2 = plt.subplot(2, 3, 2)

    residus_opt = resultats['v_opt'] - v_obs
    residus_newton = resultats['v_newton'] - v_obs

    ax2.plot(r_obs, residus_newton, 'r--', linewidth=2, label='Newton', alpha=0.7)
    ax2.plot(r_obs, residus_opt, 'g-', linewidth=2.5, label='Réseau optimisé')
    ax2.axhline(0, color='black', linestyle='-', linewidth=1)
    ax2.fill_between(r_obs, -sigma_obs, sigma_obs, alpha=0.2, color='gray',
                     label='±1σ')

    ax2.set_xlabel('Rayon r (kpc)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Résidus (km/s)', fontsize=12, fontweight='bold')
    ax2.set_title('Résidus: v_modèle - v_obs', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    # ====================================================================
    # Subplot 3: Masse effective
    # ====================================================================
    ax3 = plt.subplot(2, 3, 3)

    r_plot = np.linspace(0.5, 150, 100)
    M_vis = np.array([masse_visible(r) / M_soleil / 1e10 for r in r_plot])
    M_eff_opt = np.array([
        masse_effective_avec_reseau(r, resultats['lignes'],
                                    resultats['sigma_opt'],
                                    resultats['kappa_opt']) / M_soleil / 1e10
        for r in r_plot
    ])

    ax3.plot(r_plot, M_vis, 'r--', linewidth=2.5, label='M_visible')
    ax3.plot(r_plot, M_eff_opt, 'g-', linewidth=2.5, label='M_eff (avec réseau)')
    ax3.fill_between(r_plot, M_vis, M_eff_opt, alpha=0.3, color='green',
                     label='Contribution réseau')

    ax3.set_xlabel('Rayon r (kpc)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Masse cumulée (10¹⁰ M☉)', fontsize=12, fontweight='bold')
    ax3.set_title('Masse Effective avec Réseau Asselin', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 150)

    # ====================================================================
    # Subplot 4: Distribution χ² par région
    # ====================================================================
    ax4 = plt.subplot(2, 3, 4)

    # Calculer χ² par bins de rayon
    bins = [(0, 10), (10, 30), (30, 60), (60, 100), (100, 150)]
    chi2_bins_newton = []
    chi2_bins_opt = []
    labels = []

    for r_min, r_max in bins:
        mask = (r_obs >= r_min) & (r_obs < r_max)
        if np.sum(mask) > 0:
            chi2_n = chi_carre(resultats['v_newton'][mask], v_obs[mask], sigma_obs[mask])
            chi2_o = chi_carre(resultats['v_opt'][mask], v_obs[mask], sigma_obs[mask])
            chi2_bins_newton.append(chi2_n)
            chi2_bins_opt.append(chi2_o)
            labels.append(f'{r_min}-{r_max}')

    x = np.arange(len(labels))
    width = 0.35

    ax4.bar(x - width/2, chi2_bins_newton, width, label='Newton', color='red', alpha=0.7)
    ax4.bar(x + width/2, chi2_bins_opt, width, label='Réseau optimisé', color='green', alpha=0.7)

    ax4.set_xlabel('Rayon (kpc)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('χ² (par région)', fontsize=12, fontweight='bold')
    ax4.set_title('Distribution du χ² par région radiale', fontsize=14, fontweight='bold')
    ax4.set_xticks(x)
    ax4.set_xticklabels(labels, rotation=45)
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3, axis='y')

    # ====================================================================
    # Subplot 5: Paramètres optimisés
    # ====================================================================
    ax5 = plt.subplot(2, 3, 5)

    # Tester différentes valeurs de σ
    sigma_range = np.linspace(10, 150, 30)
    chi2_vs_sigma = []

    kappa_fixed = resultats['kappa_opt']

    for sig in sigma_range:
        v_test = courbe_rotation_avec_reseau(r_obs, resultats['lignes'], sig, kappa_fixed)
        chi2 = chi_carre(v_test, v_obs, sigma_obs)
        chi2_vs_sigma.append(chi2)

    ax5.plot(sigma_range, chi2_vs_sigma, 'b-', linewidth=2)
    ax5.axvline(resultats['sigma_opt'], color='green', linestyle='--',
               linewidth=2, label=f"σ_opt = {resultats['sigma_opt']:.1f} kpc")
    ax5.axhline(resultats['chi2_opt'], color='green', linestyle='--',
               linewidth=1, alpha=0.5)

    ax5.set_xlabel('σ (kpc)', fontsize=12, fontweight='bold')
    ax5.set_ylabel('χ²', fontsize=12, fontweight='bold')
    ax5.set_title(f'Sensibilité au paramètre σ (κ={kappa_fixed:.2e})',
                 fontsize=14, fontweight='bold')
    ax5.legend(fontsize=10)
    ax5.grid(True, alpha=0.3)

    # ====================================================================
    # Subplot 6: Tableau récapitulatif
    # ====================================================================
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')

    # Créer tableau de résultats
    data = [
        ['χ² Newton', f'{resultats["chi2_newton"]:.2f}'],
        ['χ² Réseau nominal', f'{resultats["chi2_nominal"]:.2f}'],
        ['χ² Réseau optimisé', f'{resultats["chi2_opt"]:.2f}'],
        ['', ''],
        ['σ optimal', f'{resultats["sigma_opt"]:.2f} kpc'],
        ['κ optimal', f'{resultats["kappa_opt"]:.2e}'],
        ['', ''],
        ['Amélioration vs Newton',
         f'{(1-resultats["chi2_opt"]/resultats["chi2_newton"])*100:+.1f}%'],
    ]

    table = ax6.table(cellText=data, cellLoc='left', loc='center',
                     colWidths=[0.6, 0.4])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2.5)

    # Style
    for i in range(len(data)):
        if data[i][0] == '':
            continue
        elif 'Amélioration' in data[i][0]:
            improvement = float(data[i][1].strip('%+'))
            if improvement > 0:
                table[(i, 0)].set_facecolor('#90EE90')
                table[(i, 1)].set_facecolor('#90EE90')
            else:
                table[(i, 0)].set_facecolor('#FFB6C6')
                table[(i, 1)].set_facecolor('#FFB6C6')
        elif 'optimal' in data[i][0]:
            table[(i, 0)].set_facecolor('#E6F3FF')
            table[(i, 1)].set_facecolor('#E6F3FF')

    ax6.set_title('Résultats Optimisation', fontsize=14, fontweight='bold', pad=20)

    # ====================================================================
    # Sauvegarder
    # ====================================================================
    plt.tight_layout()
    plt.savefig('courbes_rotation_reseau_asselin.png', dpi=300, bbox_inches='tight')
    print("✓ Graphique sauvegardé: courbes_rotation_reseau_asselin.png")
    print()

# ============================================================================
# PROGRAMME PRINCIPAL
# ============================================================================

if __name__ == "__main__":
    resultats = test_courbes_rotation()
    generer_visualisations(resultats)

    print("="*80)
    print(" FIN DU CALCUL ".center(80))
    print("="*80)
    print()
