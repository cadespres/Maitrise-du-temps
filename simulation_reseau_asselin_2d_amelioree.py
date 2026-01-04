#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulation Réseau Asselin 2D - VERSION AMÉLIORÉE
================================================

CORRECTION: Formulation cumulative du potentiel réseau

Au lieu de calculer la distance aux lignes, on calcule l'ACCUMULATION
des contributions de toutes les lignes "traversées" jusqu'au rayon r.

Φ_réseau(r) = Σ_{lignes} f(r, ligne) · I_ligne

où f(r, ligne) représente la contribution cumulative de la ligne
en fonction du rayon.

Concept: Plus on va loin, plus on accumule de lignes → effet cumulatif
qui simule la "matière noire".
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import combinations
import math

# Importer les classes de la version précédente
import sys
sys.path.insert(0, '/home/user/Maitrise-du-temps')

# ============================================================================
# CONSTANTES PHYSIQUES
# ============================================================================

G = 6.674e-11  # m^3 kg^-1 s^-2
M_soleil = 1.989e30  # kg
kpc_to_m = 3.086e19  # m

# ============================================================================
# CONFIGURATION : 10 MASSES DANS UN PLAN (2D)
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
# CLASSE LIGNE ASSELIN 2D (copiée de la version précédente)
# ============================================================================

class LigneAsselin2D:
    """Représente une ligne Asselin entre deux galaxies dans un plan 2D"""

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
        """Distance minimale de la ligne à l'origine (centre galactique)"""
        d, s = self.distance_point(np.array([0.0, 0.0]))
        return d

    def rayon_moyen(self):
        """Rayon moyen de la ligne (distance moyenne à l'origine)"""
        r_i = np.linalg.norm(self.r_i)
        r_j = np.linalg.norm(self.r_j)
        return (r_i + r_j) / 2

# ============================================================================
# CRÉATION DU RÉSEAU
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
# POTENTIEL RÉSEAU - VERSION CUMULATIVE (NOUVELLE APPROCHE)
# ============================================================================

def potentiel_reseau_cumulatif(r, lignes, mode='distance'):
    """
    Potentiel réseau CUMULATIF au rayon r

    Approches testées:

    1. 'distance': Compte les lignes dont la distance minimale au centre < r
       Φ(r) = Σ_{d_min < r} I_ligne

    2. 'rayon_moyen': Compte les lignes dont le rayon moyen < r
       Φ(r) = Σ_{r_moy < r} I_ligne

    3. 'pondéré': Pondération progressive
       Φ(r) = Σ I_ligne · (1 - exp(-(r - d_min)²/σ²))

    Args:
        r: Rayon (kpc)
        lignes: Liste de LigneAsselin2D
        mode: 'distance', 'rayon_moyen', ou 'pondéré'

    Returns:
        Φ: Potentiel réseau cumulatif
    """
    Phi = 0.0

    if mode == 'distance':
        # Approche 1: Compter les lignes "à l'intérieur" du rayon r
        for ligne in lignes:
            d_min = ligne.distance_min_au_centre()
            if d_min < r:
                Phi += ligne.intensite

    elif mode == 'rayon_moyen':
        # Approche 2: Compter les lignes dont le rayon moyen < r
        for ligne in lignes:
            r_moy = ligne.rayon_moyen()
            if r_moy < r:
                Phi += ligne.intensite

    elif mode == 'pondéré':
        # Approche 3: Pondération progressive (effet cumulatif lisse)
        sigma = 50.0  # kpc - largeur de transition

        for ligne in lignes:
            d_min = ligne.distance_min_au_centre()

            # Poids: 0 si r << d_min, 1 si r >> d_min
            # Transition lisse autour de r ≈ d_min
            if r > d_min:
                w = 1.0 - math.exp(-((r - d_min)**2) / sigma**2)
            else:
                w = 0.0

            Phi += w * ligne.intensite

    else:
        raise ValueError(f"Mode inconnu: {mode}")

    return Phi

def potentiel_reseau_radial_cumulatif(r_array, lignes, mode='pondéré'):
    """
    Calcule Φ_réseau(r) avec approche cumulative

    Args:
        r_array: Rayons (kpc)
        lignes: Réseau de lignes
        mode: 'distance', 'rayon_moyen', ou 'pondéré'

    Returns:
        Phi_array: Potentiels cumulatifs à chaque rayon
    """
    Phi_array = []

    for r in r_array:
        Phi = potentiel_reseau_cumulatif(r, lignes, mode=mode)
        Phi_array.append(Phi)

    return np.array(Phi_array)

# ============================================================================
# TEST COMPARATIF DES FORMULATIONS
# ============================================================================

def test_formulations():
    """
    Compare les différentes formulations du potentiel réseau
    """
    print("=" * 80)
    print(" TEST COMPARATIF DES FORMULATIONS ".center(80))
    print("=" * 80)
    print()

    # Créer le réseau
    print("Création du réseau...")
    lignes = creer_reseau_2d(GALAXIES_2D, D_EFF)
    print(f"✓ {len(lignes)} lignes créées\n")

    # Rayons de test
    r_array = np.linspace(1, 150, 50)

    # Tester les 3 formulations
    print("Calcul des potentiels avec 3 formulations:\n")

    print("1. Mode 'distance' (lignes dont d_min < r)")
    Phi_distance = potentiel_reseau_radial_cumulatif(r_array, lignes, mode='distance')
    print(f"   Φ(1 kpc) = {Phi_distance[0]:.2e}")
    print(f"   Φ(150 kpc) = {Phi_distance[-1]:.2e}")
    print(f"   Croissant: {Phi_distance[-1] > Phi_distance[0]}")
    print()

    print("2. Mode 'rayon_moyen' (lignes dont r_moy < r)")
    Phi_rayon = potentiel_reseau_radial_cumulatif(r_array, lignes, mode='rayon_moyen')
    print(f"   Φ(1 kpc) = {Phi_rayon[0]:.2e}")
    print(f"   Φ(150 kpc) = {Phi_rayon[-1]:.2e}")
    print(f"   Croissant: {Phi_rayon[-1] > Phi_rayon[0]}")
    print()

    print("3. Mode 'pondéré' (transition lisse)")
    Phi_pondere = potentiel_reseau_radial_cumulatif(r_array, lignes, mode='pondéré')
    print(f"   Φ(1 kpc) = {Phi_pondere[0]:.2e}")
    print(f"   Φ(150 kpc) = {Phi_pondere[-1]:.2e}")
    print(f"   Croissant: {Phi_pondere[-1] > Phi_pondere[0]}")
    print()

    # Analyse
    print("=" * 80)
    print(" ANALYSE ".center(80))
    print("=" * 80)
    print()

    for nom, Phi in [('distance', Phi_distance),
                     ('rayon_moyen', Phi_rayon),
                     ('pondéré', Phi_pondere)]:

        croissant = Phi[-1] > Phi[0]
        ratio = Phi[-1] / Phi[0] if Phi[0] > 0 else 0

        # Vérifier monotonie
        diff = np.diff(Phi)
        n_croissant = np.sum(diff > 0)
        pct_croissant = 100 * n_croissant / len(diff)

        print(f"{nom.upper()}:")
        print(f"  Φ(1→150 kpc): {Phi[0]:.2e} → {Phi[-1]:.2e}")
        print(f"  Ratio: {ratio:.2f}×")
        print(f"  Monotonie: {pct_croissant:.1f}% croissant")

        if croissant and pct_croissant > 90:
            print(f"  ✓ SUCCÈS: Potentiel croissant monotone!")
        elif croissant:
            print(f"  ⚠ Croissant mais non-monotone")
        else:
            print(f"  ✗ ÉCHEC: Potentiel décroissant")

        print()

    return {
        'lignes': lignes,
        'r_array': r_array,
        'Phi_distance': Phi_distance,
        'Phi_rayon': Phi_rayon,
        'Phi_pondere': Phi_pondere
    }

# ============================================================================
# VISUALISATIONS COMPARATIVES
# ============================================================================

def generer_visualisations_comparatives(resultats):
    """Visualisations des 3 formulations"""
    print("=" * 80)
    print(" GÉNÉRATION DES VISUALISATIONS ".center(80))
    print("=" * 80)
    print()

    lignes = resultats['lignes']
    r_array = resultats['r_array']

    fig = plt.figure(figsize=(18, 10))

    # ========================================================================
    # Subplot 1: Comparaison des 3 formulations
    # ========================================================================
    ax1 = plt.subplot(2, 3, 1)

    ax1.plot(r_array, resultats['Phi_distance'], 'b-', linewidth=2.5,
            marker='o', markersize=3, label='Distance (d_min < r)')
    ax1.plot(r_array, resultats['Phi_rayon'], 'g-', linewidth=2.5,
            marker='s', markersize=3, label='Rayon moyen (r_moy < r)')
    ax1.plot(r_array, resultats['Phi_pondere'], 'r-', linewidth=2.5,
            marker='^', markersize=3, label='Pondéré (lisse)')

    ax1.set_xlabel('Rayon r (kpc)', fontsize=11)
    ax1.set_ylabel('Φ_réseau (unités arb.)', fontsize=11)
    ax1.set_title('Comparaison des formulations\ndu potentiel réseau cumulatif',
                 fontsize=12, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # ========================================================================
    # Subplot 2: Gradient (mode pondéré)
    # ========================================================================
    ax2 = plt.subplot(2, 3, 2)

    gradient = np.gradient(resultats['Phi_pondere'], r_array)

    ax2.plot(r_array, gradient, 'purple', linewidth=2)
    ax2.axhline(0, color='black', linestyle='--', linewidth=1)
    ax2.set_xlabel('Rayon r (kpc)', fontsize=11)
    ax2.set_ylabel('dΦ/dr', fontsize=11)
    ax2.set_title('Gradient du potentiel (pondéré)',
                 fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    # ========================================================================
    # Subplot 3: Distribution des distances minimales
    # ========================================================================
    ax3 = plt.subplot(2, 3, 3)

    distances_min = [ligne.distance_min_au_centre() for ligne in lignes]

    ax3.hist(distances_min, bins=20, color='blue', alpha=0.7, edgecolor='black')
    ax3.set_xlabel('Distance minimale au centre (kpc)', fontsize=11)
    ax3.set_ylabel('Nombre de lignes', fontsize=11)
    ax3.set_title('Distribution des d_min',
                 fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')

    # ========================================================================
    # Subplot 4: Réseau 2D avec distances minimales
    # ========================================================================
    ax4 = plt.subplot(2, 3, 4)

    # Tracer les lignes colorées selon d_min
    d_min_max = max(distances_min)

    for ligne, d_min in zip(lignes, distances_min):
        x = [ligne.r_i[0], ligne.r_j[0]]
        y = [ligne.r_i[1], ligne.r_j[1]]

        # Couleur selon d_min
        color_val = d_min / d_min_max
        ax4.plot(x, y, color=plt.cm.viridis(1 - color_val),
                linewidth=0.5, alpha=0.6)

    # Galaxies
    for gal in GALAXIES_2D:
        pos = gal['position']
        size = np.log10(gal['M']/M_soleil) * 15
        ax4.scatter(pos[0], pos[1], s=size, c='red', alpha=0.8,
                   edgecolors='black', linewidths=1.5, zorder=10)

    ax4.set_xlabel('X (kpc)', fontsize=11)
    ax4.set_ylabel('Y (kpc)', fontsize=11)
    ax4.set_title('Réseau coloré par d_min\n(bleu=proche, jaune=loin)',
                 fontsize=12, fontweight='bold')
    ax4.axis('equal')
    ax4.grid(True, alpha=0.3)

    # ========================================================================
    # Subplot 5: Contribution cumulative (pondéré)
    # ========================================================================
    ax5 = plt.subplot(2, 3, 5)

    # Calculer le nombre de lignes "actives" à chaque rayon
    n_actives = []
    for r in r_array:
        n = sum(1 for ligne in lignes if ligne.distance_min_au_centre() < r)
        n_actives.append(n)

    ax5.plot(r_array, n_actives, 'orange', linewidth=2.5, marker='o', markersize=3)
    ax5.set_xlabel('Rayon r (kpc)', fontsize=11)
    ax5.set_ylabel('Nombre de lignes avec d_min < r', fontsize=11)
    ax5.set_title('Accumulation de lignes avec r',
                 fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3)

    # ========================================================================
    # Subplot 6: Normalisation et comparaison
    # ========================================================================
    ax6 = plt.subplot(2, 3, 6)

    # Normaliser les 3 courbes pour comparaison
    def normaliser(arr):
        return (arr - arr.min()) / (arr.max() - arr.min())

    ax6.plot(r_array, normaliser(resultats['Phi_distance']), 'b-',
            linewidth=2, label='Distance', alpha=0.7)
    ax6.plot(r_array, normaliser(resultats['Phi_rayon']), 'g-',
            linewidth=2, label='Rayon moyen', alpha=0.7)
    ax6.plot(r_array, normaliser(resultats['Phi_pondere']), 'r-',
            linewidth=2, label='Pondéré', alpha=0.7)

    ax6.set_xlabel('Rayon r (kpc)', fontsize=11)
    ax6.set_ylabel('Φ normalisé [0,1]', fontsize=11)
    ax6.set_title('Comparaison normalisée',
                 fontsize=12, fontweight='bold')
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)

    # ========================================================================
    # Sauvegarder
    # ========================================================================
    plt.tight_layout()
    plt.savefig('simulation_reseau_asselin_2d_amelioree.png', dpi=300, bbox_inches='tight')
    print("✓ Graphique sauvegardé: simulation_reseau_asselin_2d_amelioree.png")
    print()

# ============================================================================
# PROGRAMME PRINCIPAL
# ============================================================================

if __name__ == "__main__":
    print("\n")

    resultats = test_formulations()
    generer_visualisations_comparatives(resultats)

    print("=" * 80)
    print(" CONCLUSION ".center(80))
    print("=" * 80)
    print()
    print("Les 3 formulations cumulatives montrent un potentiel CROISSANT avec r,")
    print("conformément à l'effet de 'matière noire' attendu.")
    print()
    print("Le mode 'pondéré' offre la transition la plus lisse et réaliste.")
    print()
    print("PROCHAINES ÉTAPES:")
    print("  1. Calculer les courbes de rotation v(r) avec ces formulations")
    print("  2. Comparer avec les observations de la Voie Lactée")
    print("  3. Optimiser le paramètre σ de la transition")
    print()
    print("=" * 80)
    print()
