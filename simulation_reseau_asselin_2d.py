#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulation Simplifiée 2D du Réseau Asselin
==========================================

Étape 1 : Simulation Simplifiée 2D

Configuration :
- 10 masses ponctuelles dans un plan
- Positions et masses réalistes (Voie Lactée + voisines)
- d_eff = 100 kpc (fixé)

Calcul :
- Toutes les lignes d'ordre 1 (45 lignes)
- Toutes les intersections (990 max)
- Filtrer intersections réelles (critère géométrique)
- Calculer Φ_réseau(r) le long du disque

Vérification :
- Φ_réseau(r) augmente-t-il avec r ?
- Profil cohérent avec observations ?
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import combinations
import math

# ============================================================================
# CONSTANTES PHYSIQUES
# ============================================================================

G = 6.674e-11  # m^3 kg^-1 s^-2
M_soleil = 1.989e30  # kg
kpc_to_m = 3.086e19  # m

# ============================================================================
# CONFIGURATION : 10 MASSES DANS UN PLAN (2D)
# ============================================================================

# Projection des galaxies du Groupe Local dans le plan XY
GALAXIES_2D = [
    {
        'nom': 'Voie Lactée',
        'M': 8.0e10 * M_soleil,
        'position': np.array([0.0, 0.0])  # Centre du référentiel
    },
    {
        'nom': 'M31 (Andromède)',
        'M': 1.5e12 * M_soleil,
        'position': np.array([750.0, 250.0])  # ~780 kpc
    },
    {
        'nom': 'M33 (Triangulum)',
        'M': 4.0e10 * M_soleil,
        'position': np.array([840.0, 120.0])  # ~850 kpc
    },
    {
        'nom': 'Grand Nuage Magellan (LMC)',
        'M': 1.0e10 * M_soleil,
        'position': np.array([-40.0, 30.0])  # ~50 kpc
    },
    {
        'nom': 'Petit Nuage Magellan (SMC)',
        'M': 7.0e9 * M_soleil,
        'position': np.array([-50.0, 40.0])  # ~63 kpc
    },
    {
        'nom': 'Naine du Sagittaire',
        'M': 4.0e8 * M_soleil,
        'position': np.array([20.0, -15.0])  # ~26 kpc
    },
    {
        'nom': 'Naine du Sculpteur',
        'M': 2.0e8 * M_soleil,
        'position': np.array([70.0, -50.0])  # ~86 kpc
    },
    {
        'nom': 'Naine du Fourneau',
        'M': 2.0e8 * M_soleil,
        'position': np.array([-120.0, 80.0])  # ~147 kpc
    },
    {
        'nom': 'Naine de la Carène',
        'M': 1.5e8 * M_soleil,
        'position': np.array([60.0, -70.0])  # ~94 kpc
    },
    {
        'nom': 'Naine du Dragon',
        'M': 1.0e8 * M_soleil,
        'position': np.array([70.0, 50.0])  # ~92 kpc
    }
]

# Paramètre fixe
D_EFF = 100.0  # kpc

# ============================================================================
# CLASSE LIGNE ASSELIN 2D
# ============================================================================

class LigneAsselin2D:
    """Représente une ligne Asselin entre deux galaxies dans un plan 2D"""

    def __init__(self, galaxie_i, galaxie_j, idx_i, idx_j, d_eff=D_EFF):
        """
        Args:
            galaxie_i, galaxie_j: Dictionnaires avec 'M' et 'position'
            idx_i, idx_j: Indices des galaxies
            d_eff: Distance effective (kpc)
        """
        self.i = galaxie_i
        self.j = galaxie_j
        self.idx_i = idx_i
        self.idx_j = idx_j
        self.r_i = galaxie_i['position']  # 2D
        self.r_j = galaxie_j['position']  # 2D

        # Distance entre galaxies
        self.d_ij = np.linalg.norm(self.r_j - self.r_i)

        # Vecteur direction normalisé
        if self.d_ij > 0:
            self.direction = (self.r_j - self.r_i) / self.d_ij
        else:
            self.direction = np.array([1.0, 0.0])

        # Intensité Asselin
        self.intensite = self.calculer_intensite(d_eff)

    def calculer_intensite(self, d_eff):
        """
        Intensité Liaison Asselin

        I_ij = √(M_i·M_j) / d²_ij · exp(-d_ij/d_eff)
        """
        M_i = self.i['M']
        M_j = self.j['M']

        if self.d_ij < 0.1:  # Protection
            d = 0.1
        else:
            d = self.d_ij

        I = math.sqrt(M_i * M_j) / (d**2) * math.exp(-d / d_eff)
        return I

    def distance_point(self, P):
        """
        Distance d'un point P à la ligne (segment)

        Args:
            P: Position 2D (array)

        Returns:
            d: Distance minimale de P à la ligne (kpc)
            s: Paramètre de projection [0,1]
        """
        # Vecteur direction
        u = self.r_j - self.r_i

        # Vecteur vers le point
        w = P - self.r_i

        # Projection
        u_dot_u = np.dot(u, u)
        if u_dot_u < 1e-10:
            s = 0.0
        else:
            s = np.dot(w, u) / u_dot_u

        # Clamp à [0,1] (limiter au segment)
        s = max(0.0, min(s, 1.0))

        # Point le plus proche sur la ligne
        P_proj = self.r_i + s * u

        # Distance
        d = np.linalg.norm(P - P_proj)

        return d, s

    def __repr__(self):
        return f"Ligne({self.i['nom']} - {self.j['nom']}, d={self.d_ij:.1f}kpc, I={self.intensite:.2e})"

# ============================================================================
# CALCUL D'INTERSECTION DE DEUX LIGNES 2D
# ============================================================================

def intersection_lignes_2d(ligne1, ligne2):
    """
    Calcule l'intersection de deux segments de ligne dans le plan 2D

    Ligne 1: P1 = r1_i + s*(r1_j - r1_i),  s ∈ [0,1]
    Ligne 2: P2 = r2_i + t*(r2_j - r2_i),  t ∈ [0,1]

    Args:
        ligne1, ligne2: LigneAsselin2D

    Returns:
        dict avec:
            - 'existe': bool, True si intersection réelle
            - 'point': np.array([x,y]) si existe
            - 's', 't': paramètres si existe
            - 'distance': distance minimale entre les lignes
    """
    # Points et directions
    P1 = ligne1.r_i
    d1 = ligne1.r_j - ligne1.r_i

    P2 = ligne2.r_i
    d2 = ligne2.r_j - ligne2.r_i

    # Système linéaire: P1 + s*d1 = P2 + t*d2
    # => s*d1 - t*d2 = P2 - P1
    # => [d1 | -d2] [s, t]^T = P2 - P1

    # Matrice 2x2
    A = np.column_stack([d1, -d2])
    b = P2 - P1

    # Déterminant
    det = A[0,0]*A[1,1] - A[0,1]*A[1,0]

    # Si déterminant nul -> lignes parallèles
    if abs(det) < 1e-10:
        # Lignes parallèles, pas d'intersection
        # Calculer distance minimale
        # Distance d'un point de ligne2 à ligne1
        d_min, _ = ligne1.distance_point(P2)

        return {
            'existe': False,
            'point': None,
            's': None,
            't': None,
            'distance': d_min,
            'raison': 'parallèles'
        }

    # Résoudre le système
    s = (A[1,1]*b[0] - A[0,1]*b[1]) / det
    t = (A[0,0]*b[1] - A[1,0]*b[0]) / det

    # Vérifier si s et t sont dans [0,1]
    # (intersection dans les segments, pas juste les lignes infinies)
    if 0 <= s <= 1 and 0 <= t <= 1:
        # Intersection réelle !
        point = P1 + s * d1

        return {
            'existe': True,
            'point': point,
            's': s,
            't': t,
            'distance': 0.0,
            'raison': 'intersection'
        }
    else:
        # Intersection hors des segments
        # Calculer distance minimale entre segments

        # Tester 4 combinaisons: extrémités de chaque segment
        distances = [
            ligne1.distance_point(ligne2.r_i)[0],
            ligne1.distance_point(ligne2.r_j)[0],
            ligne2.distance_point(ligne1.r_i)[0],
            ligne2.distance_point(ligne1.r_j)[0]
        ]
        d_min = min(distances)

        return {
            'existe': False,
            'point': None,
            's': s,
            't': t,
            'distance': d_min,
            'raison': f'hors_segment (s={s:.2f}, t={t:.2f})'
        }

# ============================================================================
# CRÉATION DU RÉSEAU
# ============================================================================

def creer_reseau_2d(galaxies, d_eff=D_EFF):
    """
    Crée toutes les lignes Asselin d'ordre 1 entre paires de galaxies

    Args:
        galaxies: Liste de dictionnaires galaxies
        d_eff: Distance effective (kpc)

    Returns:
        lignes: Liste de LigneAsselin2D
    """
    lignes = []
    N = len(galaxies)

    for i in range(N):
        for j in range(i+1, N):
            ligne = LigneAsselin2D(galaxies[i], galaxies[j], i, j, d_eff)
            lignes.append(ligne)

    print(f"✓ Réseau créé: {len(lignes)} lignes d'ordre 1")
    print(f"  Nombre attendu: C({N},2) = {N*(N-1)//2}")

    return lignes

# ============================================================================
# CALCUL DES INTERSECTIONS
# ============================================================================

def calculer_intersections(lignes):
    """
    Calcule toutes les intersections entre paires de lignes

    Args:
        lignes: Liste de LigneAsselin2D

    Returns:
        dict avec:
            - 'toutes': liste de toutes les combinaisons testées
            - 'reelles': liste des intersections réelles (dans les segments)
            - 'hors_segment': intersections des droites infinies mais hors segments
            - 'paralleles': paires de lignes parallèles
    """
    N_lignes = len(lignes)
    N_combinaisons = N_lignes * (N_lignes - 1) // 2

    print(f"\nCalcul des intersections...")
    print(f"  Nombre de paires de lignes: C({N_lignes},2) = {N_combinaisons}")

    toutes = []
    reelles = []
    hors_segment = []
    paralleles = []

    # Tester toutes les paires
    for i, ligne1 in enumerate(lignes):
        for j in range(i+1, len(lignes)):
            ligne2 = lignes[j]

            # Éviter de calculer intersection entre lignes qui partagent une galaxie
            # (elles se touchent forcément à cette galaxie)
            if ligne1.idx_i in [ligne2.idx_i, ligne2.idx_j] or \
               ligne1.idx_j in [ligne2.idx_i, ligne2.idx_j]:
                # Intersection triviale aux extrémités
                continue

            inter = intersection_lignes_2d(ligne1, ligne2)
            inter['ligne1'] = ligne1
            inter['ligne2'] = ligne2

            toutes.append(inter)

            if inter['existe']:
                reelles.append(inter)
            elif inter['raison'] == 'parallèles':
                paralleles.append(inter)
            else:
                hors_segment.append(inter)

    print(f"\n  Résultats:")
    print(f"    - Combinaisons testées: {len(toutes)}")
    print(f"    - Intersections réelles: {len(reelles)}")
    print(f"    - Hors segment: {len(hors_segment)}")
    print(f"    - Parallèles: {len(paralleles)}")

    return {
        'toutes': toutes,
        'reelles': reelles,
        'hors_segment': hors_segment,
        'paralleles': paralleles
    }

# ============================================================================
# POTENTIEL RÉSEAU
# ============================================================================

def potentiel_reseau_2d(P, lignes, sigma=10.0):
    """
    Potentiel au point P depuis le réseau de lignes Asselin

    Φ_réseau(P) = Σ_lignes w(d_ligne) · I_ligne

    où w(d) = exp(-d²/σ²) (poids gaussien)

    Args:
        P: Position 2D (kpc)
        lignes: Liste de LigneAsselin2D
        sigma: Largeur gaussienne (kpc)

    Returns:
        Φ en unités arbitraires (contribution réseau)
    """
    Phi = 0.0

    for ligne in lignes:
        # Distance du point à la ligne
        d_ligne, s = ligne.distance_point(P)

        # Poids gaussien
        w = math.exp(-d_ligne**2 / sigma**2)

        # Contribution au potentiel
        Phi += w * ligne.intensite

    return Phi

def potentiel_reseau_radial(r_array, lignes, sigma=10.0, n_angles=8):
    """
    Calcule Φ_réseau(r) en moyennant sur plusieurs angles

    Args:
        r_array: Rayons (kpc)
        lignes: Réseau de lignes
        sigma: Largeur gaussienne (kpc)
        n_angles: Nombre d'angles pour moyenner

    Returns:
        Phi_array: Potentiels moyens à chaque rayon
    """
    Phi_array = []

    for r in r_array:
        Phi_moy = 0.0

        # Moyenner sur plusieurs angles pour avoir un profil radial
        for i in range(n_angles):
            theta = 2 * np.pi * i / n_angles
            P = np.array([r * np.cos(theta), r * np.sin(theta)])
            Phi_moy += potentiel_reseau_2d(P, lignes, sigma)

        Phi_moy /= n_angles
        Phi_array.append(Phi_moy)

    return np.array(Phi_array)

# ============================================================================
# POTENTIEL AUX INTERSECTIONS
# ============================================================================

def potentiel_aux_intersections(intersections_reelles, lignes, sigma=10.0):
    """
    Calcule le potentiel réseau aux points d'intersection

    Args:
        intersections_reelles: Liste des intersections
        lignes: Réseau de lignes
        sigma: Largeur gaussienne

    Returns:
        liste de (point, Phi)
    """
    resultats = []

    for inter in intersections_reelles:
        P = inter['point']
        Phi = potentiel_reseau_2d(P, lignes, sigma)
        resultats.append((P, Phi))

    return resultats

# ============================================================================
# TEST PRINCIPAL
# ============================================================================

def test_simulation_2d():
    """
    Test complet de la simulation 2D du réseau Asselin
    """
    print("=" * 80)
    print(" SIMULATION 2D DU RÉSEAU ASSELIN ".center(80))
    print("=" * 80)
    print()

    # Afficher configuration
    print("CONFIGURATION:")
    print(f"  Nombre de galaxies: {len(GALAXIES_2D)}")
    print(f"  Distance effective: d_eff = {D_EFF} kpc")
    print()

    print("GALAXIES:")
    for i, gal in enumerate(GALAXIES_2D):
        pos = gal['position']
        d = np.linalg.norm(pos)
        print(f"  {i}. {gal['nom']:<30} M={gal['M']/M_soleil:.2e} M☉, "
              f"pos=({pos[0]:>7.1f}, {pos[1]:>7.1f}) kpc, d={d:.1f} kpc")
    print()

    # Étape 1: Créer le réseau
    print("ÉTAPE 1: Création du réseau de lignes d'ordre 1")
    print("-" * 80)
    lignes = creer_reseau_2d(GALAXIES_2D, D_EFF)
    print()

    # Statistiques sur les lignes
    intensites = [l.intensite for l in lignes]
    distances = [l.d_ij for l in lignes]

    print("Statistiques des lignes:")
    print(f"  Intensité min: {min(intensites):.2e}")
    print(f"  Intensité max: {max(intensites):.2e}")
    print(f"  Intensité moy: {np.mean(intensites):.2e}")
    print(f"  Distance min: {min(distances):.1f} kpc")
    print(f"  Distance max: {max(distances):.1f} kpc")
    print(f"  Distance moy: {np.mean(distances):.1f} kpc")
    print()

    # Étape 2: Calculer les intersections
    print("ÉTAPE 2: Calcul des intersections")
    print("-" * 80)
    intersections = calculer_intersections(lignes)
    print()

    # Étape 3: Calculer Φ_réseau(r)
    print("ÉTAPE 3: Calcul de Φ_réseau(r) le long du disque")
    print("-" * 80)

    r_array = np.linspace(1, 150, 50)
    Phi_reseau = potentiel_reseau_radial(r_array, lignes, sigma=10.0, n_angles=16)

    print(f"  Calculé pour {len(r_array)} rayons de 1 à 150 kpc")
    print(f"  Φ_réseau(r=1 kpc) = {Phi_reseau[0]:.2e}")
    print(f"  Φ_réseau(r=50 kpc) = {Phi_reseau[24]:.2e}")
    print(f"  Φ_réseau(r=150 kpc) = {Phi_reseau[-1]:.2e}")
    print()

    # Étape 4: Vérifications
    print("ÉTAPE 4: Vérifications")
    print("-" * 80)

    # Vérifier si Φ augmente avec r
    differences = np.diff(Phi_reseau)
    n_croissant = np.sum(differences > 0)
    n_decroissant = np.sum(differences < 0)

    print(f"  Évolution de Φ_réseau(r):")
    print(f"    - Segments croissants: {n_croissant}/{len(differences)} ({100*n_croissant/len(differences):.1f}%)")
    print(f"    - Segments décroissants: {n_decroissant}/{len(differences)} ({100*n_decroissant/len(differences):.1f}%)")

    if Phi_reseau[-1] > Phi_reseau[0]:
        print(f"  ✓ Φ_réseau augmente globalement: {Phi_reseau[0]:.2e} → {Phi_reseau[-1]:.2e}")
        ratio = Phi_reseau[-1] / Phi_reseau[0]
        print(f"    Ratio Φ(r_max)/Φ(r_min) = {ratio:.2f}")
    else:
        print(f"  ✗ Φ_réseau diminue: {Phi_reseau[0]:.2e} → {Phi_reseau[-1]:.2e}")

    print()

    # Potentiel aux intersections
    if len(intersections['reelles']) > 0:
        print(f"  Potentiel aux {len(intersections['reelles'])} intersections réelles:")
        potentiels_inter = potentiel_aux_intersections(intersections['reelles'], lignes, sigma=10.0)

        Phi_values = [p[1] for p in potentiels_inter]
        print(f"    - Φ min: {min(Phi_values):.2e}")
        print(f"    - Φ max: {max(Phi_values):.2e}")
        print(f"    - Φ moy: {np.mean(Phi_values):.2e}")

    print()

    return {
        'lignes': lignes,
        'intersections': intersections,
        'r_array': r_array,
        'Phi_reseau': Phi_reseau,
        'galaxies': GALAXIES_2D
    }

# ============================================================================
# VISUALISATIONS
# ============================================================================

def generer_visualisations(resultats):
    """
    Génère les visualisations de la simulation 2D
    """
    print("=" * 80)
    print(" GÉNÉRATION DES VISUALISATIONS ".center(80))
    print("=" * 80)
    print()

    lignes = resultats['lignes']
    intersections = resultats['intersections']
    r_array = resultats['r_array']
    Phi_reseau = resultats['Phi_reseau']
    galaxies = resultats['galaxies']

    fig = plt.figure(figsize=(18, 12))

    # ========================================================================
    # Subplot 1: Vue d'ensemble du réseau 2D
    # ========================================================================
    ax1 = plt.subplot(2, 3, 1)

    # Tracer toutes les lignes
    intensites = [l.intensite for l in lignes]
    I_max = max(intensites)

    for ligne in lignes:
        x = [ligne.r_i[0], ligne.r_j[0]]
        y = [ligne.r_i[1], ligne.r_j[1]]
        alpha = min(1.0, (ligne.intensite / I_max) * 5)
        linewidth = 0.3 + (ligne.intensite / I_max) * 2
        ax1.plot(x, y, 'b-', alpha=alpha, linewidth=linewidth)

    # Tracer les galaxies
    for i, gal in enumerate(galaxies):
        pos = gal['position']
        size = np.log10(gal['M']/M_soleil) * 15
        ax1.scatter(pos[0], pos[1], s=size, c='red', alpha=0.8,
                   edgecolors='black', linewidths=1.5, zorder=10)
        ax1.text(pos[0], pos[1]-20, str(i), fontsize=8, ha='center',
                fontweight='bold')

    # Tracer les intersections réelles
    if len(intersections['reelles']) > 0:
        pts_inter = [inter['point'] for inter in intersections['reelles']]
        pts_inter = np.array(pts_inter)
        ax1.scatter(pts_inter[:,0], pts_inter[:,1], s=30, c='green',
                   marker='x', linewidths=2, zorder=15, label='Intersections')

    ax1.set_xlabel('X (kpc)', fontsize=11)
    ax1.set_ylabel('Y (kpc)', fontsize=11)
    ax1.set_title(f'Réseau Asselin 2D\n{len(lignes)} lignes, {len(intersections["reelles"])} intersections',
                 fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axis('equal')
    if len(intersections['reelles']) > 0:
        ax1.legend(fontsize=9)

    # ========================================================================
    # Subplot 2: Zoom sur le centre
    # ========================================================================
    ax2 = plt.subplot(2, 3, 2)

    # Même chose mais zoom sur [-200, 200] kpc
    for ligne in lignes:
        x = [ligne.r_i[0], ligne.r_j[0]]
        y = [ligne.r_i[1], ligne.r_j[1]]
        alpha = min(1.0, (ligne.intensite / I_max) * 5)
        linewidth = 0.5 + (ligne.intensite / I_max) * 2
        ax2.plot(x, y, 'b-', alpha=alpha, linewidth=linewidth)

    for i, gal in enumerate(galaxies):
        pos = gal['position']
        size = np.log10(gal['M']/M_soleil) * 20
        ax2.scatter(pos[0], pos[1], s=size, c='red', alpha=0.8,
                   edgecolors='black', linewidths=1.5, zorder=10)
        ax2.text(pos[0]+5, pos[1]+5, gal['nom'].split()[0],
                fontsize=7, ha='left')

    if len(intersections['reelles']) > 0:
        pts_inter = np.array([inter['point'] for inter in intersections['reelles']])
        ax2.scatter(pts_inter[:,0], pts_inter[:,1], s=50, c='green',
                   marker='x', linewidths=2, zorder=15)

    ax2.set_xlim(-200, 200)
    ax2.set_ylim(-200, 200)
    ax2.set_xlabel('X (kpc)', fontsize=11)
    ax2.set_ylabel('Y (kpc)', fontsize=11)
    ax2.set_title('Zoom sur le centre (±200 kpc)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.axis('equal')

    # ========================================================================
    # Subplot 3: Distribution des intensités
    # ========================================================================
    ax3 = plt.subplot(2, 3, 3)

    ax3.hist(intensites, bins=20, color='blue', alpha=0.7, edgecolor='black')
    ax3.set_xlabel('Intensité Asselin', fontsize=11)
    ax3.set_ylabel('Nombre de lignes', fontsize=11)
    ax3.set_title('Distribution des intensités des lignes', fontsize=12, fontweight='bold')
    ax3.set_yscale('log')
    ax3.grid(True, alpha=0.3)

    # ========================================================================
    # Subplot 4: Φ_réseau(r) - Profil radial
    # ========================================================================
    ax4 = plt.subplot(2, 3, 4)

    ax4.plot(r_array, Phi_reseau, 'purple', linewidth=2.5, marker='o',
            markersize=4, label='Φ_réseau(r)')
    ax4.set_xlabel('Rayon r (kpc)', fontsize=11)
    ax4.set_ylabel('Φ_réseau (unités arbitraires)', fontsize=11)
    ax4.set_title('Potentiel réseau radial Φ_réseau(r)', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.legend(fontsize=10)

    # Ajouter une ligne de tendance
    z = np.polyfit(r_array, Phi_reseau, 2)
    p = np.poly1d(z)
    ax4.plot(r_array, p(r_array), 'r--', linewidth=1.5, alpha=0.7,
            label=f'Tendance (degré 2)')
    ax4.legend(fontsize=9)

    # ========================================================================
    # Subplot 5: Gradient de Φ_réseau(r)
    # ========================================================================
    ax5 = plt.subplot(2, 3, 5)

    # Calculer le gradient (dérivée numérique)
    gradient = np.gradient(Phi_reseau, r_array)

    ax5.plot(r_array, gradient, 'green', linewidth=2, marker='s', markersize=3)
    ax5.axhline(0, color='black', linestyle='--', linewidth=1)
    ax5.set_xlabel('Rayon r (kpc)', fontsize=11)
    ax5.set_ylabel('dΦ/dr (unités arbitraires)', fontsize=11)
    ax5.set_title('Gradient radial du potentiel réseau', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3)

    # Zones où le gradient est positif
    positif = gradient > 0
    if np.any(positif):
        pct_positif = 100 * np.sum(positif) / len(positif)
        ax5.text(0.05, 0.95, f'{pct_positif:.1f}% croissant',
                transform=ax5.transAxes, fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # ========================================================================
    # Subplot 6: Statistiques des intersections
    # ========================================================================
    ax6 = plt.subplot(2, 3, 6)

    # Créer un tableau de statistiques
    stats_data = [
        ['Lignes d\'ordre 1', len(lignes)],
        ['Paires testées', len(intersections['toutes'])],
        ['Intersections réelles', len(intersections['reelles'])],
        ['Hors segment', len(intersections['hors_segment'])],
        ['Parallèles', len(intersections['paralleles'])],
    ]

    ax6.axis('off')
    table = ax6.table(cellText=stats_data,
                     colLabels=['Catégorie', 'Nombre'],
                     cellLoc='left',
                     loc='center',
                     colWidths=[0.7, 0.3])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    # Style
    for i in range(len(stats_data) + 1):
        if i == 0:
            table[(i, 0)].set_facecolor('#4CAF50')
            table[(i, 1)].set_facecolor('#4CAF50')
            table[(i, 0)].set_text_props(weight='bold', color='white')
            table[(i, 1)].set_text_props(weight='bold', color='white')
        else:
            if i == 3:  # Intersections réelles
                table[(i, 0)].set_facecolor('#90EE90')
                table[(i, 1)].set_facecolor('#90EE90')

    ax6.set_title('Statistiques du réseau', fontsize=12, fontweight='bold', pad=20)

    # ========================================================================
    # Sauvegarder
    # ========================================================================
    plt.tight_layout()
    plt.savefig('simulation_reseau_asselin_2d.png', dpi=300, bbox_inches='tight')
    print("✓ Graphique sauvegardé: simulation_reseau_asselin_2d.png")
    print()

# ============================================================================
# RAPPORT DÉTAILLÉ
# ============================================================================

def generer_rapport(resultats):
    """
    Génère un rapport textuel détaillé
    """
    print("=" * 80)
    print(" RAPPORT DÉTAILLÉ ".center(80))
    print("=" * 80)
    print()

    lignes = resultats['lignes']
    intersections = resultats['intersections']
    Phi_reseau = resultats['Phi_reseau']

    print("1. RÉSEAU DE LIGNES")
    print("-" * 80)
    print(f"   Nombre total de lignes: {len(lignes)}")
    print()
    print("   Top 5 lignes par intensité:")
    lignes_triees = sorted(lignes, key=lambda l: l.intensite, reverse=True)
    for i, ligne in enumerate(lignes_triees[:5]):
        print(f"     {i+1}. {ligne.i['nom']} - {ligne.j['nom']}")
        print(f"        d={ligne.d_ij:.1f} kpc, I={ligne.intensite:.2e}")
    print()

    print("2. INTERSECTIONS")
    print("-" * 80)
    print(f"   Intersections réelles: {len(intersections['reelles'])}")

    if len(intersections['reelles']) > 0:
        print()
        print("   Exemples d'intersections réelles:")
        for i, inter in enumerate(intersections['reelles'][:5]):
            l1 = inter['ligne1']
            l2 = inter['ligne2']
            pt = inter['point']
            print(f"     {i+1}. Ligne ({l1.i['nom']} - {l1.j['nom']})")
            print(f"        ∩ ({l2.i['nom']} - {l2.j['nom']})")
            print(f"        Point: ({pt[0]:.1f}, {pt[1]:.1f}) kpc")
    print()

    print("3. POTENTIEL RÉSEAU Φ_réseau(r)")
    print("-" * 80)
    print(f"   Φ(r=1 kpc) = {Phi_reseau[0]:.2e}")
    print(f"   Φ(r=50 kpc) = {Phi_reseau[len(Phi_reseau)//3]:.2e}")
    print(f"   Φ(r=100 kpc) = {Phi_reseau[2*len(Phi_reseau)//3]:.2e}")
    print(f"   Φ(r=150 kpc) = {Phi_reseau[-1]:.2e}")
    print()

    # Analyse de la croissance
    if Phi_reseau[-1] > Phi_reseau[0]:
        print("   ✓ Comportement attendu: Φ_réseau(r) AUGMENTE avec r")
        print(f"     Facteur de croissance: {Phi_reseau[-1]/Phi_reseau[0]:.2f}×")
        print()
        print("   INTERPRÉTATION:")
        print("     • Plus on s'éloigne du centre, plus on accumule de")
        print("       contributions des lignes Asselin")
        print("     • Cohérent avec l'effet de 'masse manquante' observé")
        print("     • Le réseau crée un potentiel effectif croissant")
    else:
        print("   ✗ Comportement inattendu: Φ_réseau(r) DIMINUE avec r")
        print()
        print("   À INVESTIGUER:")
        print("     • Vérifier le poids gaussien σ")
        print("     • Analyser la distribution des lignes")

    print()

    print("4. CONCLUSION")
    print("-" * 80)
    print("   La simulation 2D du réseau Asselin démontre:")
    print()
    print(f"   • {len(lignes)} lignes gravitationnelles créées entre 10 galaxies")
    print(f"   • {len(intersections['reelles'])} intersections réelles détectées")
    print(f"   • Potentiel réseau Φ_réseau(r) calculé sur 50 rayons")

    if Phi_reseau[-1] > Phi_reseau[0]:
        print()
        print("   ✓ SUCCÈS: Le potentiel réseau croît avec la distance,")
        print("     reproduisant l'effet de 'matière noire' sans masse cachée!")

    print()
    print("=" * 80)
    print()

# ============================================================================
# PROGRAMME PRINCIPAL
# ============================================================================

if __name__ == "__main__":
    print("\n")

    # Exécuter la simulation
    resultats = test_simulation_2d()

    # Générer les visualisations
    generer_visualisations(resultats)

    # Générer le rapport
    generer_rapport(resultats)

    print("=" * 80)
    print(" FIN DE LA SIMULATION ".center(80))
    print("=" * 80)
    print()

    print("Fichiers générés:")
    print("  • simulation_reseau_asselin_2d.png - Visualisations complètes")
    print()
