#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Intersections d'Ordre 2+ et Renforcement Non-Linéaire
=====================================================

Objectif:
- Identifier les points où 3+ lignes se croisent (intersections d'ordre 2+)
- Calculer le renforcement non-linéaire à ces points
- Évaluer l'impact sur le potentiel réseau total

Concept:
- Intersection ordre 1: 2 lignes se croisent
- Intersection ordre 2: 3 lignes se croisent
- Intersection ordre n: n lignes se croisent

Hypothèse:
Le renforcement aux intersections d'ordre élevé est NON-LINÉAIRE
Φ_inter(n) ∝ n^α où α > 1 (amplification)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import combinations
import math

# ============================================================================
# CONSTANTES
# ============================================================================

G = 6.674e-11
M_soleil = 1.989e30
kpc_to_m = 3.086e19

# ============================================================================
# GALAXIES 3D
# ============================================================================

GALAXIES_3D = [
    {'nom': 'Voie Lactée', 'M': 8.0e10 * M_soleil, 'position': np.array([0.0, 0.0, 0.0])},
    {'nom': 'M31 (Andromède)', 'M': 1.5e12 * M_soleil, 'position': np.array([750.0, 250.0, 100.0])},
    {'nom': 'M33 (Triangulum)', 'M': 4.0e10 * M_soleil, 'position': np.array([840.0, 120.0, -50.0])},
    {'nom': 'Grand Nuage Magellan (LMC)', 'M': 1.0e10 * M_soleil, 'position': np.array([-40.0, 30.0, -20.0])},
    {'nom': 'Petit Nuage Magellan (SMC)', 'M': 7.0e9 * M_soleil, 'position': np.array([-50.0, 40.0, -15.0])},
    {'nom': 'Naine du Sagittaire', 'M': 4.0e8 * M_soleil, 'position': np.array([20.0, -15.0, 10.0])},
    {'nom': 'Naine du Sculpteur', 'M': 2.0e8 * M_soleil, 'position': np.array([70.0, -50.0, 30.0])},
    {'nom': 'Naine du Fourneau', 'M': 2.0e8 * M_soleil, 'position': np.array([-120.0, 80.0, -40.0])},
    {'nom': 'Naine de la Carène', 'M': 1.5e8 * M_soleil, 'position': np.array([60.0, -70.0, 25.0])},
    {'nom': 'Naine du Dragon', 'M': 1.0e8 * M_soleil, 'position': np.array([70.0, 50.0, -40.0])}
]

D_EFF = 100.0

# ============================================================================
# CLASSE LIGNE 3D
# ============================================================================

class LigneAsselin3D:
    """Ligne Asselin 3D"""

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
            self.direction = np.array([1.0, 0.0, 0.0])

        self.intensite = self.calculer_intensite(d_eff)

    def calculer_intensite(self, d_eff):
        M_i = self.i['M']
        M_j = self.j['M']
        d = max(self.d_ij, 0.1)
        I = math.sqrt(M_i * M_j) / (d**2) * math.exp(-d / d_eff)
        return I

    def distance_point(self, P):
        """Distance d'un point P à la ligne"""
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

    def point_le_plus_proche(self, autre_ligne):
        """
        Trouve le point le plus proche entre deux lignes 3D

        Returns:
            (P1, P2, distance): Points les plus proches et leur distance
        """
        # Lignes paramétriques:
        # L1(s) = r1_i + s * (r1_j - r1_i)
        # L2(t) = r2_i + t * (r2_j - r2_i)

        p1 = self.r_i
        d1 = self.r_j - self.r_i

        p2 = autre_ligne.r_i
        d2 = autre_ligne.r_j - autre_ligne.r_i

        # Vecteur entre les deux origins
        w = p1 - p2

        # Coefficients
        a = np.dot(d1, d1)
        b = np.dot(d1, d2)
        c = np.dot(d2, d2)
        d = np.dot(d1, w)
        e = np.dot(d2, w)

        denom = a * c - b * b

        if abs(denom) < 1e-10:
            # Lignes parallèles
            s = 0.0
            t = e / c if abs(c) > 1e-10 else 0.0
        else:
            s = (b * e - c * d) / denom
            t = (a * e - b * d) / denom

        # Clamp à [0,1]
        s = max(0.0, min(1.0, s))
        t = max(0.0, min(1.0, t))

        # Points les plus proches
        P1 = p1 + s * d1
        P2 = p2 + t * d2

        distance = np.linalg.norm(P1 - P2)

        return P1, P2, distance

    def __repr__(self):
        return f"L{self.idx_i}-{self.idx_j}"

# ============================================================================
# RÉSEAU 3D
# ============================================================================

def creer_reseau_3d(galaxies, d_eff=D_EFF):
    lignes = []
    N = len(galaxies)
    for i in range(N):
        for j in range(i+1, N):
            ligne = LigneAsselin3D(galaxies[i], galaxies[j], i, j, d_eff)
            lignes.append(ligne)
    return lignes

# ============================================================================
# DÉTECTION DES INTERSECTIONS D'ORDRE 2+
# ============================================================================

def trouver_point_convergence(lignes_list, tolerance=5.0):
    """
    Trouve le point de convergence approximatif d'un groupe de lignes

    Args:
        lignes_list: Liste de lignes
        tolerance: Distance maximale pour considérer une convergence (kpc)

    Returns:
        (point, convergence_ok): Point moyen et bool indiquant si convergence
    """
    # Calculer le barycentre des points les plus proches par paires
    points = []

    for i, ligne1 in enumerate(lignes_list):
        for ligne2 in lignes_list[i+1:]:
            P1, P2, dist = ligne1.point_le_plus_proche(ligne2)

            if dist < tolerance:
                # Prendre le milieu
                point_inter = (P1 + P2) / 2
                points.append(point_inter)

    if len(points) == 0:
        return None, False

    # Point moyen
    point_convergence = np.mean(points, axis=0)

    # Vérifier que toutes les lignes passent près de ce point
    distances = []
    for ligne in lignes_list:
        d, s = ligne.distance_point(point_convergence)
        distances.append(d)

    # Toutes les distances doivent être < tolerance
    convergence_ok = all(d < tolerance for d in distances)

    return point_convergence, convergence_ok

def detecter_intersections_ordre2_plus(lignes, tolerance=10.0):
    """
    Détecte toutes les intersections d'ordre 2 et plus

    Args:
        lignes: Liste de lignes 3D
        tolerance: Distance max pour considérer intersection (kpc)

    Returns:
        dict: {ordre: [(point, lignes, intensité_totale), ...]}
    """
    print("\n" + "="*80)
    print(" DÉTECTION DES INTERSECTIONS D'ORDRE 2+ ".center(80))
    print("="*80)
    print()
    print(f"Tolérance: {tolerance} kpc")
    print()

    intersections_par_ordre = {}

    # Tester toutes les combinaisons de 3+ lignes
    N = len(lignes)

    print(f"Test des combinaisons de lignes...")

    # Ordre 2: 3 lignes
    print(f"\n  Ordre 2 (3 lignes): C({N},3) = {len(list(combinations(range(N), 3)))} combinaisons")

    count_ordre2 = 0
    intersections_ordre2 = []

    for combo in combinations(lignes, 3):
        point, ok = trouver_point_convergence(combo, tolerance)

        if ok and point is not None:
            # Calculer intensité totale
            intensite_tot = sum(l.intensite for l in combo)
            intersections_ordre2.append((point, combo, intensite_tot))
            count_ordre2 += 1

    intersections_par_ordre[2] = intersections_ordre2
    print(f"    Trouvé: {count_ordre2} intersections")

    # Ordre 3: 4 lignes
    print(f"\n  Ordre 3 (4 lignes): C({N},4) = {len(list(combinations(range(N), 4)))} combinaisons")

    count_ordre3 = 0
    intersections_ordre3 = []

    for combo in combinations(lignes, 4):
        point, ok = trouver_point_convergence(combo, tolerance)

        if ok and point is not None:
            intensite_tot = sum(l.intensite for l in combo)
            intersections_ordre3.append((point, combo, intensite_tot))
            count_ordre3 += 1

    intersections_par_ordre[3] = intersections_ordre3
    print(f"    Trouvé: {count_ordre3} intersections")

    # Ordre 4+: 5+ lignes (limiter pour éviter explosion combinatoire)
    max_ordre = min(6, N)

    for ordre in range(4, max_ordre):
        n_lignes = ordre + 1
        n_combos = len(list(combinations(range(N), n_lignes)))

        if n_combos > 10000:
            print(f"\n  Ordre {ordre} ({n_lignes} lignes): Trop de combinaisons ({n_combos}), skip")
            continue

        print(f"\n  Ordre {ordre} ({n_lignes} lignes): C({N},{n_lignes}) = {n_combos} combinaisons")

        count = 0
        intersections_ordre = []

        for combo in combinations(lignes, n_lignes):
            point, ok = trouver_point_convergence(combo, tolerance)

            if ok and point is not None:
                intensite_tot = sum(l.intensite for l in combo)
                intersections_ordre.append((point, combo, intensite_tot))
                count += 1

        intersections_par_ordre[ordre] = intersections_ordre
        print(f"    Trouvé: {count} intersections")

    print()
    print("="*80)
    print()

    return intersections_par_ordre

# ============================================================================
# RENFORCEMENT NON-LINÉAIRE
# ============================================================================

def calculer_renforcement_non_lineaire(intersections_par_ordre, alpha=1.5):
    """
    Calcule le renforcement non-linéaire aux intersections

    Hypothèse: Φ_inter(n) ∝ (Σ I_lignes) · n^α

    Args:
        intersections_par_ordre: Dict des intersections
        alpha: Exposant non-linéaire (α>1 → renforcement)

    Returns:
        contribution_totale: Contribution au potentiel total
    """
    print("CALCUL DU RENFORCEMENT NON-LINÉAIRE")
    print("-"*80)
    print(f"Exposant α = {alpha}")
    print()

    contribution_totale = 0.0
    details = {}

    for ordre, intersections in intersections_par_ordre.items():
        n_lignes = ordre + 1  # Ordre 2 = 3 lignes, etc.

        contribution_ordre = 0.0

        for point, lignes, intensite_tot in intersections:
            # Renforcement: n^α
            facteur_renforcement = n_lignes**alpha

            # Contribution: (somme intensités) · n^α
            contrib = intensite_tot * facteur_renforcement

            contribution_ordre += contrib

        contribution_totale += contribution_ordre

        details[ordre] = {
            'n_intersections': len(intersections),
            'n_lignes': n_lignes,
            'contribution': contribution_ordre,
            'facteur_moyen': n_lignes**alpha
        }

        print(f"  Ordre {ordre} ({n_lignes} lignes):")
        print(f"    Intersections: {len(intersections)}")
        print(f"    Facteur renforcement: {n_lignes**alpha:.2f}×")
        print(f"    Contribution: {contribution_ordre:.2e}")
        print()

    print(f"Contribution totale: {contribution_totale:.2e}")
    print()

    return contribution_totale, details

# ============================================================================
# TEST COMPLET
# ============================================================================

def test_intersections_ordre2():
    """Test complet des intersections d'ordre 2+"""
    print("\n" + "="*80)
    print(" ANALYSE DES INTERSECTIONS D'ORDRE 2+ ".center(80))
    print("="*80)
    print()

    # Créer réseau
    print("Création du réseau 3D...")
    lignes = creer_reseau_3d(GALAXIES_3D, D_EFF)
    print(f"✓ {len(lignes)} lignes créées")
    print()

    # Détecter intersections
    tolerance_test = [5.0, 10.0, 20.0]

    resultats_tolerances = {}

    for tol in tolerance_test:
        print(f"\n{'='*80}")
        print(f" TEST AVEC TOLÉRANCE = {tol} kpc ".center(80))
        print(f"{'='*80}")

        intersections = detecter_intersections_ordre2_plus(lignes, tolerance=tol)

        # Calculer renforcement avec différents α
        alphas_test = [1.0, 1.5, 2.0]

        contributions = {}
        for alpha in alphas_test:
            print(f"\n--- α = {alpha} ---")
            contrib, details = calculer_renforcement_non_lineaire(intersections, alpha)
            contributions[alpha] = (contrib, details)

        resultats_tolerances[tol] = {
            'intersections': intersections,
            'contributions': contributions
        }

    # Résumé
    print("\n" + "="*80)
    print(" RÉSUMÉ ".center(80))
    print("="*80)
    print()

    for tol, res in resultats_tolerances.items():
        print(f"Tolérance = {tol} kpc:")
        intersections = res['intersections']

        total_intersections = sum(len(inter_list) for inter_list in intersections.values())
        print(f"  Total intersections: {total_intersections}")

        for ordre, inter_list in intersections.items():
            if len(inter_list) > 0:
                print(f"    Ordre {ordre}: {len(inter_list)} intersections")

        print()

    return resultats_tolerances, lignes

# ============================================================================
# VISUALISATIONS
# ============================================================================

def generer_visualisations_ordre2(resultats, lignes):
    """Visualisations des intersections d'ordre 2+"""
    print("="*80)
    print(" GÉNÉRATION DES VISUALISATIONS ".center(80))
    print("="*80)
    print()

    fig = plt.figure(figsize=(20, 12))

    # Utiliser tolérance = 10 kpc pour visualisation
    res = resultats[10.0]
    intersections = res['intersections']

    # ================================================================
    # Subplot 1: Distribution par ordre
    # ================================================================
    ax1 = plt.subplot(2, 3, 1)

    ordres = sorted(intersections.keys())
    counts = [len(intersections[o]) for o in ordres]

    ax1.bar(ordres, counts, color='blue', alpha=0.7, edgecolor='black')
    ax1.set_xlabel('Ordre d\'intersection', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Nombre d\'intersections', fontsize=12, fontweight='bold')
    ax1.set_title('Distribution des intersections par ordre\n(tolérance = 10 kpc)',
                 fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.set_xticks(ordres)
    ax1.set_xticklabels([f'{o}\n({o+1} lignes)' for o in ordres])

    # ================================================================
    # Subplot 2: Renforcement vs α
    # ================================================================
    ax2 = plt.subplot(2, 3, 2)

    alphas = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0]
    contributions_vs_alpha = []

    for alpha in alphas:
        contrib, _ = calculer_renforcement_non_lineaire(intersections, alpha)
        contributions_vs_alpha.append(contrib)

    ax2.plot(alphas, contributions_vs_alpha, 'g-', linewidth=2.5, marker='o', markersize=6)
    ax2.set_xlabel('Exposant α', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Contribution totale', fontsize=12, fontweight='bold')
    ax2.set_title('Effet de l\'exposant non-linéaire α',
                 fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.axvline(1.5, color='red', linestyle='--', alpha=0.5, label='α = 1.5')
    ax2.legend(fontsize=10)

    # ================================================================
    # Subplot 3: Contribution par ordre (α=1.5)
    # ================================================================
    ax3 = plt.subplot(2, 3, 3)

    contrib_alpha15, details = res['contributions'][1.5]

    ordres_avec_contrib = []
    contribs = []

    for ordre in sorted(details.keys()):
        if details[ordre]['n_intersections'] > 0:
            ordres_avec_contrib.append(ordre)
            contribs.append(details[ordre]['contribution'])

    ax3.bar(ordres_avec_contrib, contribs, color='purple', alpha=0.7, edgecolor='black')
    ax3.set_xlabel('Ordre d\'intersection', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Contribution au potentiel', fontsize=12, fontweight='bold')
    ax3.set_title('Contribution par ordre (α = 1.5)',
                 fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.set_xticks(ordres_avec_contrib)
    ax3.set_xticklabels([f'{o}\n({o+1} lignes)' for o in ordres_avec_contrib])

    # ================================================================
    # Subplot 4: Positions 3D des intersections ordre 2
    # ================================================================
    ax4 = plt.subplot(2, 3, 4, projection='3d')

    if 2 in intersections and len(intersections[2]) > 0:
        points_ordre2 = np.array([inter[0] for inter in intersections[2]])

        ax4.scatter(points_ordre2[:,0], points_ordre2[:,1], points_ordre2[:,2],
                   c='green', s=50, alpha=0.6, label='Ordre 2 (3 lignes)')

    # Galaxies
    for gal in GALAXIES_3D:
        pos = gal['position']
        ax4.scatter(pos[0], pos[1], pos[2], c='red', s=100,
                   alpha=0.8, edgecolors='black')

    ax4.set_xlabel('X (kpc)', fontsize=10)
    ax4.set_ylabel('Y (kpc)', fontsize=10)
    ax4.set_zlabel('Z (kpc)', fontsize=10)
    ax4.set_title('Positions 3D des intersections d\'ordre 2',
                 fontsize=12, fontweight='bold')
    ax4.legend(fontsize=9)

    # ================================================================
    # Subplot 5: Comparaison tolérances
    # ================================================================
    ax5 = plt.subplot(2, 3, 5)

    tolerances = sorted(resultats.keys())
    totaux_inter = []

    for tol in tolerances:
        inter = resultats[tol]['intersections']
        total = sum(len(inter_list) for inter_list in inter.values())
        totaux_inter.append(total)

    ax5.plot(tolerances, totaux_inter, 'b-', linewidth=2.5, marker='s', markersize=7)
    ax5.set_xlabel('Tolérance (kpc)', fontsize=12, fontweight='bold')
    ax5.set_ylabel('Nombre total d\'intersections', fontsize=12, fontweight='bold')
    ax5.set_title('Sensibilité à la tolérance',
                 fontsize=13, fontweight='bold')
    ax5.grid(True, alpha=0.3)

    # ================================================================
    # Subplot 6: Tableau résumé
    # ================================================================
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')

    # Résumé pour tol=10 kpc, α=1.5
    data = [
        ['Paramètre', 'Valeur'],
        ['', ''],
        ['Tolérance', '10.0 kpc'],
        ['Exposant α', '1.5'],
        ['', ''],
    ]

    for ordre in sorted(details.keys()):
        if details[ordre]['n_intersections'] > 0:
            data.append([f'Ordre {ordre} ({details[ordre]["n_lignes"]} lignes)',
                        f'{details[ordre]["n_intersections"]} inter.'])

    data.append(['', ''])
    data.append(['Total inter.', f'{sum(len(intersections[o]) for o in intersections)}'])
    data.append(['Contribution totale', f'{contrib_alpha15:.2e}'])

    table = ax6.table(cellText=data, cellLoc='left', loc='center',
                     colWidths=[0.6, 0.4])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    # Style
    table[(0, 0)].set_facecolor('#4CAF50')
    table[(0, 1)].set_facecolor('#4CAF50')
    table[(0, 0)].set_text_props(weight='bold', color='white')
    table[(0, 1)].set_text_props(weight='bold', color='white')

    ax6.set_title('Résumé Intersections Ordre 2+', fontsize=12, fontweight='bold', pad=20)

    # Sauvegarder
    plt.tight_layout()
    plt.savefig('simulation_intersections_ordre2.png', dpi=300, bbox_inches='tight')
    print("✓ Graphique sauvegardé: simulation_intersections_ordre2.png")
    print()

# ============================================================================
# PROGRAMME PRINCIPAL
# ============================================================================

if __name__ == "__main__":
    resultats, lignes = test_intersections_ordre2()
    generer_visualisations_ordre2(resultats, lignes)

    print("="*80)
    print(" FIN DE L'ANALYSE ".center(80))
    print("="*80)
    print()
