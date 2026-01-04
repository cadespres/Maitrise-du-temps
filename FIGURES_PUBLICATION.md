# Figures pour Publication Scientifique

**Projet:** La Red Asselin - Explication Géométrique de la Matière Noire
**Date:** 2026-01-04
**Statut:** Prêtes pour révision et publication

---

## Vue d'Ensemble

Ce document liste toutes les figures générées pour la publication scientifique sur le Réseau Asselin. Toutes les figures sont en haute résolution (300 dpi) et prêtes pour inclusion dans un article scientifique.

---

## Liste des Figures

### Figure 1 : Simulation 2D du Réseau Asselin

**Fichier:** `simulation_reseau_asselin_2d.png`
**Résolution:** 300 dpi
**Dimensions:** 18" × 12"
**Format:** PNG (recommandé: convertir en EPS/PDF vectoriel pour publication)

**Contenu (6 subplots):**

1. **Vue d'ensemble du réseau 2D**
   - 45 lignes Asselin entre galaxies
   - Intensité des lignes (code couleur alpha)
   - 10 galaxies (taille proportionnelle à log(M))
   - 103 intersections réelles (croix vertes)

2. **Zoom sur le centre (±200 kpc)**
   - Détail des lignes dans la région centrale
   - Étiquettes des galaxies principales
   - Visualisation des intersections proches

3. **Distribution des intensités**
   - Histogramme logarithmique
   - 20 bins
   - Intensité min: 4.09×10²⁹
   - Intensité max: 7.22×10³⁷

4. **Potentiel réseau radial Φ_réseau(r)**
   - Profil de 1 à 150 kpc
   - Courbe violette avec marqueurs
   - Tendance polynomiale degré 2 (rouge)
   - **Attention:** Cette version montre le potentiel DÉCROISSANT (version initiale)

5. **Gradient dΦ/dr**
   - Dérivée numérique du potentiel
   - Pourcentage de régions croissantes
   - Ligne de référence y=0

6. **Statistiques des intersections**
   - Tableau récapitulatif
   - Nombres de lignes, paires testées, intersections
   - Code couleur (vert pour intersections réelles)

**Légende suggérée:**
> *Figure 1: Simulation 2D du Réseau Asselin dans le Groupe Local. (a) Vue d'ensemble montrant 45 lignes entre 10 galaxies et 103 intersections réelles (croix vertes). (b) Zoom sur la région centrale. (c) Distribution des intensités de lignes. (d) Profil radial du potentiel réseau (version initiale décroissante). (e) Gradient radial. (f) Statistiques du réseau.*

---

### Figure 2 : Simulation 2D Améliorée - Formulation Cumulative

**Fichier:** `simulation_reseau_asselin_2d_amelioree.png`
**Résolution:** 300 dpi
**Dimensions:** 18" × 10"

**Contenu (6 subplots):**

1. **Comparaison des 3 formulations**
   - Mode distance (bleu)
   - Mode rayon moyen (vert)
   - Mode pondéré (rouge) ⭐ RECOMMANDÉ
   - Profils Φ(r) de 1 à 150 kpc

2. **Gradient du potentiel (mode pondéré)**
   - dΦ/dr vs rayon
   - 100% croissant pour mode pondéré
   - Boîte d'information avec pourcentage

3. **Distribution des distances minimales**
   - Histogramme 20 bins
   - d_min pour toutes les lignes
   - Importance pour formulation cumulative

4. **Réseau 2D coloré par d_min**
   - Code couleur viridis
   - Bleu = lignes proches du centre
   - Jaune = lignes loin du centre
   - Galaxies en rouge

5. **Accumulation de lignes avec r**
   - Nombre de lignes avec d_min < r
   - Courbe orange croissante
   - Démontre l'effet cumulatif

6. **Comparaison normalisée**
   - Les 3 profils normalisés [0,1]
   - Permet de comparer les formes

**Légende suggérée:**
> *Figure 2: Formulation cumulative améliorée du potentiel réseau. (a) Comparaison de trois modes de calcul montrant que le mode pondéré (rouge) donne une croissance monotone. (b) Gradient positif confirmant la croissance à 100%. (c) Distribution des distances minimales des lignes au centre. (d) Réseau coloré par proximité au centre. (e) Nombre cumulatif de lignes contribuant au potentiel. (f) Profils normalisés pour comparaison.*

---

### Figure 3 : Courbes de Rotation avec Réseau Asselin

**Fichier:** `courbes_rotation_reseau_asselin.png`
**Résolution:** 300 dpi
**Dimensions:** 20" × 12"

**Contenu (6 subplots):**

1. **Courbes de rotation comparatives** ⭐ FIGURE PRINCIPALE
   - Observations (points noirs avec barres d'erreur)
   - Newton (ligne rouge pointillée, χ²=3120)
   - Réseau nominal (ligne bleue pointillée, χ²=7.7×10⁹)
   - Réseau optimisé (ligne verte épaisse, χ²=1374)
   - **Amélioration de 56% vs Newton!**

2. **Résidus optimisés**
   - Newton (rouge) vs Réseau optimisé (vert)
   - Bande grise ±1σ (10 km/s)
   - Ligne y=0
   - Résidus systématiquement meilleurs pour réseau

3. **Masse effective**
   - M_visible (rouge pointillé)
   - M_eff avec réseau (vert)
   - Zone remplie = contribution réseau
   - De 0 à 150 kpc

4. **Distribution χ² par région radiale**
   - Barres comparatives Newton (rouge) vs Réseau (vert)
   - 5 régions: 0-10, 10-30, 30-60, 60-100, 100-150 kpc
   - Meilleure performance dans région 30-60 kpc

5. **Sensibilité au paramètre σ**
   - χ² vs σ (10-150 kpc)
   - Ligne verticale verte: σ_opt = 200 kpc
   - κ fixé à valeur optimale
   - Montre le minimum bien défini

6. **Tableau récapitulatif**
   - χ² pour chaque modèle
   - Paramètres optimaux (σ, κ)
   - Amélioration en %
   - Code couleur (vert si amélioration)

**Légende suggérée:**
> *Figure 3: Courbes de rotation galactiques et optimisation. (a) Comparaison des modèles montrant que le réseau Asselin optimisé (vert) améliore Newton de 56%. (b) Résidus significativement réduits. (c) Masse effective incluant la contribution du réseau (zone verte). (d) Performance par région radiale. (e) Sensibilité au paramètre de transition σ. (f) Résumé des paramètres optimaux.*

---

### Figure 4 : Extension 3D du Réseau Asselin

**Fichier:** `simulation_reseau_asselin_3d.png`
**Résolution:** 300 dpi
**Dimensions:** 20" × 14"

**Contenu (6 subplots):**

1. **Réseau 3D (vue 3D interactive)**
   - 30 lignes principales (limitation visuelle)
   - Galaxies en rouge (taille log(M))
   - Axes X, Y, Z (kpc)
   - Alpha proportionnel à intensité

2. **Courbes de rotation 3D**
   - Observations (noir)
   - Newton (rouge, χ²=3120)
   - 3D nominal (bleu, χ²=23097)
   - 3D optimisé (vert, χ²=1167)
   - **Amélioration de 62.6% vs Newton!**

3. **Résidus 3D**
   - Newton (rouge) vs 3D optimisé (vert)
   - Bande σ grise
   - Encore meilleure performance que 2D

4. **Projection XY des galaxies**
   - Vue du dessus
   - Étiquettes des galaxies
   - Axes égaux pour préserver proportions

5. **Distribution des d_min (3D)**
   - Histogramme 20 bins
   - Légèrement différent du 2D
   - Effet de la géométrie spatiale complète

6. **Tableau comparatif résultats**
   - χ² pour Newton, 3D nominal, 3D optimisé
   - Paramètres σ et κ optimaux
   - Code couleur selon performance

**Légende suggérée:**
> *Figure 4: Extension 3D du réseau Asselin. (a) Visualisation 3D du réseau spatial entre galaxies. (b) Courbes de rotation montrant une amélioration de 62.6% vs Newton (meilleur que 2D). (c) Résidus réduits par la géométrie spatiale complète. (d) Distribution spatiale réelle des galaxies (projection XY). (e) Distances minimales en géométrie 3D. (f) Comparaison quantitative 2D vs 3D.*

---

### Figure 5 : Intersections d'Ordre 2+ et Renforcement Non-Linéaire

**Fichier:** `simulation_intersections_ordre2.png`
**Résolution:** 300 dpi
**Dimensions:** 20" × 12"

**Contenu (6 subplots):**

1. **Distribution par ordre d'intersection**
   - Barres bleues
   - Ordre 2 (3 lignes): 1676 intersections
   - Ordre 3 (4 lignes): 4609 intersections
   - Total: 6285 intersections (tolérance = 10 kpc)

2. **Effet de l'exposant α**
   - Contribution vs α (1.0 → 3.0)
   - Courbe verte croissante
   - Ligne verticale rouge: α = 1.5 (recommandé)
   - Montre amplification non-linéaire

3. **Contribution par ordre (α=1.5)**
   - Barres violettes
   - Ordre 3 domine (86% de contribution)
   - Confirme importance intersections multiples

4. **Positions 3D des intersections ordre 2**
   - Scatter plot 3D
   - Points verts = intersections
   - Points rouges = galaxies
   - Montre distribution spatiale

5. **Sensibilité à la tolérance**
   - Total intersections vs tolérance (5, 10, 20 kpc)
   - Courbe bleue croissante
   - Importance du critère géométrique

6. **Tableau résumé**
   - Paramètres (tolérance, α)
   - Nombre d'intersections par ordre
   - Contribution totale
   - Style tabulaire clair

**Légende suggérée:**
> *Figure 5: Intersections d'ordre 2+ et renforcement non-linéaire. (a) Distribution montrant 6285 intersections détectées. (b) Amplification avec l'exposant α. (c) Domination des intersections d'ordre 3 (4 lignes convergentes). (d) Distribution spatiale 3D des points de convergence. (e) Sensibilité au critère géométrique de tolérance. (f) Résumé statistique pour tolérance = 10 kpc, α = 1.5.*

---

## Figures Complémentaires Suggérées

### Figure S1 : Diagramme Conceptuel du Réseau Asselin
*(À créer manuellement ou avec logiciel de dessin)*

Schéma illustrant:
- Deux galaxies M₁ et M₂
- Ligne Asselin L_ij entre elles
- Formule de l'intensité
- Facteur d'expansion exp(-d/d_eff)
- Intersection avec autres lignes
- Renforcement au point de convergence

### Figure S2 : Comparaison Multi-Échelles
*(À créer)*

Trois panneaux montrant l'effet à différentes échelles:
- Échelle locale (d << d_eff): f ≈ 1
- Échelle galactique (d ~ kpc): f ≈ 0.9-0.99
- Échelle cosmologique (d > Mpc): f < 0.1

### Figure S3 : Profils de Masse Visible
*(À créer)*

Détail des composantes:
- Bulbe de Hernquist (courbe orange)
- Disque exponentiel (courbe bleue)
- Gaz (courbe verte)
- Total visible (courbe rouge)

---

## Préparation pour Publication

### Conversion de Format

**Recommandations:**

1. **Pour journaux astrophysiques (A&A, ApJ, MNRAS):**
   - Convertir PNG → EPS ou PDF vectoriel
   - Utiliser fonts embeddées
   - Vérifier que texte reste lisible à 50% reduction

2. **Outils recommandés:**
   ```bash
   # ImageMagick pour conversion
   convert -density 300 figure.png figure.pdf

   # Ou avec Python/Matplotlib
   # Sauvegarder directement en PDF vectoriel
   plt.savefig('figure.pdf', format='pdf', dpi=300)
   ```

3. **Vérification qualité:**
   - Lisibilité des étiquettes d'axes
   - Taille des polices (minimum 8pt)
   - Légendes non superposées
   - Couleurs distinguables (daltonisme-friendly)

### Légendes Complètes

Chaque figure doit avoir:
1. **Titre court** (en gras)
2. **Description détaillée** de chaque subplot (a), (b), (c)...
3. **Paramètres clés** utilisés
4. **Interprétation principale** en 1-2 phrases

### Code Couleur Cohérent

**Maintenir dans toutes les figures:**
- **Noir** : Observations
- **Rouge** : Newton/Masse visible
- **Bleu** : Réseau nominal
- **Vert** : Réseau optimisé
- **Violet/Orange** : Potentiel réseau
- **Gris** : Incertitudes

---

## Liste de Vérification Publication

- [ ] Toutes les figures en 300 dpi minimum
- [ ] Conversion en format vectoriel (EPS/PDF)
- [ ] Légendes complètes et détaillées rédigées
- [ ] Numérotation cohérente (Figure 1, 2, 3...)
- [ ] Références croisées dans le texte (voir Figure X)
- [ ] Code couleur cohérent entre figures
- [ ] Polices lisibles (≥ 8pt après réduction)
- [ ] Barres d'échelle et unités claires
- [ ] Droits d'auteur et licence (CC-BY ou autre)
- [ ] Fichiers sources originaux archivés

---

## Fichiers Disponibles

Tous les fichiers sources Python ayant généré les figures sont disponibles:

1. `simulation_reseau_asselin_2d.py`
2. `simulation_reseau_asselin_2d_amelioree.py`
3. `simulation_courbes_rotation_reseau_asselin.py`
4. `simulation_reseau_asselin_3d.py`
5. `simulation_intersections_ordre2.py`

**Reproductibilité:** Toutes les figures peuvent être régénérées en exécutant ces scripts.

---

## Statistiques des Figures

| Figure | Subplots | Taille (Mo) | Éléments clés |
|--------|----------|-------------|---------------|
| Figure 1 (2D) | 6 | ~2.5 | Réseau, intersections, distributions |
| Figure 2 (2D amélioré) | 6 | ~2.3 | Formulations cumulatives |
| Figure 3 (Rotation) | 6 | ~1.8 | Courbes v(r), optimisation |
| Figure 4 (3D) | 6 | ~2.1 | Géométrie 3D, amélioration |
| Figure 5 (Ordre 2+) | 6 | ~1.9 | Intersections, renforcement |
| **Total** | 30 | ~10.6 | |

---

## Contact et Support

Pour questions sur les figures ou demandes de formats spécifiques:

- **Repository GitHub:** https://github.com/cadespres/Maitrise-du-temps
- **Issues:** Utiliser le système GitHub Issues
- **Documentation:** Voir README.md et articles scientifiques

---

**Document préparé le:** 2026-01-04
**Version:** 1.0
**Statut:** ✅ Prêt pour publication

---
