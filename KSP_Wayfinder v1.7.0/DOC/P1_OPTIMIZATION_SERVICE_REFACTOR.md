# P1 — extraction prudente d'OptimizationService

## Objectif

Réduire progressivement la classe `Wayfinder` sans modifier son API publique,
les formats SQLite, ni le comportement numérique de l'optimiseur.

La stratégie retenue est la composition avec une façade de compatibilité : les
utilisateurs continuent d'appeler `Wayfinder`, tandis que les responsabilités
techniques migrent vers de petits services testables.

## Tranche 1

Le module `_OptimizationService.py` reçoit uniquement six helpers sans état :

- calcul conservateur du nombre de workers ;
- limitation automatique du nombre d'îles ;
- construction des topologies Pygmo ;
- suivi de l'origine des identifiants de population ;
- télémétrie de migration et de topologie ;
- équilibrage du budget SADE / recuit simulé.

Les anciens noms restent disponibles sur `Wayfinder` et délèguent au service.
Les appels internes et externes existants restent donc compatibles.

## Explicitement hors périmètre

- aucune modification des algorithmes Pygmo ;
- aucune modification de la construction de l'entonnoir ;
- aucune modification du fitness, de l'éjection ou du décodage ;
- aucune modification du schéma SQLite ou des jobs ;
- aucune modification des signatures publiques.

## Tranche 2

La politique de construction du plan d'entonnoir a ensuite été déplacée dans
`OptimizationService`. La façade conserve `_funnel_stage_plan()` avec la même
signature et délègue au service ; l'ancien corps a été supprimé.

Avant le déplacement, une empreinte SHA-256 du JSON normalisé a été capturée
pour chacune des neuf variantes supportées : `legacy`, `local`, `hybrid`, les
deux variantes phase-elites et les quatre noms scout/archive. Les tests
vérifient désormais ces neuf références, y compris les budgets et paramètres
imbriqués.

## Vérification

- tests de comportement existants conservés ;
- test explicite de délégation de la façade ajouté ;
- suite complète après la tranche 2 : 54 tests et 9 sous-tests réussis ;
- compilation de `WayfinderCore` et `Tests` réussie.

## Tranche 3

Les helpers de diversité et d'archive ont été déplacés micro-étape par
micro-étape :

1. sélection de seeds exacts, entièrement pure ;
2. embedding de phase, avec `bodies_by_name` fourni explicitement ;
3. sélection d'élites et mise à jour d'archive, avec leurs dépendances injectées
   sous forme de callbacks.

Cette injection évite toute référence du service vers `Wayfinder` et conserve
la possibilité de remplacer l'embedding dans un test ou une spécialisation. Les
wrappers historiques restent en place sur la façade.

Une incohérence découverte pendant l'extraction a été corrigée séparément :
`select_exact_diverse_seeds(..., count=0)` renvoie désormais `[]`, tandis que
`count=1` conserve le meilleur candidat. Les valeurs négatives et une liste de
candidats vide sont également couvertes.

La suite atteint maintenant 57 tests et 9 sous-tests. `_Wayfinder.py` est passé
d'environ 3 734 à 3 485 lignes au cours des trois tranches.

## Tranche 4

Les métadonnées nécessaires à la reconstruction d'un run historique sont
désormais persistées au niveau SQL `runs` en schéma v15 :

- `effective_optimizer_seed`, y compris quand le job n'avait pas de seed
  explicite ;
- `funnel_config_json`, instantané canonique de la configuration complète du
  funnel ;
- `funnel_config_hash`, pour grouper et comparer les runs sans décoder tout le
  JSON ;
- `code_revision`, rempli best-effort depuis Git ;
- `planet_pack_hash`, calculé depuis les constantes orbitales chargées.

`OptimizationService` expose aussi la correspondance canonique entre stratégie
publique et stratégie exacte de funnel, ainsi qu'un constructeur de config de
run. Un test dédié vérifie qu'un run historique peut être relancé avec la même
seed effective et le même hash de configuration.

À l'issue des tranches actuelles, la suite atteint 104 tests et 36 sous-tests.

## Tranche 5

`StageConfig` et `FunnelConfig` introduisent une représentation typée mais
transitoire du funnel. Ces objets ne changent pas encore l'exécution Pygmo :
ils servent à produire exactement les dictionnaires canoniques déjà consommés
par `Wayfinder` et persistés en JSON dans SQLite.

La règle de compatibilité reste stricte : `StageConfig.to_dict()` et
`FunnelConfig.to_dict()` doivent reproduire bit-à-bit les structures historiques
validées par les hashes de tests. Le moteur continue donc de lire des dicts,
tandis que la couche configuration dispose désormais d'un point d'ancrage plus
propre pour la future GUI.

La validation de `StageConfig` couvre maintenant les invariants minimaux utiles
à une édition GUI sûre : tailles positives, population compatible SADE,
algorithmes et modèles d'éjection connus, topologie reconnue si elle est
définie, taux de migration non négatif, et options structurées pour
`adaptive_stop` / `annealing`. Cette validation reste compatible avec les dicts
historiques : les hashes de plans existants restent inchangés.

## État après les expériences L0

Le funnel configurable supporte désormais la cascade de pression et une option
expérimentale 16+4 à anneaux séparés. Les invariants de cette topologie sont
validés par `StageConfig`; sa partition, son pont et sa télémétrie sont
persistés pour permettre le replay. Les sorties volumineuses de benchmark sont
ignorées par Git, tandis que les graphes retenus vivent dans `DOC/assets`.

La couche indépendante en amont de L0 est maintenant posée : scout
tree/Tisserand puis filtre Lambert réel sur le premier arc. Elle ne touche ni au
moteur Pygmo ni à la persistance du funnel. La tranche suivante doit mesurer le
recall et la réduction de jobs sur plusieurs fenêtres avant de transformer les
points Lambert survivants en boîtes T0/TOF et en jobs SQLite.
