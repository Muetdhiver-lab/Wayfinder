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

La suite atteint maintenant 58 tests et 9 sous-tests.

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

## Tranche suivante proposée

Ne pas déplacer immédiatement `_run_sqlite_funnel_job`, qui mélange encore
orchestration, persistance et construction Pygmo. La prochaine étape devrait
d'abord identifier une petite frontière d'exécution testable — probablement la
construction d'une population ou d'un algorithme d'île — sans ouvrir ni écrire
de base SQLite.

Critère d'arrêt : tout besoin de modifier simultanément une politique, une
signature ou un résultat numérique reporte l'extraction à une branche dédiée.
