# App Shiny sur GitHub Pages

App Shiny (R) exportée en version autonome navigateur avec [shinylive](https://posit-dev.github.io/r-shinylive/), déployée automatiquement sur GitHub Pages via GitHub Actions.

## Comment ça marche

GitHub Pages n'héberge que des fichiers statiques : l'app est donc exécutée côté navigateur via [WebR](https://docs.r-wasm.org/webr/latest/) (R compilé en WebAssembly). Aucun serveur R n'est nécessaire, mais le premier chargement télécharge l'environnement R (~20-40 Mo), ce qui prend quelques secondes.

## Structure

```
.
├── app/
│   ├── app.R            # Le code de l'app Shiny (c'est ici qu'on modifie)
│   ├── dataSetAll.csv   # Genotypes (separateur ';', colonne 'phenotype' = USH|RP)
│   ├── domains.csv      # Domaines proteiques USH2A (UniProt O75445, isoforme 1, 5202 aa)
│   └── favicon.png      # Icône de l'onglet
└── .github/
    └── workflows/
        └── deploy.yml   # Build + déploiement automatique vers GitHub Pages
```

## Domaines protéiques

Les bandeaux colorés le long des axes du contact-map proviennent de `app/domains.csv`, généré depuis l'annotation [UniProt O75445 (usherin)](https://www.uniprot.org/uniprotkb/O75445) : les répétitions Fibronectin type-III consécutives ont été agrégées en blocs. Pour ajuster les frontières, éditer directement `name,start,end` en residues 1-5202.

## Paramètres d'URL

Le formulaire peut être pré-rempli via l'URL : `https://char-al.github.io/shinytest/?pheno=USH&r1=759&r2=934`

## Mise en place du dépôt

1. Créer un dépôt GitHub (public ou privé — Pages fonctionne pour les deux si GitHub Actions est la source de déploiement).
2. Pousser ce contenu sur la branche `main`.
3. Dans le dépôt : **Settings > Pages > Source** = **GitHub Actions**.
4. Le workflow `Deploy Shiny app to GitHub Pages` se déclenche à chaque push sur `main`.
5. L'app est disponible sur `https://<owner>.github.io/<repo>/`.

## Développement en local

```r
install.packages(c("shiny", "shinylive"))

# Prévisualiser l'app comme elle apparaîtra en ligne
shinylive::export("app", "site")
serve_site <- function() shiny::runApp("site")
```

L'app peut aussi être prévisualisée classiquement (plus rapide) :

```r
shiny::runApp("app")
```

## Notes

- Packages R supplémentaires : les ajouter dans le step `Install R packages` du workflow. Seuls les packages installés localement dans WebR seront disponibles dans l'app — vérifier la compatibilité des packages non-High vers WebR.
- Les libs shinylive sont Exportées dans `site/` (non versionné) et régénérées à chaque build — ne pas les commiter, `.gitignore` s'en charge.
