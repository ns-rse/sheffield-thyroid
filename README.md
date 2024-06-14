# Sheffield Thyroid Cancer Study

This repository contains the [Quarto](https://quarto.org) manuscript for work on the Sheffield Thyroid Cancer Study
which is the PhD work of [Ovie Edafe (mdp21oe)](https://github.com/mdp21oe).


## GitHub Pages

Because the raw data can not be included in this repository the manuscript can not be rendered via a GitHub Action.

In order to publish the pages you therefore need to do so locally.

1. Complete work on branches.
2. Check the pages render correctly using the _Render_ button in RStudio (shown when you are editing the `index.qmd`
   file).
3. Git commit and push changes to GitHub.
4. Create a Pull Request and merge the changes into the `main` branch.
5. Switch to the `main` branch locally and pull the merged changes down.
6. in RStudio in the terminal run the following...

``` bash
quarto publish
```

7. Follow the prompts and the site should be rendered and pushed to the `gh-pages` which is setup to be the basis of the
   GitHub hosted pages which can be viewed at
   [ns-rse.github.io/sheffield-thyroid](https://ns-rse.github.io/sheffield-thyroid/). This link can be shared with
   others to see the state of progress.
