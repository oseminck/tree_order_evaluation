# tree_order_evaluation
## About
Code presented in the communication *Evaluating Hierarchical Clustering Methods for Corpora with Chronological Order* presented at the conference without proceedings EADH 2021 by Olga Seminck, Philippe Gambette, Dominique Legallois and Thierry Poibeau.

The figure below illustrates what the script `tree_order_evaluation` does: it reorders the children of the internal nodes of the tree in order to minimise the number of conflicts between the order of the leaves of the tree and the lexicographic order of their labels (for example there is a conflict between leaves `1874a_Nouveaux_contes_a_Ninon` and `1865_La_confession_de_Claude`: the first one appears before the second one in the order of the leaves from top to bottom whereas the first one is ranked after the second one in the lexicographic - or chronologic - order). 

![Illustration of the input and output of tree_order_evaluation.py!](/figures/fig-ZolaEADH.jpg "Illustration of the input and output of tree_order_evaluation.py")

Note that the code currently ranks according to the lexicographic order of the leaf labels: we sometimes added letter `a` or `b` right after the year in the file name in order to minimise the number of conflicts in the same solution, when considering novels whose first publication year is identical. The trees in the figure above were drawn automatically from their Newick format using [Dendroscope](https://uni-tuebingen.de/fr/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/dendroscope/).

## Acknowledgement
This work was funded in part by the French government under the management of the Agence Nationale de la Recherche as part of the “Investissements d’avenir” program, references ANR-19-P3IA-0001 (PRAIRIE 3IA Institute) and ANR-16-IDEX-0003 (I-Site Future, programme “Cité des dames, créatrices dans la cité”).

## Data
Content of the `tree` folder:
* `inputBalzac.txt`, `inputZolaRougonMacquart.txt`: first test files generated from corpora of French novels by Honoré de Balzac and the Rougon Macquart cycle of 20 novels
* `inputCounterExample1.txt`, `inputCounterExample2.txt`: examples of distinct trees minimizing distinct criteria, 24 conflicts and an optimal number of 5 leaves to delete (e.g. 1, 3, 6, 7 and 8) for inputCounterExample1.txt versus an optimal number of 17 conflicts (1 conflicting with 2, 4, 5, 6 and 7; 3 with 4, 5, 6 and 7; 6 and 7 with 4 and 5; 10 with 9 and 11; 11 with 8 and 9) and 6 leaves to delete (e.g. 1, 3, 6, 7, 10, 11).
* `inputZola.txt`: test file generated from the corpus of French novels by Émile Zola extracted from [corpus CIDRE](https://www.ortolang.fr/market/corpora/cidre)
* `inputMoisl2020.txt`: tree provided in Figure 15 on page 16 of H. Moisl (2020): “[How to visualize high-dimensional data: a roadmap](https://doi.org/10.46298/jdmdh.5594)”, *Journal of Data Mining and Digital Humanities*, Special Issue on Visualisations in Historical Linguistics: 1–19.
* `inputSchoech2012.txt`: rooted version of the tree provided by C. Schöch at https://dragonfly.hypotheses.org/43

## Using the code
Please put the input tree in the Newick format in a file named input.txt in the same folder as the `script tree_order_evaluation.py`. 

The output file `output.txt` will then be generated. 

On the second line, it will contain the ordered version of this tree with the minimum number of conflicts, followed by the number of conflicts on the third line, followed by the minimum number of leaves to delete to get the remaining leaves in alphabetical order, followed by the set of leaves to keep.

The value of variable `testNb` may be changed in order to change the number of random orders generated to evaluate the probability of getting the observed number of conflicts (or less) by chance.

## Versions
* v1.0: computation of the minimum number of conflicts and significance evaluation by random simulations
* v2.0: added the computation of the minimum number of leaves to delete
* v2.1: corrected a bug (child permutation initialization) and saved the time taken to compute parameters
* v2.2: added a code to compute parameters for all files of a folder