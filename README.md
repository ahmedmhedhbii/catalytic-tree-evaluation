# Tree Evaluation is in Space $\mathcal{O}(\log n \cdot \log \log n)$ Implementation

Voici la première version de l'implémentation de l'algorithme dans la section 4 (la "version d'échauffement"). Le code est en Python et sera expliqué (pas de manière détaillée pour l'instant mais plus tard).

## Fonctions et Explications

### La fonction `initialize_tree_and_catalyst()`
Initialise une instance de $TreeEval_{k, h}$, les registres d'une manière aléatoire avec des valeurs du $GF(2 ^{\lceil \log ( 2 \lceil \log k \rceil + 2) \rceil})$ comme proposé dans le papier et aussi $\omega$, qui est une racine primitive de notre corps qui va être toujours $GF(2)$. [Explication dans ce lien](https://math.stackexchange.com/questions/4104767/is-x-always-a-primitive-element-of-textgf2m).

### La fonction `initialize_field(k)`
Initialise le corps fini $GF(2 ^{\lceil \log ( 2 \lceil \log k \rceil + 2) \rceil})$

### La fonction `clean_computation(h)`
Appelle la fonction `recursive_clean` qui est l'équivalente de $P_u$ et prend une copie du registre de la racine, et fait le calcul pour donner le résultat final.

### La fonction `recursive_clean(node, level, h)`
C'est la fonction $P_u$ comme dans le papier, implémentée telle quelle.

### La fonction `q_u_i(y, z, f_u, field, i)`
C'est la fonction $q_{u,i}(y, z)$ comme dans le papier, implémentée telle quelle.

### La fonction `e_poly(y, beta, field)`
C'est la fonction $e(y, \beta) = \prod_{i=1}^{\lceil \log k \rceil} (1 - y_i + (2y_i - 1)\beta_i)$ comme dans le papier, implémentée telle quelle.

### La fonction `field_sum(bits, field)`
C'est pour faire le calcul pour le résultat final (la différence d'avant et après avoir fait le calcul).

## Remarques
- Les bits se lisent de gauche à droite par exemple $11_2 = 1011$ et non pas $1101$.
- `f'{m:0{n}b}'[::-1]` est une chaîne de caractères qui donne le nombre binaire `m` de taille `n` (avec des zéros à gauche).
- Le code n'est pas encore bien organisé (je vais faire ça plus tard), merci pour votre compréhension.
