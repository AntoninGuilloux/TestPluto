## Des idées en vrac pour la visualisation


### copier coller de la discussion slack

Ce qu'on cherche à construire, c'est une partie finie de l'arbre donnée par l'automate. A chaque feuille, on a donc une chaine finie de caractères p, un cycle c, et un triangle, comme tu l'as codé.
Il est facile de retrouver à partir du triangle:
- l'intervalle L_p: c'est un des intervalles définit par le triangle.
- la transformation associée à la chaine (formule ci dessous). Appliquer cette transformation inverse revient à changer la racine de l'arbre, ce qui peut être pratique si notre parcours de l'arbre est trop déséquilibré, comme ça semble être le cas en pratique.

Ca semble donc être la bonne structure de donnée: un arbre, dont les arêtes sont étiquetées par le caractère associé à la symmétrie, et des feuilles décorées par (p,c,L_p,triangle_p).
Pour les noeuds internes, je ne sais pas s'il faut stocker toute cette info (p,L_p,triangle_p) ou bien seulement partiellement, à voir.

Remarque que l'arbre est "planaire: pour chaque noeud, il y a un ordre cyclique sur la famille (le parent, les enfants) suivant l'endroit ou on pointe dans le cercle  ->il faut implémenter cette fonction: étant donné un des intervalles de feuille, rendre le suivant/le précédent. C'est un algorithme simple de parcours d'arbre planaire.

Comme ça, quand on a un intervalle J, on descend dans dans l'arbre jusqu'à soit trouver une feuille telle que L est inclus dans J: c'est gagné ; soit trouver la feuille L qui contient le milieu de J. On applique alors la transformation inverse associée, pour se ramener à la situation où on part de l'origine.

**La formule**: Si on a le triangle initial T0 $(-1,1,e^{-i \pi/3},etc)$ et que le nouveau triangle est $(a_-,a_+,b_-,etc)$, associé à un mot pair. L'homographie est celle qui envoie $(-1,1,e^{-i pi/3})$ sur $(a_-,a_+,b_-)$.
Pour la formule, sauf erreur de ma part, on note
- $u = \frac{e^{i\pi/3}-1}{e^{i\pi/3}+1}$,
- $v = \frac{b_- - a_+}{b_- - a_-}$
- $w(t) = \frac{t+1}{t-1}$ (pour t un complexe)
Alors l'homographie envoie $t$ sur
$$ s = \frac{a_+ u w(t) - a_-u}{w(t)u-v}$$
Cerise sur le gateau: il semble que w(t) explose quand t tend vers 1, mais si on change le choix des trois points parmi les 6 qui décrivent les triangles, on fait une variante qui permet d'éviter que le dénominateur n'explose.




##### Anecdotique, probablement inutile

- si x, y, z sont des nombres complexes de module 1, alors le produit P(x,y,z) = (x/y-1)(y/z-1)(z/x-1) est imaginaire pur et la partie imaginaire est positive ssi ils sont bien ordonnées. autrement dit, between(x,y,z) = sign(imag(P(x,y,z)))