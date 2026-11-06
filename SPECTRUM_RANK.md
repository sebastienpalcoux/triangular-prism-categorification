# The minimum rank for the spectrum criteria

The minimum is four for either criterion, and the same ring attains it for both. The lower bound uses only the six nonzero boundary coefficients and the two integers measuring the sizes of the two pentagon sums.

## The lower bound

Ostrik’s ring-theoretic parametrization gives one non-self-dual rank-three case, the group ring of the cyclic group of order three, categorified by graded vector spaces. Rank one and rank two are self-dual. In a self-dual fusion ring multiplication is commutative, and Frobenius reciprocity makes $t(a,b,c)=N_{ab}^{c}$ symmetric. [Ostrik, Proposition 3.1](https://doi.org/10.17323/1609-4514-2015-15-2-373-396)

Write $(a,b,c,u,v,w,p,q,r)=(i_1,\ldots,i_9)$. Admissibility means that $t(a,u,w),t(b,u,v),t(c,v,w),t(a,p,r),t(b,p,q),t(c,q,r)$ are all positive. Put

$$
S=\sum_h N_{up}^h N_{vq}^h N_{wr}^h,\qquad O=N_{ab}^c.
$$

A zero spectrum obstruction requires $(S,O)=(0,1)$, and a one spectrum obstruction requires $(S,O)=(1,0)$.

**Invertible labels.** Simultaneous permutation of the three pairs $(u,p),(v,q),(w,r)$, with the corresponding permutation of the labels of the joining sides, and interchange of the two ends preserve admissibility, $S$ and $O$. Thus an invertible label $g$ can be placed at $u$ or at $a$. Self-duality implies $g^2=1$.

- If $u=g$, then $w=ga$, $b=gv$. Hence $O=N_{vw}^c>0$, and the $h=gp$ term in $S$ is $N_{bp}^qN_{pr}^a>0$.
- If $a=g$, then $w=gu$, $r=gp$ and $O=\delta_{b,gc}$. Both $b$ and $gc$ occur in $uv$ and in $pq$. With the simple basis orthonormal, associativity and reciprocity give
  $$
  S=\sum_h (N_{up}^h)^2 N_{vq}^h
  \ge \langle up,vq\rangle=\langle uv,pq\rangle
  \ge 2-\delta_{b,gc}=2-O.
  $$

Both cases preclude the required pair. If there is just one noninvertible $R$, all nine labels must therefore be $R$. Then $O=N_{RR}^R>0$ by admissibility, and $1$ contributes to $S$.

**Two noninvertibles.** The remaining case has basis $1,X,Y$ and table

$$
X^2=1+mX+kY,\qquad XY=kX+\ell Y,\qquad
Y^2=1+\ell X+nY,\qquad k^2+\ell^2=1+kn+\ell m.
$$

If $k=0$, then $\ell(\ell-m)=1$, so $\ell=1$, $m=0$ and $X$ is invertible. The analogous argument gives $\ell>0$, and the displayed identity implies $m+n>0$. Thus, if $m>0$, all three nonunit products contain $X$, and if $n>0$, they all contain $Y$. We have $S>0$.

If $O=0$, symmetry of $t$ and positivity of $k,\ell$ force $a=b=c=X$ with $m=0$, after exchanging $X,Y$ if necessary. Now $n>0$. Each of the triples $(u,v,w)$ and $(p,q,r)$ contains at most one $X$, since every pair in each triple must occur at a vertex with label $X$. If a spectral product is $X^2$, the other two are $Y^2$, and the three products share $1$ and $Y$. Otherwise every spectral product is $XY$ or $Y^2$, and they share $X$ and $Y$. Thus $S\ge2$, as required.

## A rank-four witness for each criterion

In Proposition 6.10 the four basis elements are numbered from 1 to 4. In the verifier they are numbered from 0 to 3. In either order, write them as $(1,x,y,z)$. Their products are

| Product | Value | Product | Value |
|---|---|---|---|
| $x^2$ | $1+x+z$ | $y^2$ | $1+x+y$ |
| $xy$ | $y+z$ | $yz$ | $x+2z$ |
| $xz$ | $x+y+2z$ | $z^2$ | $1+2x+2y+3z$ |

The table is on line 33531 of the decompressed rank-four census, with the historical identifier recorded in the classification crosswalk. [Census](https://github.com/sebastienpalcoux/fusion-ring-census/blob/afe3ffa0a371daf5131f733c9459c2eb85e18c83/results/rank4/tables.txt.gz) · [Crosswalk](https://github.com/sebastienpalcoux/grothendieck-ring-classification/blob/ebaeb62e80eaeb9082dd92d5645c218745923809/supplements/rank4/data/through32_mapping.json.gz)

For the zero spectrum criterion take $(i_1,\ldots,i_9)=(y,y,x,x,y,z,y,y,y)$. The spectral products $xy,y^2,zy$ have empty common support, $N_{yy}^x=1$, and the two lists of alternative sums are $(1,2,3)$ and $(3,1,1)$. For the one spectrum criterion take $(i_0,\ldots,i_9)=(y,z,y,y,x,y,y,y,y,x)$. The spectral products $xy,y^2,yx$ have the single common constituent $y$, each with multiplicity one, $N_{yz}^{y}=0$, and the three lists of alternative sums are $(3,1,1)$, $(1,1,3)$ and $(3,2,1)$. All six boundary coefficients equal one for each witness.

## Independent finite check

Run `python3 verify_spectrum_rank.py`. The verifier evaluates the displayed rank-four certificates, checks the fusion-ring axioms and counts twelve witnesses of each kind. For the lower bound it truncates every coefficient to $0,1,2$, with $2$ denoting any value at least two. All tests use only positivity or equality to $0$ or $1$ of sums of products. Hence truncated arithmetic preserves every test and every necessary associativity equality. Reciprocity makes the nonunit coefficient tensor symmetric, leaving $0,1,4$ free entries in ranks one through three. Exactly $1,3,50$ truncated tables, respectively, satisfy truncated associativity. The exhaustive test on all of them finds no witness, without any bound on the original multiplicities.
