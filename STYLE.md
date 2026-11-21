# STYLE.md — Boundary-controlled algebraic vorticity order in Shakti artificial spin ice

Rules applied to the whole manuscript (main.tex). New text must follow them.

## Terminology (use / do not use)
| Use | Do not use | Note |
|---|---|---|
| tetromino | L-tile, L-shape (as a noun) | "L-shaped tetromino" on first use |
| vorticity $v$ | circulation number, chirality | chirality = L vs J tetromino |
| unhappy vertex | frustrated vertex | |
| plaquette (horizontal / vertical) | cell (for Shakti) | "cell" = chessboard cell only |
| placement I / II | type I / II placement | |
| half vortex, full vortex (noun); half-vortex value (modifier) | half-vortex (noun) | hyphen only before a noun |
| arctic circle | Arctic circle | lower case throughout |
| frozen region / liquid region | ordered / disordered region | |
| block-aligned (edges along block boundaries) | | vs. edges through the blocks |
| basketweave tiling | brick-weave | flat states of the phase key |

## Abbreviations
| Abbreviation | Expansion | First defined |
|---|---|---|
| GFF | Gaussian free field | Sec. IV A |
| MC | (not used; write "Monte Carlo") | |
| ASI | (not used; write "artificial spin ice") | |

## Notation decisions
- US spelling (center, behavior, color).
- Cross references: Fig.~\ref{}, Figs.~..., Eq.~\eqref{} (with parentheses), Sec.~\ref{}, Table~\ref{}, Ref.~\cite{}; "Figure"/"Equation" only at sentence start.
- Ranges with en dash (`--`), e.g. $n=36$--$39$.
- Uncertainties as $a\pm b$.
- Approximation symbols: $\approx$ numerical value, $\simeq$ asymptotic/approximate equality, $\sim$ scaling, $\propto$ proportional.
- Subscript words in roman: $\Psi_{\rm sat}$, $T_{\rm on}$, $\lambda_{\rm d}$, $A_{\rm G}$.
- Display equations end with the sentence punctuation.

## Symbols
| Symbol | Meaning | Defined |
|---|---|---|
| $v$ | tetromino vorticity ($\pm1$) | Sec. II B |
| $h$, $\langle h\rangle$ | Thurston height, mean height | Sec. III A |
| $r(p)$ | residue of $h$ mod 4 at vertex $p$ (always written with its argument) | Eq. (1) |
| $s(d)$, $s_{\rm H}$, $s_{\rm V}$ | period-4 prefactor | Eq. (3) |
| $\varepsilon(p)$ | $(-1)^{\lfloor h/4\rfloor}$ | Sec. III A |
| $V$ | total vorticity | Sec. III C |
| $n$ | order of the Aztec diamond; $n=4m-1, 4m$ | Sec. II C, IV B |
| $L$; $W\times H$ | chessboard square side; rectangle (short side $H$) | Sec. II C |
| $H$ (holes) | hole size (only in Sec. V D and Fig. 13) | Sec. V D |
| $N$ | linear size in scaling laws | Sec. IV A, Eq. (10) |
| $x$, $x_v$ | scaling dimensions (half and full vortex) | Eq. (6), Sec. VII |
| $g$ | GFF variance coefficient | Sec. IV A |
| $A$ | local superstructure amplitude of $\langle v\rangle$ | Fig. 4 caption, Eq. (10) |
| $A_{\rm G}$ | Glaisher–Kinkelin constant | Eq. (9) |
| $C$ | asymptotic constant of the central amplitude | Sec. IV B |
| $\rho_\Omega$, $c_\Omega$, $\Omega$, $w$ | conformal radius, prefactor, region, position | Eq. (10) |
| $\Psi$, $\Psi_{\rm sat}$ | phase field, its frozen-region value | Eq. (11), Fig. 5 |
| $\ell$ | phase index of the flat states, $\arg\Psi=\ell\pi/4$ | Sec. V B |
| $k$ | winding sector of a hole | Sec. V D |
| $\bm q$, $\bm q_0$, $\bm\kappa$ | wave vector, peak position, offset | Eq. (17) |
| $J$ | unhappy-vertex energy | Sec. V E |
| $T_{\rm on}$, $\lambda$ | onset temperature, edge-plaquette fugacity | Eq. (15), App. D |
| $Q$ | boundary matrix of the Pfaffian | App. D |
| $\lambda_{\rm d}$ | defect fugacity | Sec. VII |
| $D$, $\mu$, $a$ | dipolar energy $\mu_0\mu^2/4\pi a^3$, island moment, lattice constant | Sec. VII |
| $T_c$ | dipolar locking temperature | Sec. VII |
| $z$ | vertex coordination number | Sec. II A |
