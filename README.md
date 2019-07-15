# Nearest neighbor matching on several popular distance measures

The repo contain the Stata program `matchME` for matching observations based on several commonly used distance measures.

## Installation

To load `matchME`, include the following line in your do file:

```
qui do "https://raw.githubusercontent.com/goshevs/matchME/master/matchME.ado"

```

### Syntax

```
	matchME varlist, dist(string) [cutoff(real) smp(integer) force]
```
<br>

**Required arguments**


| argument    | description            |
|-------------|------------------------|
| *varlist*   | the set of variables to match on |
| *dist*      | distance measure. Currently supposted distances are `Euclieadean', `nEuclidean` (normalized Euclidean) and `Mahalanobis` |

<br>

**Optional arguments**


| argument       | description            |
|----------------|------------------------|
| *cutoff*       | maximum distance within which to search for matching observations; default is no cutoff |
| *smp*          | size of matching pool (i.e. the number of candidates to match from); default is 3|
| *force*        | match beyond reciprocal matches |

## Examples

Please, see example.do in this repo