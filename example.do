********************************************************************************
*** Example

*** Variables to match on
local myMatchVars "ageM sexM canyoureadM canyouwriteM ageSD sexSD canyoureadSD canyouwriteSD"


*** Euclidean distance
matchME `myMatchVars', dist(Euclidean) smp(5) cutoff(0.5) force

*** Normalized Euclidean distance
matchME `myMatchVars', dist(nEuclidean) smp(20) force

*** Mahalanobis distance
matchME `myMatchVars', dist(Mahalanobis) smp(7)

