.BG D
.FN survexp.us
.TL
Census data sets for the expected survival and person years functions
.PP
.AG us
total United States population, by age and sex, 1960 to 1980.
.AG uswhite
United States white population, by age and sex, 1950 to 1980.
.AG usr
United States population, by age, sex and race, 1960 to 1980.  Race is
white, nonwhite, or black.  For 1960 and 1970 the black population values
were not reported separately, so the nonwhite values were used.
.AG mn
total Minnesota population, by age and sex, 1970 and 1980.
.AG mnwhite
Minnesota white population, by age and sex, 1960 to 1980.
.AG fl
total Florida population, by age and sex, 1970 and 1980.
.AG flr
Florida population, by age, sex and race, 1970-1980.  Race is
white, nonwhite, or black.  For 1970 the black population values
were not reported separately, so the nonwhite values were used.
.AG az
total Arizona population, by age and sex, 1970 and 1980.
.AG azr
Arizona population, by age, sex and race, 1970-1980.  Race is
white versus nonwhite.  For 1970 the nonwhite population values
were not reported separately.
In order to make the rate table be a matrix, the 1980 values were
repeated.  (White and non-white values are quite different).
.PP
Each of these tables contains the daily hazard rate for a matched subject
from the population, defined as -log(1-q)/365.24 where q is the 1 year
probability of death as reported in the original tables.
For age 25 in 1970, for instance, p = 1-q is is the probability that a subject
who becomes 25 years of age in 1970 will achieve his/her 26th birthday.
The tables are recast in terms of hazard per day entirely for computational
convenience.
(The fraction .24 in
the denominator is based on 24 leap years per century.)
.PP
Each table is stored as an array, with additional attributes, and
can be subset and manipulated as standard S arrays.
Interpolation between calander years is done "on the fly" by the survexp
routine.  Judging from past experience, the 1990 data should become
available in 1995 or 96.
.PP
Some of the deficiencies, e.g. 1970 Arizona non-white, are a result of
local conditions.  The data probably exists, but we don't have a copy
it in the library.
.PP In November 1994 all of the tables were augmented to contain extrapolated
values for 1990 and 2000.  The details can be found in technical report 55.
.KW survival
.WR


