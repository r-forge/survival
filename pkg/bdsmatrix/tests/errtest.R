#
# Test out the error message in depth.s
#   add a couple of children to test1, that form a circl
#
kindepth(id = 1:15,
      mom.id=c( 0,15, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13, 13),
      dad.id=c( 0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10, 10))

# This should result in the "Impossible loop in pedigree" error
