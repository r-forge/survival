#
#  a script to help the comparison of testall.v4 ( version 3.4)
#   and testall.out (Splus verion 5) by deleteing unwanted stuff
#
tr \' \" <testall.save | sed -e 's/observations//' > temp1
sed -e 's/>[ >]*//' -e 's/^+/ /' -e 's/observations//' < testall.out | tr \' \" > temp2
diff -b  temp1 temp2 > zed

# At this point "zed" contains the diffs.  Almost all of them are
#   changes in how S echos the input lines.

# Now get rid of even more stuff: programming lines, blanks, and comments
sed -e'/<-/d' -e'/[{}+~\!=]/d' -e'/Call/d' -e'/^$/d'< temp1 > temp3
sed -e'/<-/d' -e'/[{}+~\!=]/d' -e'/Call/d' -e'/^$/d'< temp2 > temp4
sed -e's/^[ >]*//' -e'/^#/d'< temp3 > temp1
sed -e's/^[ >]*//' -e'/^#/d'< temp4 > temp2

diff -w temp1 temp2 > zed2
