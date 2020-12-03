#!/bin/csh
#

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
cd ../tmp
#ls -lgt > x  --- edit x leaving those files to be removed
awk '{print $8}' x > in

set tot = `wc -l in | awk '{print $1}'`
if ($tot > 960) then
   @ n = 1
   @ rem = $tot
   while ($rem > 960)
           @ h = $n * 960
           head -$h in | tail -960 > in.$n
           @ n += 1
           @ rem -= 960
   end
   tail -$rem in > in.$n
else
   cp in in.1
endif

set lists = (in.[1-9]*)
foreach list ($lists)
   set cgs = (`cut -f1 $list`)
   foreach cg ($cgs)
        rm $cg
   end
end
unalias rm
rm x in in.[1-9]*
df -sk
exit
