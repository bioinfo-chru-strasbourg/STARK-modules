#!/bin/csh
#
#		logo_web.csh <pid> <type>
#		<type> = ps | pdf | gif
#	Makes postscript logos & converts them to gif format
#  The files colors, marks, wave and default.amino.frq must be in
#  the current directory
#	

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

set pid = $1
set tmp = "../tmp"
setenv GS_LIB /opt/sfw/share/ghostscript
setenv LD_LIBRARY_PATH /usr/openwin/lib:/usr/dt/lib:/usr/lib:/usr/ucblib:/opt/sf
w/lib:/usr/local/lib
set path = ($path /opt/sfw/bin)

# Assume all logos have been made and post-processed 
# and the postcript files exist

#-------------------------------------------------------------------------
#	See if any logo files were created
#	If so, merge them all into one big postscript file
echo what huh
set logos = ($tmp/logo.$pid.*.ps)
echo wefa
if ($#logos < 1) then
	echo "Content-type: text/html"
	echo ""
	echo "Logo error: No logos produced"
	exit
endif

#-------------------------------------------------------------------------
#	Show the logos
if ($#argv < 2) then
   set type = "gif"
else
   set type = $2
endif

if ($type == "ps") then
   echo "Content-type: application/postscript"
   echo ""
   cat $tmp/logo.$pid*.ps
   # combine all postscript files into one
#   rm logo.$$* >& /dev/null
else
   if ($type == "pdf") then
      foreach logo ($logos)
#			      0-492 was the original logos
	convert -crop 648x300+0-492 ps:$logo pdf:$logo.pdf
#	rm -f $logo
      end
      echo "Content-type: application/pdf"
      echo ""
      set logos = ($tmp/logo.$pid*.pdf)
      montage +frame +shadow +label -geometry 648x300\!+0+0 -tile 1x$#logos $logos $tmp/$pid.logos.pdf
      cat $tmp/$pid.logos.pdf
      rm $tmp/$pid.logos.pdf $tmp/logo.$pid*.pdf >& /dev/null
   else if ($type == "gif") then
      #	Crop off blank bottom & convert to gif format
      #	Logos are about 663 pixels wide by 858 pixels high
      foreach logo ($logos)
	convert -crop 648x300+0-492 ps:$logo $logo.gif
#	rm -f $logo
        echo "Content-type: image/gif"
	echo ""
	cat $logo.gif
	exit (0)
  	end

      #-------------------------------------------------------------------------
      #	Group them on a page; it may be a long page
      set logos = (logo.$pid*.gif)
      montage +frame +shadow +label -geometry 648x300\!+0+0 -tile 1x$#logos $logos $tmp/$pid.logos.gif
      rm -f $logo.$pid*.gif

      if (-e $tmp/$pid.logos.gif) then
         #---------------Now show the combined logos ----------------------
         echo "Content-type: image/gif"
         echo "" 
         cat $tmp/$pid.logos.gif
      else
         echo "Content-type: text/html"
         echo ""
         echo "Logo error: No gif file produced"
      endif

#-------------------------------------------------------------------------

   endif #end if pdf 
endif # end if ps 
exit(0)
