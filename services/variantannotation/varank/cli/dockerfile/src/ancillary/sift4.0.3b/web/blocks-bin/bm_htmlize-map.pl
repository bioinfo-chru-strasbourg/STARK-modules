#!/usr/bin/perl
#
#	bm_htmlize-map.pl - /blocks-bin/ version
#		Can't use btest's version becuase wants "../tmp"
# Htmlize a block map drawing (from block_vis) and insert links to sequences
# (in the NCBI Entrez WWW database) and blocks (in file file.blocks).
# Script can work either as filter or just with arguments.
#
# Invocation: htmlize-map.pl file.vis file.maps
#             htmlize-map.pl file.maps < file.vis
#
# Written by Ross Morgan-Linial

# Get the block number
$blocknum = pop @ARGV;

# HTML header tags
print '<HTML>', "\n";
print '<TITLE>Block map</TITLE>', "\n"; 
print '<BODY>', "\n";
print '<PRE>', "\n";

# Various constant bits of code, to make them easy to change.
# The second command-line argument gets incorporated into $blockurl.
$tagstart = '<a href="';
$tagend = '">';
$blockurl = '/blocks-bin/showblock.sh?' . $blocknum . '+';
$sequenceurl = 'http://www3.ncbi.nlm.nih.gov:80/htbin-post/Entrez/query?uid=';
$options = '&form=6&db=p&Dopt=g';

# Loop over all the input lines
while (<>) {
    # Change '&', '<' and '>' for HTML 
    s/\&/\&amp;/g;
    s/\</\&lt;/g;
    s/\>/\&gt;/g;

    # Match a line containing actual data.
    # This is highly sensitive to the block_vis output format.
    if (/^(\s*)([\w\d]+)(\s+\(?[ \d]+\)?\s*)([^ \n]*)(.*)/) {
        # Print the first part (sequence name & length)
#       print $1, $tagstart, $sequenceurl, $2, $options, $tagend, $2, '</a>', $3;
        print "$1 $2 $3";

        # Save the map visualization in the magic variable $_
        $_ = $4;

        # Save the extra stuff
        $extra = $5;

        # Bracket the blocks with WWW links
        # We need to set $family, below, before doing this
#       s"(([A-Z])\2*)"$tagstart$blockurl$family$2$tagend$1</a>"g;

        # Reintroduce the newline removed in the regexp above
        print $_, $extra, "\n";
    }
    else { 
        # The block family accession need to be found before the above is done.
        # The accesion is identified as 7 non-whitespace chars at the 
        # begining of a line followed by ": " string.
        if (/^(\S{7}): /) {
            $family = $1;
        }

        # Print the line unchanged
        print;
    }
}

# Close HTML tags
print '</PRE>', "\n";
print '</BODY>', "\n";
print '</HTML>', "\n";

exit(0);