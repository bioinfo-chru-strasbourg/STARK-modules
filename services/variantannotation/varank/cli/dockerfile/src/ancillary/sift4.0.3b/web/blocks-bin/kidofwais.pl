#!/usr/bin/perl

#	kidofwais.pl is executed by getblock.sh to look up blocks
#	in the index. It runs getblock to add a title.
# 3/15/06 Use getblock.pl instead of getblock.sh

$cgibin = "blocks-bin";       # added by BJA

# turn this on (=1) if you want the block descriptions to be shown.  This
# slows down the retrieval a LOT.  The reason is that for each block it calls
# getblock to get the ID and AC lines.   
# If you don't want the descriptions turn this off (=0);
#   -- BJA

$show_block_description = 1;


# kidofwais.pl -- WAIS search interface
#   (derived from sonofwais.pl, derived from wais.pl)
#
# a url of "http://your.www.server/cgibindir/kidofwais.pl" will search the
# default WAIS index; a url of the form:
# "http://your.www.server/cgibindir/kidofwais.pl/another_wais_src"
# will search the WAIS index called "another_wais_src".

# NOTE: To really make the best use of this script, create a "Source table".
# You can read about a "Source table" at
# "http://ewshp2.cso.uiuc.edu/multiindex_post.txt", and
# you can see a sample "Source table" at
# "http://ewshp2.cso.uiuc.edu/Source_table" .

# do we need to debug the returned lines from waisq? If so, write into log.
$DEBUG = 0;     # set to 1 to turn on debugging code; set next to logfile
$debugLOG = "/var/info/www/httpd/logs/cgi_log";


# where is your waiq binary?
$waisq ="./waisq";

# where are your source files?
$waisd = "../index";

# what database do you want to search? (Default)
$default_src = "blocks";
$blocksdb = "../data-blocks/ID_DE.dat";
$printsdb = "../data-prints/ID_DE.dat";

# what is the opening title you want to present to users (default)
$openingTitle = "BLOCKS Database Keyword Search";

# after searching, what do you want the title to be? (default)
$closingTitle = $openingTitle;

# Use a table to look up title (and other info) to use on Search page?
# If 1, then should supply filename for table.
$use_Source_table = 0;
$Source_table = "/var/info/www/wais-sources/Source_table";

# Specify the directory where your WWW docs reside
# (This is the same path you subtracted when you waisindexed (using -t URL)
# if you are indexing your whole site. If not, then this gives enough info
# to find the top of the docs tree, and the "Source table" will give enough
# info to find a particular set of docs and its index. Or you can make these
# point to a particular directory within your docs tree, and a particluar index,
# but then this copy of "kidofwais" can only be used to search that particular
# index. So do take a look at the "Source table" concept I created to go along
# with this script -- see "http://ewshp2.cso.uiuc.edu/Source_table" .)
$wwwDocpath = "../www/";

# Specify the url for this WWW server
# (The same comments as just above for $wwwDocpath apply.)
$serverURL = "http://blocks.fhcrc.org/";

# maximum number of hits to return
$max_hits = 60;

# Should we use the "highlighting script" for filetypes that "make sense" for
# it? (1 is yes, 0 is no)
$use_hilite = 0;   # not using since the "document" is from a CGI script.

# specify the www url to the "highlighting script" (used for .html and .txt)
# whether this is used is controlled by $use_hilite flag.
$hilite_script = "http://blocks.fhcrc.org/$cgibin/print_hit_bold.pl/";

# specify the "first hit anchor" to be used with hiliting
$anchor = "#first_hit";

# who maintains this service?
$maintainer = "<I>webmaster\@blocks.fhcrc.org</I>";

# and when was it last modified.
$modified = "Nov 1997";

# you shouldn't have to edit anything below this line, except if you want to change the help text

sub extractTitle {
	# try and get the <title> ... </title> field from file
	# only try to find it in the first 5 lines, and then give up
	local($fl) = @_;
	local($intitle) = 0;
	local($title) = "$theFile: No title, please let $maintainer know.";
	local($linenum) = 1;
	local($_);

	# read the file and extract the title
	# this is a combination of code from Eric Lease Morgan and Sean Dowd
	open (FP, "$fl") || return "File $theFile can't be read. Please contact
$maintainer.";

	while (<FP>) {
		chop;
		last if ($linenum > 5);
		$linenum ++;
		if (/<TITLE\s?>(.*)<\/TITLE\s?>/i) { # all on one line
			$title = $1;
			last;
		}
		elsif (/<TITLE\s?>(.*)$/i) {	# on multiple lines
			$title = $1;
			$intitle = 1;
		}
		elsif (/^(.*)<\/TITLE\s?>/i) {	# finish of multiple lines
			$title = "$title$1";
			$intitle = 0;
			last;
		}
		elsif ($intitle) {		# add to title, and keep going
			$title = "$title$_";
		}
	}
	close (FP);
	$title =~ s/^\s*//; # remove whitespace at front
	$title =~ s/\s*$//; # remove whitespace at end
	$title = "$theFile: Empty title, please let $maintainer know." unless $title;
	return $title;
} # end sub extractTitle

sub extractTableTitle { # Get file titles from table. Add a call to this
	# routine for any filetype that you decide to create titles for.
	# Currently, only the PDF type is set up to call this. (See the
	# subroutine &type_file). If nothing is found, it returns the current
	# $doc_title unchanged (which is the filename relative to doc root).
	# Table read into array at first reference to it.
	local ($fl) = @_;
	return $doc_title if ($file_title_table eq "");
	local($_,$name_to_find,$table_entry,$filename,$filetitle);
	# Change the next line to "$name_to_find = $theFile" if you want to
	# use the whole path relative to doc root as the name to lookup in the
	# table. The current code uses just the filename w/o path info.
	$name_to_find = substr($fl, rindex($fl, '/') + 1); # "basename"
	if ($current_file_title_array ne $file_title_table) { # Read in new
							      # table to array
		undef %title_array;	# erase current array
		$current_file_title_array = $file_title_table;
		open (TABLE_TITLE, $file_title_table) || return $doc_title;
		while ($table_entry = <TABLE_TITLE>) {
			chop;
			next if $table_entry =~ /^\s*#/;  # skip comments
			next if $table_entry =~ /^\s*$/;  # skip blank lines
			next if $table_entry !~ /~/;  # skip non-tilde lines
			($filename, $filetitle) = split(/~/, $table_entry, 2);
			$title_array{$filename} = $filetitle;
		}
		close (TABLE_TITLE);
	}
	if ($filetitle = $title_array{$name_to_find}) { return $filetitle;}
	else { return $doc_title;}
}

sub send_index {
    print "Content-type: text/html\n\n";
    
    print "<HEAD>\n<TITLE>$openingTitle</TITLE></HEAD>\n";
    print "<BODY>\n<H2><IMG SRC=\"/blocks/icons/small-blocks.xbm\">", $openingTitle, "</H2>\n";
    print "<ISINDEX prompt=\"Enter keywords and press RETURN:   \"><P>\n";

    print "[<A HREF=\"/blocks/getblock.html\">Get Block by number</A>] ";
    print "[<A HREF=\"/blocks\">Blocks Home</A>] ";
    print "[<A HREF=\"/blocks/blocks_search.html\">Block Searcher</A>] ";
    print "[<A HREF=\"/$cgibin/LAMA_search.sh\">LAMA Searcher</A>] ";
    print "[<A HREF=\"/blocks/make_blocks.html\">Block Maker</A>] ";
    print "<P><HR>\n";

    print "<H3>Help</H3>";
    print "To search, simply enter a query -- type the words you want to ";
    print "search for in the '<em>search term box</em>' (or follow the ";
    print "instructions from your client as to how to enter a search). If you enter two or ";
    print "more words, than any document that contains ANY of the words will ";
    print "be found. (For example, if you search for <B>ice cream</B>, then ";
    print "any document with either <B>ice</B> or <B>cream</B> will be on the ";
    print "list.) If you want only documents containing BOTH <B>ice</B> AND ";
    print "<B>cream</B>, then use <B>AND</B> between the words; your query ";
    print "would be: <B>ice and cream</B>. See the list below for a full ";
    print "description of searching options and their effect.<p> ";
    print "The result of this search will be a list of the documents which ";
    print "satisfy your search request (or a reply that indicates that NO ";
    print "documents match). The list is ordered by 'relevancy rank' (highest ";
    print "first), which very roughly corresponds to the number of times the ";
    print "word you searched for occurred in the document versus the size of ";
    print "the document.<p>\n";
    print "The following is a list of all the search options available:\n";
    print "<DL>\n";
    print "<DT><B>Right-hand truncation</B> (stemming) queries -- Expand Your Search\n";
    print "<DD>The query 'astro*' will find documents containing the words";
    print " 'astronomy' as well as 'astrophysics'. If you are not sure of the\n";
    print " spelling or form of the word used in the documents, use this.<p>\n";
    print "<DT>'<B>And</B>' queries -- Narrow Your Search\n";
    print "<DD>The query 'red and blue' will find the <B>intersection</B> of all";
    print " the documents containing the words 'red', and 'blue'.";
    print "The use of 'and' limits your retrieval.<p>\n";
    print "<DT>'<B>Or</B>' queries -- Expand Your Search\n";
    print "<DD>The query 'red or blue' will find the <B>union</B> of all the";
    print " documents containing the words 'red' and 'blue'.";
    print "The use of 'or' increases your retrieval.<p>\n";
    print "<DT>'<B>Not</B>' queries -- Narrow Your Search\n";
    print "<DD>The query 'red not green' will find all the documents containing";
    print " the word 'red', and <B>excluding</B> the documents containing the word 'green'.";
    print "The use of 'not' limits your retrieval.<p>\n";
    print "<DT><B>Nested</B> queries -- Complicated Searches\n";
    print "<DD>The query '(red and green) or blue not pink' will find the union of all";
    print " the documents containing the words 'red', and 'green'. It will then add (union)";
    print " all documents containing the word 'blue'. Finally, it will exclude all documents";
    print " containing the word 'pink'.\n";
    print "</DL>\n";
    print "<HR>";
    print "This page is maintained by $maintainer, and it was last modified on $modified.<p>\n";
    print "</BODY>\n";
}

sub type_file {
	# Set file type based on file extension; also has the "side effect" of
	# modifying the $url_to_use if $use_hilite flag is turned on and the
	# file type is appropriate to use it.
	# You can add other types if you want very easily. Also, if you have
	# created a "file_title_table" and want it to be accessed for a 
	# particular filetype, just add "$doc_title = &extractTableTitle ($_);"
	# to the appropriate do{} structure below. See the "pdf" file type
	# below for an example.
	local($filename) = @_;
	local($type) = "";
	local($_);
	SUFFIX: for ($filename) {
		/\.html$/i	&& do {
				 $type = "HTML file";
				 $doc_title = &extractTitle ($the_full_File);
				 if ($use_hilite) {
				  $url_to_use = $hilite_script . $filename .
						"?" . $query_plus . $anchor;
				 }
				 last SUFFIX;
				};
		/\.te?xt$/i	&& do {
				 $type = "text file";
				 if ($use_hilite) {
				  $url_to_use = $hilite_script . $filename .
						"?" . $query_plus . $anchor;
				 }
				 last SUFFIX;
				};
		/\.gif$/i	&& do {
				 $type = "GIF graphic";
				 last SUFFIX;
				};
		/\.ps$/i	&& do {
				 $type = "PostScript file";
				 last SUFFIX;
				};
		/\.pdf$/i	&& do {
				 $type = "PDF";
				 # try to get a better file title
				 $doc_title = &extractTableTitle ($_);
				 last SUFFIX;
				};
		/\.jpg$/i	&& do {
				 $type = "JPEG graphic";
				 last SUFFIX;
				};
		/\.mpg$/i	&& do {
				 $type = "MPEG movie";
				 last SUFFIX;
				};
		/\.Z$/i		&& do {
				 $type = "compressed file";
				 last SUFFIX;
				};
		/\.gz$/i	&& do {
				 $type = "compressed file";
				 last SUFFIX;
				};
		/\.au$/i	&& do {
				 $type = "Sun audio file";
				 last SUFFIX;
				};
		/\.hqx$/i	&& do {
				 $type = "Binhex file";
				 last SUFFIX;
				};
		/\.tar$/i	&& do {
				 $type = "tar'red file";
				 last SUFFIX;
				};
#		$type = "Unknown type"; # "fall thru" default case
		$type = "Block"; # "fall thru" default case
		do {
		    #$type = "HTML file";
		    #$doc_title = "BLOCK: $filename";
		    if ($use_hilite) {
			$url_to_use = $hilite_script . $filename .
			    "?" . $query_plus . $anchor;
		    }
		    last SUFFIX;
		};
	} # end suffix
	return $type;
} # end sub type_file

sub byscores { $scores[$b] <=> $scores[$a];} # descending numeric sort routine

sub print_it {	# Print out the hit list, or if multiple sources, save all
		# the info to be printed in a array (where all the info for
		# "one" hit becomes one array element), only to be printed
		# once all sources are searched. This is so we can sort the
		# combined hit lists into descending order by score.
	local($line) = @_;
	if ($line eq "</UL>\n") {	# All sources have been searched, we
					# can now sort and print the entire
					# hit list array.
		print @hit_list[ sort byscores $[..$#hit_list ];
		print "$line";
	} elsif ($line eq "</LI>\n") {
		# End of one "hit" listing, save in the hit array.
		$output_hit .= $line;
		push ( @hit_list, $output_hit );
		push ( @scores, $score );
		$output_hit = "";
	} else {	# more stuff for the same hit
		$output_hit .= $line;
	}
}

sub do_wais {

    $src = $default_src;
    # if 'PATH_INFO' has a non-null value, then use it as the name of the
    # WAIS source to search, otherwise will default to $default_src.
    $path_extension = $ENV{'PATH_INFO'};
    if ($path_extension =~ /^\/(.+)$/) { $src = $1; }
    $kidofwaisREF = $ENV{'SCRIPT_NAME'} . "/$src";
    @Sources = ( $src );	# Initialize array of sources to search
    if ( $use_Source_table ) {	# Read in the Source table info into
				# associative array.
	open (INDEX_TITLES, $Source_table) || last;
	while (<INDEX_TITLES>) {
		chop;
		next if /^\s*#/;	# skip comments
		next if /^\s*$/;	# skip blank lines
		($src_name, $remainder) = split(/~/, $_, 2);
		$src_array{$src_name} = $remainder;
	}
	close (INDEX_TITLES);
	($src_title, $src_multiple, $src_prefix, $file_title_table, $sources,
	   $go_to_url, $go_to_title) = split(/~/, $src_array{$src});
	if ($src_title ne "") {
		$openingTitle = "Search of $src_title";
		$closingTitle = "Search results from $src_title";
	}
	$src_multiple && (@Sources = split(/,/, $sources)); # Store the sources
    }							    # to be searched.

    do { &send_index; return; } unless defined @ARGV; # No search terms yet.

    local(@query) = @ARGV;
    local($pquery) = join(" ", @query);
    # NCSA's HTTPD puts backslashes in front of "funny" or "dangerous"
    # characters in the input supplied thru argv. In the case of search terms
    # for WAIS, this can screw up the search (parens and "*" get backslashed
    # and then don't work correctly). So remove the backslashes, AND the
    # potentially "dangerous" characters ( ; ` ! ).
    $pquery =~ tr/!\;\`\\//d;		# just in case, get rid of ;`! and \
    @query = split(' ',$pquery);	# and recreate query word array
    $query_plus = join("+", @query);

    print "Content-type: text/html\n\n"; # Start the "html" doc to be returned

    print "<HEAD>\n<TITLE>$closingTitle</TITLE></HEAD>\n";
    print "<BODY>\n<H2><IMG SRC=\"/blocks/icons/small-blocks.xbm\">", $closingTitle, "</H2>\n";
    print "\n<ISINDEX><P>\n";

#@@@@@
    print "<A HREF=\"/blocks/getblock.html\">[Getblock]</A>   ";
    print "<A HREF=\"/\">[Return to BLOCKS Home Page]</A><P>";

    print "Note that you can enter a new query in the <em>search term box</em> \n";
    print "from this screen/page without having to go back.\n";
    if ($use_Source_table && $go_to_url) {
	print "Another option is to \n";
	print "<A HREF=\"$go_to_url\">go to the $go_to_title</A>.\n";
    }

    local($hits, $score, $headline, $lines, $bytes, $type, @types, $date);
    $DEBUG && do { open (LOG, ">>$debugLOG") || die "can't open log";};

    foreach $src (@Sources) {	# Search each indicated index for the terms
      ($src_prefix, $file_title_table) = (split(/~/, $src_array{$src}))[2,3];
      open(WAISQ, "-|") || exec ($waisq, "-c", $waisd, "-s", $waisd, "-m",
		$max_hits, "-f", "-", "-S", "$src.src", "-g", @query);
      while (<WAISQ>) {
	$DEBUG && print LOG $_;
        /:score\s+(\d+)/ && ($score = $1);
        /:number-of-lines\s+(\d+)/ && ($lines = $1);
        /:number-of-bytes\s+(\d+)/ && ($bytes = $1);
        /:type "(.*)"/ && (push (@types, $1));
        /:headline "(.*)"/ && ($headline = $1);#%%%%
        /:date "(\d+)"/ && ($date = $1, $hits++, &docdone);
      }
      close(WAISQ);
      $total_hits += $hits;
      $hits = 0;
    }
    #&print_it ("</UL>\n");	# signal to print out hit array if we've been
				# building it (for multiple sources).
    if ($total_hits == 0) {
        print "<P>\n";
	print "<B>No items found that match your search query.</B> You can enter \n";
	print "another query in the <em>search term box</em> if you want \n";
	print "to try searching for something else. If you would like to \n";
	print "see a description of the search options again, you can \n";
	print "<A HREF=\"$kidofwaisREF\">go back to the main search page</A> \n";
	print "for this index.\n";
	print "<HR>\n";
    } elsif ($total_hits >= $max_hits) {
	print "<p>The following are the first $total_hits items that match \n";
	print "your query <B>\`$pquery\'</B>. Note that there may be more \n";
	print "items that match that are not shown (the search is limited \n";
	print "to $max_hits matches). You might want to further qualify \n";
	print "your search (use AND and NOT) to limit the matches. \n";
	print "<UL>\n";
	&print_it ("</UL>\n");
    } else {
	print "<p>The following $total_hits item(s) match your query \n";
	print "<B>\`$pquery\'</B>:\n";
	print "<UL>\n";
	&print_it ("</UL>\n");
    }

    print "</BODY>\n";	# End the "html" doc being returned
}

sub docdone {	# Called for each "hit" returned by waisq
    local($endfile,$path_to_file,$file_proper,$file_ext,$multi_type,$alt_count);
    if ($headline =~ /Search produced no result/) {
	if ($src_multiple) {  # don't print source listing if
	   $hits--;           # multi-index search
	}
	else {
           print "</UL><P>\n";
	   print "<B>No items found that match your search query.</B> You can enter \n";
	   print "another query in the <em>search term box</em> if you want \n";
	   print "to try searching for something else. If you would like to \n";
	   print "see a description of the search options again, you can \n";
	   print "<A HREF=\"$kidofwaisREF\">go back to the main search page</A> \n";
	   print "for this index.\n";

	   # Hack to exit after the above is printed (since we are using only
           # one DB this should be OK) -- BJA
           exit(0);  # added by BJA
	}
    } elsif (($headline =~ /^Information on database:/) ||
			($headline =~ /^Catalog for database:/)) {
	$hits--;
    } else {	# this is a "real" hit
	($endfile, $path_to_file) = split(' ', $headline);
        if ($path_to_file ne "") { # Not indexed with -t url, so headline of
				   # form:  "filename /path/to/file/"
		# Multitype indexed files will probably have this form, as at
		# least I can't get "-t url" to coexist with "-M type,type"
		$the_full_File = $theFile = $path_to_file.$endfile;
#		$theFile =~ s/^.*$wwwDocpath//i; # changed by BJA
		$theFile = $endfile;
#		$url_to_use = $serverURL.$theFile;  # changed by BJA
		$url_to_use = $serverURL.$cgibin.'/getblock.pl?'.$endfile;
	} else { # should have been indexed as "-t url", so headline of form
		 # http://your_server_url/path/to/actual/file
	
		# get the string to munge
		$theFile = $url_to_use = $headline;

		# parse out the file name (remove the server URL from the front)
		$theFile =~ s/^.*$serverURL//i;
	
		# concatenate the "wwwDocpath" variable with the file name
		$the_full_File = $wwwDocpath.$theFile;
	}

	$last_period = rindex($theFile, ".");	# need filename without .ext if
	if ($last_period > 0) {			# it turns out to be multitype
		$file_proper = substr($theFile, 0, $last_period);
		$file_ext = substr($theFile, $last_period + 1);
	}
	$doc_title = $theFile;
	$type  = &type_file ($theFile); # also modifies $url_to_use if flag
					# $use_hilite is set and right filetype
	$src_multiple && ($doc_title = "$src_prefix $doc_title");
	if ($bytes < 1000) { $calc_bytes = "&lt 1 Kbyte"; }
	else { $calc_bytes = int(($bytes + 500)/1000) . " Kbytes"; }

	# this is what is printed on each line -- BJA
        &print_it ("<LI><DL><DT><A HREF=\"$url_to_use\">$doc_title</A>");
        # changed the output to remove the type specification -- BJA
	&print_it (&block_info($doc_title));
	&print_it ("</DL>");
	&print_it ("\n");

	if (($#types > 0) && ($file_proper ne "")) {  # Multitype indexing
	  $alt_count = 0;			      # offer alternatives
	  foreach $multi_type (@types) {
		next if $multi_type eq $file_ext; # skip, already listed
		&print_it ("<BR>...Alternate Types Available: ") if $alt_count == 0;
		$alt_count ++;
		$theFile = "$file_proper.$multi_type";
		$url_to_use = $serverURL.$theFile;
		$type  = &type_file ($theFile);
		&print_it ("<A HREF=\"$url_to_use\">$type</A>, ");
	  } # end foreach $multi_type
	}
	&print_it ("</LI>\n");
    }
    $score = $headline = $lines = $bytes = $type = $date = '';
    $file_proper = $file_ext = '';
    @types = ();
}

sub block_info {
    local($block) = @_;
    local($tmp, $grp_rec, $grp_ac, $grp_id, $grp_de);

    $tmp = $grp_ac = $grp_id = $grp_de = "";
    if ($show_block_description) {

        if ($block =~ m/^PR/) {
           open(GREP, "grep \"^$block\" $printsdb |");
        }
        else  {
           open(GREP, "grep \"^$block\" $blocksdb |");
        }
        while($grp_rec = <GREP>)
        {
           ($grp_ac, $grp_id, $grp_de) = split(/\s+/, $grp_rec);
        }
        close(GREP);
        if ($grp_ac eq $block)
	{ $tmp = " <DD> <STRONG>".$grp_id."</STRONG> ".$grp_de; }
    }
    $tmp;
}

open (STDERR,"> /dev/null");
eval '&do_wais';

