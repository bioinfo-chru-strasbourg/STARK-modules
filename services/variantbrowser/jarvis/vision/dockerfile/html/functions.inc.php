<?PHP

# Release 0.9.2
# Date 01/04/2015

# RELEASE NOTE
# 0.9.2.1-01/04/2015: ListSAMPLE updated. List of Sample for a run directly from the database

function QC($sample_id,$metrics=array(),$FASTQC_metrics_name=array("FASTQC")) {

	require "connect.php";

	$QC=array();
	$metrics_query="";
	if ($metrics!=array()) {
		$metrics_query=" AND QC.metrics IN ('".join("','",$metrics)."') ";
	};#if
	#print $metrics_query."<BR>";

	$query="SELECT QC.id as QC_id, QC.ref AS QC_ref, QC.type, QC.metrics AS QC_metrics, vcf.aligner AS vcf_aligner, vcf.pipeline, QC.QC AS QC_QC, sample.ref AS sample_ref, manifest.ref AS manifest_ref
		FROM vcf
		INNER JOIN QC_association ON (QC_association.object_id=vcf.id AND object_type='vcf')
		INNER JOIN QC ON (QC.id=QC_association.QC_id)
		INNER JOIN sample ON (sample.id=vcf.sample_id)
		INNER JOIN manifest ON (manifest.id=sample.manifest_id)
		WHERE vcf.sample_id IN ($sample_id)
		   $metrics_query
		GROUP BY vcf.sample_id, vcf.aligner, QC.metrics #QC.id
		"; # , QC.QC #AND project.ref='$project'
	#print "$query<BR>";

	$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
   	$summary= "<table class=' '>\n";

	while ($row = mysql_fetch_assoc($result)) {
		#echo "<pre>"; print_r($row); print "</pre><BR>";
		$QC_id=$row['QC_id'];
		$QC[$QC_id]=$row;

		$QC_metrics=$row['QC_metrics'];
		#print "$QC_metrics<BR>";
		$QC_QC=$row['QC_QC'];
		$vcf_aligner=$row['vcf_aligner'];
		#print $vcf_aligner;

		$QCFile_array=array();
		$QCFile_array_tmp=explode("\n",$QC_QC);
		#echo "<pre>"; print_r($QCFile_array_tmp); print "</pre><BR>";
		foreach ($QCFile_array_tmp as $k=>$l) { $QC[$QC_id]["QCarray"][$k]=explode("\t",trim($l)); };

		#if ($QC_metrics==$FASTQC_metrics_name) {
		#print_r($FASTQC_metrics_name); print_r($QC_metrics); echo "<BR>";
		if (in_array($QC_metrics,$FASTQC_metrics_name) || $QC_metrics==$FASTQC_metrics_name) {
			#print("TEST");
			$QCFASTQC_id=$QC_id;
			for ($i=4; $i<=9; $i++) {
				$QCData_array[$QC[$QC_id]["QCarray"][$i][0]]=$QC[$QC_id]["QCarray"][$i][1];
			};#for
			$module_name="";
			$module_for_warnings=array("Basic Statistics","Per base sequence quality","Per sequence quality scores","Per base N content","Sequence Length Distribution");

			foreach ($QC[$QC_id]["QCarray"] as $k=>$l) {
				# Warning data
				if (in_array(strtolower($l[1]),array("fail","warn")) && (in_array(str_replace(">","",$l[0]),$module_for_warnings) || 0)) { # pass ok #Basic Statistics
					$warnings["Raw Data"][]=strtoupper($l[1])." ".str_replace(">","",$l[0]);
				};#if
				# Module Splitting
				$l0_splitted=preg_split("//",$l[0]);
				# Start or End module
				if ($l0_splitted[1]==">" && $l0_splitted[2]==">") {
					if ($l[0]==">>END_MODULE") {
					} else {
						$module_name=str_replace(">","",$l[0]);
						$QC[$QC_id]["QCmodule"][$module_name]="";
					};#if
				} else {
					if ($module_name!="") {
						$QC[$QC_id]["QCmodule"][$module_name].=implode("\t",$l)."\n";
					};#if
				};#if
				#
			};#foreach

		};#if

		if ($QC_metrics=="HsMetrics") {

			# Found lines of header and values
			$line_index_names=6;
			$line_index_values=7;
			foreach ($QC[$QC_id]["QCarray"] as $line_index=>$line_content) {
				#print "LINE CONTENT".$line_content[0]."<pre>"; print_r($line_content);  print "</pre><BR>";
				if ($line_content[0][0] != "#" && $line_content[0] != "") { # first line without #
					$line_index_names=$line_index;
					$line_index_values=$line_index+1;
					break;
				};#if
			};#foreach

			$QCFile_array_names=$QC[$QC_id]["QCarray"][$line_index_names];
			$QCFile_array_values=$QC[$QC_id]["QCarray"][$line_index_values];
			#print "$vcf_aligner<pre>"; print_r($QCFile_array_names); print_r($QCFile_array_values);  print "</pre><BR>";
			foreach ($QCFile_array_names as $variable_key=>$variable_name) {
				$variable_value=$QCFile_array_values[$variable_key];
				$QCMetrics_array[$vcf_aligner][$variable_name]=$variable_value;
				#print "$variable_name]=$variable_value<BR>";
			};#foreach
		};#if
		#echo "<PRE>"; print_r($QCMetrics_array); echo "</PRE>";

		$summary.= "<tr>"
			."<td>$QC_id</td>"
			."<td>$QC_metrics</td>"
		."</tr>";
	};#while
	#print_r($QC);
	return array("QC"=>$QC,"QCFASTQC_id"=>$QCFASTQC_id,"QCMetrics_array"=>$QCMetrics_array,"warnings"=>$warnings,"module_for_warnings"=>$module_for_warnings,"QCarray"=>$QCarray);

};#function QC

function merge_filters($filter1,$filter2,$filter_array) {
	# if filter 1 is based on filter 2, ad filter 2 not based on filter1
	#if (isset($filter_array[$filter1]["basedon"]) && $filter_array[$filter2]["basedon"]!=$filter1) {
		# if filter 2 is based on another filter
		#echo "$filter1,$filter2,".$filter_array[$filter1]["basedon"]."<BR>";
		if (isset($filter_array[$filter2]["basedon"]) && isset($filter_array[$filter_array[$filter2]["basedon"]]) && $filter_array[$filter2]["basedon"]!=$filter1) 			{
			$filter_array=merge_filters($filter2,$filter_array[$filter2]["basedon"],$filter_array);
		};#if
		if (isset($filter_array[$filter1]["basedon"]) && isset($filter_array[$filter_array[$filter1]["basedon"]]) && $filter_array[$filter1]["basedon"]!=$filter2) 			{
			$filter_array=merge_filters($filter1,$filter_array[$filter1]["basedon"],$filter_array);
		};#if
		# merge arrays
		$filter_array[$filter1]=array_merge($filter_array[$filter2],$filter_array[$filter1]);
	#};#if
	return $filter_array;
}

function basedon_filters($filter_array) {

	foreach ($filter_array as $filter => $filter_criteria) {
		#echo "FILTER: $filter<BR>";
		if (isset($filter_array[$filter]["basedon"]) && isset($filter_array[$filter_array[$filter]["basedon"]])) {
			$filter_array=merge_filters($filter,$filter_array[$filter]["basedon"],$filter_array);
		};#if
	};#foreach
	return $filter_array;

};#function merge_filters

function load_filter($filters,$filter_array,$filter_option_array=array()) {

	foreach ($filters as $filter) {
		foreach ($filter_array[$filter] as $annotation => $annotation_filter) {
			#echo "$annotation => $annotation_filter<BR>";
			if (is_array($annotation_filter)) {
				foreach ($annotation_filter as $annotation_filter_key => $annotation_filter_value) {
					#echo "&nbsp;&nbsp;&nbsp;&nbsp;$annotation_filter_key => $annotation_filter_value<BR>";
					#$filter_option_array[$annotation]="$annotation:$annotation_filter_value";
					$filter_option_array[]="$annotation:$annotation_filter_value";
				};#foreach
			} else {
				#$filter_option_array[$annotation]="$annotation:$annotation_filter";
				$filter_option_array[]="$annotation:$annotation_filter";
			};#if
		};#foreach
	};#foreach
	#print_r($filter_option_array);
	return $filter_option_array;

};#function load_filter




function write_php_ini($array, $file)
{
    $res = array();
    foreach($array as $key => $val)
    {
        if(is_array($val))
        {
            $res[] = "[$key]";
            foreach($val as $skey => $sval) {
		if (is_array($sval) ) {
			foreach($sval as $sskey => $ssval) {
				$res[] = "$skey [] = ".(is_numeric($ssval) ? $ssval : '"'.$ssval.'"');
			};#foreach
		} else {
			$res[] = "$skey = ".(is_numeric($sval) ? $sval : '"'.$sval.'"');
		};#if
	    };#foreach
        }
        else $res[] = "$key = ".(is_numeric($val) ? $val : '"'.$val.'"');
    }
    safefilerewrite($file, implode("\r\n", $res));
}

function safefilerewrite($fileName, $dataToSave)
{
  if ($fp = fopen($fileName, 'w'))
  {
      $startTime = microtime();
      do
      { $canWrite = flock($fp, LOCK_EX);
         // If lock not obtained sleep for 0 - 100 milliseconds, to avoid collision and CPU load
         if(!$canWrite) usleep(round(rand(0, 100)*1000));
      } while ((!$canWrite)and((microtime()-$startTime) < 1000));

      //file was locked so now we can store information
      if ($canWrite)
      { fwrite($fp, $dataToSave);
       flock($fp, LOCK_UN);
      }
      fclose($fp);
     }

}

function ScanDirectory($Directory){
  	$output="";
  	$MyDirectory = opendir($Directory); # or die('Erreur');
	while($Entry = @readdir($MyDirectory)) {
		if(is_dir($Directory.'/'.$Entry)&& $Entry != '.' && $Entry != '..') {
                         $output.= '<ul>'.$Directory;
			$output.=ScanDirectory($Directory.'/'.$Entry);
                        $output.= '</ul>';

		}
		else {
			$output.= '<li>'.$Entry.'</li>';
                }
	}
  	closedir($MyDirectory);
	return $output;
}

function ListRUN($Directory,$ListSample,$ListAligner,$ListCaller,$ListFiles){
  	$output="";
  	$MyDirectory = opendir($Directory); # or die('Erreur');
	while($Entry = @readdir($MyDirectory)) {
		if(is_dir($Directory.'/'.$Entry)&& $Entry != '.' && $Entry != '..') {
            		$output.= '<ul>'.$Entry;
			if ($ListSample) {
				$output.= '<ul>';
				$output.=ListSAMPLE($Directory.'/'.$Entry,$ListAligner,$ListCaller,$ListFiles);
				$output.= '</ul>';
				#return $output;
			};#if
			$output.= '</ul>';
		}
		else {
			#echo '<li>'.$Entry.'</li>';
        };#if
	}

	closedir($MyDirectory);
	return $output;
}

function ListSAMPLE($Directory,$RUN,$SAMPLES=array()){
  	require "connect.php";
	require "config.php";

	$output="";
	$RUNLink="run.php?action_todo=show_run&run=$RUN";
	$RUNList.="<DIV class='brick run' onclick=\"javascript:document.location='$RUNLink'\">$RUN</DIV>";

	$query="SELECT sample.ref AS sample_ref
			FROM sample
			INNER JOIN run ON (run.id=sample.run_id)
			WHERE run.ref='$RUN'
			";
		#print "$query<BR>";
		$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		while ($row = mysql_fetch_assoc($result)) {
			$SAMPLE=$row["sample_ref"];
			#print_r($SAMPLES);
			if (empty($SAMPLES) || in_array($SAMPLE,$SAMPLES)) {
			$SAMPLELink="sample.php?action_todo=show_sample&run=$RUN&sample=$SAMPLE";
			$SAMPLEList.="<DIV class='brick sample' onclick=\"javascript:document.location='$SAMPLELink'\" title='RUN: $RUN'>$SAMPLE</div>";
			};#if
		};#while


	/*
	$MyDirectory = opendir($Directory); # or die('Erreur');
	while($SAMPLE = @readdir($MyDirectory)) {

		if(is_dir($Directory.'/'.$SAMPLE) && $SAMPLE != '.' && $SAMPLE != '..') {

			if (empty($SAMPLES) || in_array($SAMPLE,$SAMPLES)) {
			$SAMPLELink="sample.php?action_todo=show_sample&run=$RUN&sample=$SAMPLE";
			$SAMPLEList.="<DIV class='brick sample' onclick=\"javascript:document.location='$SAMPLELink'\" title='RUN: $RUN'>$SAMPLE</div>";
			};#if
        	};#if
	};#while
	closedir($MyDirectory);
	*/
	$break="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";
	return $RUNList.$SAMPLEList.$break;
	#return $output;
}

function ListSAMPLE_old($Directory,$RUN,$SAMPLES=array()){
  	$output="";
	$MyDirectory = opendir($Directory); # or die('Erreur');
	$RUNLink="run.php?action_todo=show_run&run=$RUN";
	$RUNList.="<DIV class='brick run' onclick=\"javascript:document.location='$RUNLink'\">$RUN</DIV>";
	while($SAMPLE = @readdir($MyDirectory)) {

		if(is_dir($Directory.'/'.$SAMPLE) && $SAMPLE != '.' && $SAMPLE != '..') {

			if (empty($SAMPLES) || in_array($SAMPLE,$SAMPLES)) {
			$SAMPLELink="sample.php?action_todo=show_sample&run=$RUN&sample=$SAMPLE";
			$SAMPLEList.="<DIV class='brick sample' onclick=\"javascript:document.location='$SAMPLELink'\" title='RUN: $RUN'>$SAMPLE</div>";
			};#if
        	};#if
	};#while
	$break="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";
	closedir($MyDirectory);
	return $RUNList.$SAMPLEList.$break;
	#return $output;
}

function ListAligner($Directory,$ListCaller,$ListFiles){
  	$output="";
	echo $Directory;
  	$MyDirectory = opendir($Directory); # or die('Erreur');
	while($Entry = @readdir($MyDirectory)) {
		if(is_dir($Directory.'/'.$Entry)&& $Entry != '.' && $Entry != '..') {
            		$output.= '<ul>'.$Entry;
			if ($ListCaller) {
				$output.= '<ul>';
				$output.=ListCaller($Directory.'/'.$Entry,$ListFiles);
				$output.= '</ul>';
			};#if
			$output.= '</ul>';
		}
		else {
			#echo '<li>'.$Entry.'</li>';
        };#if
	}
	closedir($MyDirectory);
	return $output;
}

function recurse_copy($src,$dst) {
    $dir = opendir($src);
    #print "SRC $src DEST $dst<BR>\n";
    @mkdir($dst);
    while(false !== ( $file = readdir($dir)) ) {
	#print "$file   ";
        if (( $file != '.' ) && ( $file != '..' )) {
            if ( is_dir($src . '/' . $file) ) {
		#print " is DIR";
                recurse_copy($src . '/' . $file,$dst . '/' . $file);
            }
            else {
		#print ": ".$src . '/' . $file." => ".$dst . '/' . $file." ";
                if (copy($src . '/' . $file,$dst . '/' . $file)) {
			#print " COPIED";
		} else {
			#print " NOT COPIED (error)";
		};#if
            }
        }
	#print "<BR>\n";
    }
    closedir($dir);
}

function ShowSAMPLE_VariantList($Directory,$RUN,$SAMPLE,$Aligner="MiSeq-Aligner",$Caller="MiSeq-Caller",$Annotator="trakxs",$vcf_file=""){
  	global $dir_analysis, $dir_bin;
	global $hardfiltering, $orderby, $ascdesc, $limit, $annotation_list, $filters, $filter_array, $annotation_list_default;
	global $tabfileforannoationlist;
	if (!is_array($filters)) {
		$filters[0]=$filters;
	};#if
	$filter_option_array=array();
	$filter_option="";
	$filter_option_array=load_filter($filters,$filter_array,$filter_option_array);
	if (!empty($filter_option_array)) {
		$filter_option="--filters='".join(",",$filter_option_array)."'";
		# criteria filter file [TODO]
		#print "$filter_option<BR>";
	};#if
	$VCFFile="$SAMPLE.annotated.vcf";
	if ($Caller=="union") {
		$Caller="";
		$Annotator="";
		$VCFFile="$SAMPLE.annotated.union.vcf";
	};#if

	#$html_file="$Directory/$RlUN/$SAMPLE/DATA/$Aligner/$Caller/$Annotator/$SAMPLE.html";
	if ($vcf_file == "") {
		$vcf_file="$Directory/$RUN/$SAMPLE/DATA/$Aligner/$Caller/$Annotator/$VCFFile";
	};#if
	#print $vcf_file;
	if (!is_file($vcf_file)) {
		return "<DIV class='html'>No VCF file. Please contact your favorite bioinformatician to perform the analysis<BR><BR></DIV>";
	};#if
	#$annotation_list="FilterScore,FilterFlag,genesymbol,hgvs";
	$annotation_list_option=join(",",$annotation_list);
	#$orderby="FilterScore";
	$hardfiltering ? $remove_filtered="--remove_filtered" : $remove_filtered="--noremove_filtered";
	#echo $html_file;
  	#if (is_file($vcf_file)) {
	$tmpfnamebase = tempnam("tmp/", "variant_file_");

	$limits=array("HTML"=>$limit);
	#$verbose_debug=" --verbose --debug ";

	# list of annotations
	$tabfileforannoationlist="$tmpfnamebase.listannotation";
	#$command="$dir_bin/scripts/VCF.pl --input_file=$vcf_file --output_file=$tabfileforannoationlist --annotation_list='$annotation_list_default,ALL' --output_format=TAB --noremove_filtered --verbose --debug";
	#$command="$dir_bin/scripts/VCF.pl --input_file=$vcf_file --output_file=$tabfileforannoationlist --annotation_list='$annotation_list,$annotation_list_default$annotation_list_mandatoy,ALL' --output_format=TAB --noremove_filtered --verbose --debug";
	$limit_to_detect_annotations=" --limit=1000" ;
	#print "$annotation_list_option |, $annotation_list_default | $annotation_list_mandatoy,all<BR>";
	#print "$vcf_file\n";
	$command="$dir_bin/scripts/VCF.pl --input_file='$vcf_file' --output_file='$tabfileforannoationlist' --annotation_list='$annotation_list_option,$annotation_list_default$annotation_list_mandatoy,all' --output_format=TAB --noremove_filtered $verbose_debug $limit_to_detect_annotations ";
	#$command="$dir_bin/scripts/VCF.pl --input_file='$vcf_file' --output_file='$tabfileforannoationlist' --output_format=TAB --noremove_filtered $verbose_debug $limit_to_detect_annotations ";
	#echo "$command<BR>\n";
	$output_exec = shell_exec($command);
	#print "<pre>$output_exec</pre>";
	#print $_SERVER['PHP_SELF']."<BR>\n"; #name of the PHP script in the URL
	#print $_SERVER['SCRIPT_FILENAME']."<BR>\n"; #name of the PHP script
	#print "$tmpfnamebase<BR>\n";


	foreach (array("HTML"=>"div","HTMLTABLE"=>"html","VCF"=>"vcf","TAB"=>"txt") as $export=>$extension) {
	#foreach (array("HTML"=>"div") as $export=>$extension) {

		#$tmpfname = "$tmpfnamebase.".($export=="TAB"?"txt":strtolower($export));
		#$tmpfname = "$tmpfnamebase.".($export=="HTMLTABLE"?"html":strtolower($export));
		$limit_option=" ";
		if (isset($limits[$export])) {
			$limit_option=" --limit=".$limits[$export]." ";
		};#if
		$tmpfname = "$tmpfnamebase.$extension";
		#echo "$export $tmpfname<BR>";
		#$command="$dir_bin/scripts/VCF.pl --input_file='$vcf_file' --output_file='$tmpfname' --annotation_list='$annotation_list_option' --output_format=$export $remove_filtered $filter_option --order_by=$orderby,$ascdesc $limit_option $verbose_debug ";
		$command="$dir_bin/scripts/VCF.pl --input_file='$vcf_file' --output_file='$tmpfname' --annotation_list='$annotation_list_option' --output_format=$export $remove_filtered $filter_option --order_by=$orderby,$ascdesc $limit_option $verbose_debug ";
		echo "$command<BR>\n";
		$output_exec = shell_exec($command);
		print "<pre>$output_exec</pre>";
	};#foreach

	# OUTPUT
	$ExportFiles.="<DIV class='brick file' onclick=\"javascript:var newWin = window.open('tmp/".basename($tmpfnamebase).".vcf','popup','menubar=yes, status=no, scrollbars=no, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\">VCF</DIV> ";
	$ExportFiles.="<DIV class='brick file' onclick=\"javascript:var newWin = window.open('tmp/".basename($tmpfnamebase).".txt','popup','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\">TRAKXS</DIV> ";
	$ExportFiles.="<DIV class='brick file' onclick=\"javascript:var newWin = window.open('tmp/".basename($tmpfnamebase).".html','popup','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\">HTML</DIV> ";
	#$output.="<DIV class='brick file'><A HREF='http://www.broadinstitute.org/igv/projects/current/igv.php?file=http://192.168.59.62/CPS2876II06.bam,http://192.168.59.62/CPS2876II06.annotated.vcf&genome=hg19&merge=true' target='BAM'>IGV</A></DIV>";
	$files="\\\\\\\\172.31.2.41\\\\ngs\\\\analysis\\\\$RUN\\\\$SAMPLE\\\\DATA\\\\$Aligner\\\\$SAMPLE.bam,\\\\\\\\172.31.2.41\\\\ngs\\\\analysis\\\\$RUN\\\\$SAMPLE\\\\DATA\\\\$Aligner\\\\$Caller\\\\trakxs\\\\$SAMPLE.annotated.vcf";
	#$output.="<DIV class='brick file'><A HREF='http://www.broadinstitute.org/igv/projects/current/igv.php?file=$files&genome=hg19&merge=true' target='BAM' onclick=\"javascript:alert('Be sure you already open the \\\\\\\\172.31.2.41\\\\ngs\\\\ folder in a windows explorer!!!')\">IGV</A></DIV>";
	$output.=$ExportFiles;
	$output.="<DIV class='brick file' onclick=\"javascript:var newWin = alert('Be sure you already open the \\\\\\\\172.31.2.41\\\\ngs\\\\ folder in a windows explorer!!! (if you are allowed...)'); window.open('http://www.broadinstitute.org/igv/projects/current/igv.php?file=$files&genome=hg19&merge=true','popup','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\">IGV [under development...]</DIV>";
	$output.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";
	#$output.="<DIV class='html' style='clear:left;'>".file_get_contents("$tmpfnamebase.htmltable")."<BR></DIV>";
	#$output.="<DIV class='' style='clear:left;'>".file_get_contents("$tmpfnamebase.html")."<BR></DIV>";
	$output.="<DIV>".file_get_contents("$tmpfnamebase.div")."<DIV>";
	#$output.="".file_get_contents("$tmpfnamebase.html")."";

	return array($output,$tmpfnamebase,$ExportFiles);
}

function implode_with_key($assoc, $inglue = '>', $outglue = ',') {
    $return = '';

    foreach ($assoc as $tk => $tv) {
        $return .= $outglue . $tk . $inglue . $tv;
    }

    return substr($return, strlen($outglue));
}

function VCFLinetoArray ($variant_line,$variant_header) {
# Structure of output
# Array with keys: CHROM POS... and INFO SAMPLES
#    Array INFO
#    Array SAMPLES an Array of qual
	$variant_line_array=explode("\t",$variant_line);
	$variant_array= array();
	$variant_array["CHROM"]=$variant_line_array[0];
	$variant_array["POS"]=$variant_line_array[1];
	$variant_array["ID"]=$variant_line_array[2];
	$variant_array["REF"]=$variant_line_array[3];
	$variant_array["ALT"]=$variant_line_array[4];
	$variant_array["QUAL"]=$variant_line_array[5];
	$variant_array["FILTER"]=$variant_line_array[6];
	$variant_array["INFO"]=array();
	$variant_array["FORMAT"]=$variant_line_array[8];
	$variant_array["SAMPLES"]=array();

	#INFO
	foreach (explode(";",$variant_line_array[7]) as $annotation_nb=>$annotation) {
		#print "$annotation_nb=>$annotation<BR>";
		$annotation=trim($annotation);
		$annotation_array=explode("=",$annotation,2);
		$annotation_source=trim($annotation_array[0]);
		$annotation_value=trim($annotation_array[1]);
		$variant_array["INFO"][$annotation_source]=$annotation_value;
	};#foreach

	#print "$variant_line<BR>";

	#SAMPLES
	$format_array=explode(":",$variant_array["FORMAT"]);
	#print "<PRE>"; print_r($variant_array["FORMAT"]); print "</PRE>";
	foreach ($variant_line_array as $col_nb=>$col) {
		#print "$col_nb <BR>";
		if ($col_nb>=9) {
			#print "$col OK <BR>";
			$sample_name=$variant_header[$col_nb];
			#print "$col_nb $sample_name<BR>";
			#$variant_array["SAMPLES"][$sample_name]["VCFCOL"]=$col;
			$col_array=explode(":",$col);
			foreach ($col_array as $col_qual_nb=>$col_qual) {
				#print "   $col_qual_nb=>$col_qual<BR>";
				$col_qual_name=$col_qual_nb;
				if ($format_array[$col_qual_nb]!="") {
					$col_qual_name=$format_array[$col_qual_nb];
				};#if
				#$variant_array["SAMPLES"][$sample_name][$col_qual_name]=$col_qual;
				#print "   [$sample_name]|$col_qual_name]=$col_qual<BR>";
				$variant_array["SAMPLES"][$sample_name][$col_qual_name]=$col_qual;
			};#foreach
		};#if
	};#foreach
	ksort($variant_array["SAMPLES"]);
	#print "<PRE>"; print_r($variant_array["SAMPLES"]); print "</PRE>";


	# snpEff
	#print $variant_array["INFO"]["snpeff"]; print "<BR>";
	#if (array_key_exists ( "snpeff" , $variant_array["INFO"] )) { $snpeff=$variant_array["INFO"]["snpeff"];}
	if (array_key_exists ( "ANN" , $variant_array["INFO"] )) {
		$snpeff=$variant_array["INFO"]["ANN"];
	} elseif (array_key_exists ( "snpeff" , $variant_array["INFO"] )) {
		$snpeff=$variant_array["INFO"]["snpeff"];
	};
	$snpeff_symbol=explode ( "|" , explode ( "," , $snpeff )[0])[3];
	$snpeff_transcript=explode ( ":" , explode ( "|" , explode ( "," , $snpeff )[0])[6])[0];
	$snpeff_cnomen=explode ( ":" , explode ( "|" , explode ( "," , $snpeff )[0])[9])[0];
	$snpeff_pnomen=explode ( ":" , explode ( "|" , explode ( "," , $snpeff )[0])[10])[0];
	$snpeff_hgvs="$snpeff_symbol:$snpeff_transcript:$snpeff_cnomen";
	if ($snpeff_pnomen != "") { $snpeff_hgvs="$snpeff_hgvs:$snpeff_pnomen"; };
	#print "TRUC: "; print_r($snpeff_symbol); print " | ";  print_r($snpeff_hgvs); print "<BR>";
	#print "<BR>"; print "<BR>";


	# Processing Annotation
	# CHROMPOS
	$variant_array["CHROMPOS"]=$variant_array["CHROM"].":".$variant_array["POS"]."|".$variant_array["REF"].">".$variant_array["ALT"];
	# Position
	$variant_array["POSITION"]=$variant_array["CHROM"].":".$variant_array["POS"]."-".$variant_array["POS"];
	# Position
	$variant_array["REGION"]=str_replace("chr","",$variant_array["CHROM"]).":".($variant_array["POS"]-20)."-".($variant_array["POS"]+20);
	# Variant_id
	#$variant_array["VARIANT_ID"]=$variant_array["INFO"]["variant_id"];
	$variant_array["VARIANT_ID"]=$variant_array["ID"];


	# Symbol
	if ($variant_array["INFO"]["Symbol"]!="") {
		$Symbol_split=explode("\(",$variant_array["INFO"]["Symbol"]);
	} elseif ($variant_array["INFO"]["geneSymbol"]!="") {
		$Symbol_split=explode("\(",$variant_array["INFO"]["geneSymbol"]);
	} elseif ($variant_array["INFO"]["GI"]!="") {
		$Symbol_split=explode(",",$variant_array["INFO"]["GI"]);
	} elseif ($snpeff_symbol!="") {
		$Symbol_split=explode(",",$snpeff_symbol);
	};#if
	$Symbol=$Symbol_split[0];
	$Symbol_hgvs=str_replace("\)","",$Symbol_split[1]);
	#$variant_array["SYMBOL"]=$variant_array["INFO"]["Symbol"];
	$variant_array["SYMBOL"]=$Symbol;


	# HGVS
	if ($variant_array["INFO"]["hgvs"]!="") {
		$hgvs_input=$variant_array["INFO"]["hgvs"];
	} elseif ($snpeff_hgvs!="") {
		$hgvs_input=$snpeff_hgvs;
	};

	#$hgvs_full_split=explode(",",$variant_array["INFO"]["hgvs"]);
	$hgvs_full_split=explode(",",$hgvs_input);

	$hgvs_full=(trim($hgvs_full_split[0])!="")?$hgvs_full_split[0]:$Symbol_hgvs;
	$hgvs_split=explode(",",$hgvs_full);
	$hgvs=$hgvs_split[0];
	#print($hgvs);
	#print_r($hgvs_full_split);
	$variant_array["HGVS"]=$hgvs;
	$variant_array["HGVS"]=$variant_array["INFO"]["hgvs"];

	#PNOMEN
	$hgvs_pnomen="";
	foreach (explode(":",explode(",",$variant_array["HGVS"])[0]) as $key=>$hgvs_info) {
		#print $hgvs_info;
		$hgvs_pnomen=(explode(".",$hgvs_info)[0]=="p")?$hgvs_info:"";
	};
	$variant_array["PNOMEN"]=$hgvs_pnomen;
	#print "<BR>$hgvs_pnomen<BR>";
	#$pnomen=
	#print $variant_array["HGVS"];

	# SPLICING REGION SPECIFICITY
	#  Ensembl
	$ensembl_split=explode("\(",$variant_array["INFO"]["Ensembl"]);
	#print " ".$variant_array["INFO"]["Ensembl"]." -> $ensembl_split <BR>";
	$variant_array["ENSEMBL"]=$ensembl_split[0];
	#  COSMIC
	#$cosmic_split=split(array(":"),$variant_array["INFO"]["COSMIC"]);
	#print " ".$variant_array["INFO"]["COSMIC"]." -> ".$cosmic_split[0]." <BR>";
	$cosmic_preg_split=preg_split("/[:,]/",$variant_array["INFO"]["COSMIC"]);
	#print " ".$variant_array["INFO"]["COSMIC"]." -> "; print_r($cosmic_preg_split); print " <BR>";
	$variant_array["COSMIC"]=$cosmic_preg_split[1];
	# HGVS Ensembl
	$hgvs_ensembl_split=explode(",",$variant_array["INFO"]["outcomeinfos"]);
	$variant_array["HGVS_ENSEMBL"]=$hgvs_ensembl_split[0];
	# filterScore
	#print_r($variant_array["INFO"]);
	if ($variant_array["INFO"]["PZScore"]!="") {
		$variant_array["FILTER_SCORE"]=$variant_array["INFO"]["PZScore"];
	} else {
		$variant_array["FILTER_SCORE"]=$variant_array["INFO"]["FilterScore"];
	};#if
	# filterFlag
	$variant_array["FILTER_FLAG"]=$variant_array["INFO"]["PZFlag"];
	# filterComment
	$variant_array["FILTER_COMMENT"]=$variant_array["INFO"]["PZComment"];
	# filterInfos
	$variant_array["FILTER_INFOS"]=$variant_array["INFO"]["PZInfos"];
	# location
	$variant_array["LOCATION"]=$variant_array["INFO"]["location"];
	# outcome
	$variant_array["OUTCOME"]=$variant_array["INFO"]["outcome"];
	# location
	$variant_array["DBSNP"]=$variant_array["INFO"]["dbSNP"];
	# outcome
	$variant_array["DBSNP_NONFLAGGED"]=$variant_array["INFO"]["dbSNPNonFlagged"];
	# Quality
	#GQ=99 | DP=975 | AD=14,2430 | GT=1/1
	#$variant_array["GQ"]=$variant_array["INFO"]["GQ"];
	#$variant_array["DP"]=$variant_array["INFO"]["DP"];
	#$variant_array["AF"]=$variant_array["INFO"]["AlleleFrequency"];
	#$variant_array["AD"]=$variant_array["INFO"]["AD"];
	#$variant_array["GT"]=$variant_array["INFO"]["GT"];

	# ALLELE FREQUENCY calculation
	foreach ($variant_array["SAMPLES"] as $sample_name=>$qualities) {
		$AlleleFrequency="";
		if ($qualities["VAF"]!="" && $qualities["VAF"]!="." ) {
			$AlleleFrequency=$qualities["VAF"];
		} elseif ($qualities["FA"]!="" && $qualities["FA"]!="." ) {
			$AlleleFrequency=$qualities["FA"];
		} elseif ($qualities["FREQ"]!="" && $qualities["FREQ"]!="." ) {
			if (strpos($qualities["FREQ"],"%")) {
				$AlleleFrequency=str_replace("%","",$qualities["FREQ"])/100;
			};#if
		} elseif ($qualities["AD"]) {
			#print $qualities["AD"]."<BR>";
			$qualities_array=explode(",",$qualities["AD"]);
			$depth_total=0;
			$depth_ref=0;
			$depth_alt=0;
			foreach ($qualities_array as $key=>$depth) {
				$depth_total+=$depth;
				if ($key==0) {
					$depth_ref=$depth;
				} else {
					$depth_alt+=$depth;
				};#if
			};#foreach
			#print "$depth_total $depth_ref $depth_alt<BR>";
			$AlleleFrequency=$depth_alt/$depth_total;

		};#if
		$variant_array["SAMPLES"][$sample_name]["AlleleFrequency"]=number_format($AlleleFrequency,4);
	};#foreach


	# GENERAL QUALITY VALUES
	$nb_sample=0;
	$GQ_sum=0;
	$GQ_min="";
	$GQ_nb=0;
	$DP_sum=0;
	$DP_min=0;
	$DP_nb=0;
	$AlleleFrequency_sum=0;
	$AlleleFrequency_min=0;
	$AlleleFrequency_nb=0;
	foreach ($variant_array["SAMPLES"] as $sample_name=>$qualities) {
		$nb_sample++;
		$GQ_sum+=$qualities["GQ"];
		#$GQ_min=($qualities["GQ"]<$GQ_min && $GQ_min!="")?$qualities["GQ"]:$GQ_min;
		$GQ_min=($qualities["GQ"]!="" && $qualities["GQ"]!="." && ($qualities["GQ"]<$GQ_min || $GQ_min==""))?$qualities["GQ"]:$GQ_min;
		$GQ_nb+=($qualities["GQ"]!="");
		$DP_sum+=$qualities["DP"];
		$DP_min=($qualities["DP"]!="" && $qualities["DP"]!="." && ($qualities["DP"]<$DP_min || $DP_min==""))?$qualities["DP"]:$DP_min;
		$DP_nb+=($qualities["DP"]!="");
		$AlleleFrequency_sum+=$qualities["AlleleFrequency"];
		$AlleleFrequency_min=($qualities["AlleleFrequency"]<$AlleleFrequency_min || $AlleleFrequency_min=="")?$qualities["AlleleFrequency"]:$AlleleFrequency_min;
		$AlleleFrequency_nb+=($qualities["AlleleFrequency"]!="");
	};#foreach
	$GQ_avg=$GQ_sum/$GQ_nb;
	$DP_avg=$DP_sum/$DP_nb;
	$AlleleFrequency_avg=$AlleleFrequency_sum/$AlleleFrequency_nb;

	if ($variant_array["INFO"]["VAF_median"]!="") {
		$AlleleFrequency_avg=$variant_array["INFO"]["VAF_median"];
	}; # if

	$variant_array["GQ"]=$GQ_min;
	$variant_array["DP"]=$DP_min;
	$variant_array["AlleleFrequency"]=$AlleleFrequency_avg;
	$variant_array["AD"]=$variant_array["SAMPLES"][""]["AD"];
	$variant_array["GT"]=$variant_array["SAMPLES"][""]["GT"];

	foreach ($variant_array["INFO"] as $info=>$info_value) {
		if ($variant_array[strtoupper($info)]=="") {
			$variant_array[strtoupper($info)]=$info_value;
		};#if
	};#foreach



	# Return
	return $variant_array;
}


function find_hgvs ($hgvs_list,$gene_symbol="",$user_groups="",$customNM_infos="") {
	#require "connect.php";
	require "config.php";

	#print "$hgvs_list,$gene_symbol<BR>";
	#print_r($user_groups); print_r($customNM_infos);

	# $user_groups and $customNM_infos in config.php
	if ($gene_symbol=="") {
		$hgvs_split=explode(",",$hgvs_list);
		$hgvs_split2=explode(":",$hgvs_split[0]);
		$gene_symbol=$hgvs_split2[0];
	};#if
	#print $gene_symbol;

	# CustomHGVS
	$customNM_V="";
	#$customHGVS=$hgvs_list;
	$customHGVS="";
	if (isset($user_groups)) {
		foreach ($user_groups as $group_ref) {
			#print "$group_ref<BR>";
			$hgvs_full_split=explode(",",$hgvs_list);
			#print "customNM_V=$customNM_V<BR>";
			if (isset($customNM_infos[$group_ref][$gene_symbol])) {
				$customNM_array_V=$customNM_infos[$group_ref][$gene_symbol]; # ARRAY!!!
				foreach ($hgvs_full_split as $anhgvs) {
					if ($anhgvs!="") {
						foreach ($customNM_array_V as $customNM_V=>$customNM_infos) {
							#print "$anhgvs $customNM_V<BR>";
							if (preg_match("/$customNM_V/i",$anhgvs)) {
								$sep=($customHGVS!="")?",":"";
								$customHGVS.="$sep$anhgvs";
								#print "#customHGVS=$customHGVS<BR>";
								#break 2;
							};#if
						};#foreach
					};#if
				};#foreach
			} else {
				#$hgvs_full=(trim($hgvs_full_split[0])!="")?$hgvs_full_split[0]:$Symbol_hgvs;
				#$hgvs_split=split(",",$hgvs_full);
				#$hgvs=$hgvs_split[0];
				#$customHGVS=$hgvs;
				#break;
			};#if
		};#foreach
	};#if
	if ($customHGVS=="") {
		#$customHGVS=$hgvs_list;
		$hgvs_full=(trim($hgvs_full_split[0])!="")?$hgvs_full_split[0]:$Symbol_hgvs;
		$hgvs_split=explode(",",$hgvs_full);
		$hgvs=$hgvs_split[0];
		$customHGVS=$hgvs;
	};#if
	return $customHGVS;

}#function

function VCFtoHTML($Variants_VCFContent,$sample_id=0,$annotation_list=array("ALL"),$options=array()) {

	require "connect.php";
	require "config.php";

	if ($project_id=="") {
		$project_id=1;
	};#if

	# found positions
	$all_positions=array();
	foreach (split("\n",$Variants_VCFContent) as $line_nb=>$line) {
		$line_split=explode("\t",$line);
		if ($line_split[0][0] != "#") {
			#echo $line_split[0].":".$line_split[1]."<BR>";
			$all_positions[$line_split[0]][$line_split[1]]=1;
		}#if
	};#foreach

	# SNAPSHOT
	if ($options["vcf_ids"]!="" || $options["report"]!="") {
		$variant_files=$options["vcf_ids"];
		#print_r($options["vcf_ids"]);
		if ($options["vcf_ids"]!="") {
			$where_clause="vcf.id IN ($variant_files)";
		} elseif ($options["report"]!="") {
			$where_clause="vcf.sample_id=$sample_id";
		};#if

		# SNAPSHOT
		$query="SELECT snapshot.id AS snapshot_id, snapshot.ref AS snapshot_ref, snapshot.chrom AS snapshot_chrom, snapshot.pos AS snapshot_pos, vcf.aligner AS vcf_aligner, vcf.bam, vcf.ref
			FROM vcf
			INNER JOIN snapshot_vcf ON (snapshot_vcf.vcf_id=vcf.id)
			INNER JOIN snapshot ON (snapshot.id=snapshot_vcf.snapshot_id)
			WHERE $where_clause
			";
		#print "$query<BR>";

		$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));

		while ($row = mysql_fetch_assoc($result)) {
			$snapshot_id=$row["snapshot_id"];
			$snapshot_ref=$row["snapshot_ref"];
			$snapshot_chrom=$row["snapshot_chrom"];
			$snapshot_pos=$row["snapshot_pos"];
			$vcf_aligner=$row["vcf_aligner"];
			$snapshot_array[$snapshot_chrom][$snapshot_pos][$vcf_aligner]=$snapshot_id;
			#print "<pre>"; print_r($row); print "</pre>";
			#$URL="blob_image.php?id=$snapshot_id";
			#echo "<span class=' ' onclick=\"javascript: var newWin = window.open('$URL','SNAPSHOT','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=650, height=800, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>$snapshot_ref</span><BR>";
		};#while

	};#if


	# QC
	if (0 && $options["vcf_ids"]!="" || $options["report"]!="") {
		#echo "QC_depthbed";
		$variant_files=$options["vcf_ids"];
		#print_r($options["vcf_ids"]);
		if ($options["vcf_ids"]!="") {
			$where_clause="vcf.id IN ($variant_files)";
		} elseif ($options["report"]!="") {
			$where_clause="vcf.sample_id=$sample_id";
		};#if

		# QC
		$query="SELECT distinct QC.id, QC.ref, QC.QC
			FROM vcf
			INNER JOIN QC_association ON (QC_association.object_id=vcf.id AND QC_association.object_type='vcf')
			INNER JOIN QC ON (QC.id=QC_association.QC_id)
			WHERE QC.metrics='depthbed'
			  AND QC.type='AlignedBAM'
			  AND $where_clause";
		#print "$query<BR>";

		$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		while ($row = mysql_fetch_assoc($result)) {
			#print "<PRE>"; print_r($row);  print "</PRE>";
			$QC_id=$row["id"];
			$QC_ref=$row["ref"];
			$QC_QC=$row["QC"];
			$QC_QC_split=explode("\n",$QC_QC);
			#print "<PRE>"; print_r($row);  print "</PRE>";

			foreach ($QC_QC_split as $line_nb=>$QC_QC_line) {
				#print "$line_nb=>$QC_QC_line<BR>\n";
				$QC_QC_line_split=explode("\t",$QC_QC_line);
				$QC_QC_line_chrom=str_replace("chr","",$QC_QC_line_split[0]);
				$QC_QC_line_pos=$QC_QC_line_split[1];
				if ($all_positions[$QC_QC_line_split[0]][$QC_QC_line_pos] > 0 ) {
					print "# Depth position ".$QC_QC_line_split[0].$QC_QC_line_pos."<BR>";
					$QC_QC_line_depth=$QC_QC_line_split[2];
					$depthbed_array[$QC_QC_line_chrom][$QC_QC_line_pos][$QC_id]=$QC_QC_line_depth;
					$depthbed_QCidref_array[$QC_id]=$QC_ref;
				};#if
			};#foreach
			#print "<pre>"; print_r($row); print "</pre>";
			#$URL="blob_image.php?id=$snapshot_id";
			#echo "<span class=' ' onclick=\"javascript: var newWin = window.open('$URL','SNAPSHOT','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=650, height=800, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>$snapshot_ref</span><BR>";
		};#while
		#die();
		foreach ($depthbed_array as $chrom=>$poss) {
			foreach ($poss as $pos=>$QCs) {
				$depthbed_stats_array[$chrom][$pos]["min"]=min($QCs);
				$depthbed_stats_array[$chrom][$pos]["max"]=max($QCs);
				$depthbed_stats_array[$chrom][$pos]["average"]=array_sum($QCs)/count($QCs);
				$depthbed_stats_array[$chrom][$pos]["stddev"]=standard_deviation($QCs);
				#$depthbed_stats_array[$chrom][$pos]["list"]=join(",",$QCs);
				$depthbed_stats_array[$chrom][$pos]["comment"]="Depth statistics:";
				$depthbed_stats_array[$chrom][$pos]["comment"].="&nbsp;&nbsp;&nbsp;Min=".$depthbed_stats_array[$chrom][$pos]["min"];
				$depthbed_stats_array[$chrom][$pos]["comment"].="&nbsp;&nbsp;&nbsp;Max=".$depthbed_stats_array[$chrom][$pos]["max"];
				$depthbed_stats_array[$chrom][$pos]["comment"].="&nbsp;&nbsp;&nbsp;Average=".$depthbed_stats_array[$chrom][$pos]["average"];
				$depthbed_stats_array[$chrom][$pos]["comment"].="&nbsp;&nbsp;&nbsp;Stadard Deviation=".$depthbed_stats_array[$chrom][$pos]["stddev"];

				foreach ($QCs as $QC_id=>$depth) {
					$depthbed_stats_array[$chrom][$pos]["QCrefs"].=$depthbed_QCidref_array[$QC_id].",";
				};#foreach
			};#foreach
		};#foreach

		#print "<PRE>"; print_r($depthbed_stats_array);  print "</PRE>";

	};#if


	$timestart=microtime(true);
	if (trim($sample_id) != "") {
		$query="SELECT sample.manifest_id, count(sample.id) AS nb_sample
			FROM sample
			WHERE sample.manifest_id IN
				(SELECT manifest.id
				FROM manifest
				INNER JOIN sample AS sample ON (manifest.id=sample.manifest_id)
				WHERE sample.id=$sample_id)
			GROUP BY sample.manifest_id
			";
		#print "<pre>$query</pre><BR>";

		$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));

		$manifest_id=0;
		while ($row = mysql_fetch_assoc($result)) {
			#echo "<pre>"; print_r($row); echo "</pre>";
			$manifest_id=$row["manifest_id"];
			$nb_sample_manifest=$row["nb_sample"]-1;
		};#while
	};#if
	$timestop=microtime(true);
	$timeDBexportManifestNbSample=$timestop-$timestart;
	#echo "Time DBexport Manifest NbSample: $timeDBexportManifestNbSample<BR>";


	$Variants_HTMLContent="";
	$nb_variant=0;

	# First read of the VCF
		if ($options["vcf_ids"]!="" && 1) {

			$variant_files=$options["vcf_ids"];

			# List of variants
			$query="SELECT variant.id
			FROM variant
			INNER JOIN variant_vcf ON (variant_vcf.variant_id=variant.id)
			WHERE variant_vcf.vcf_id IN ($variant_files)";
			#print "<pre>$query</pre><BR>";
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			$variant_ids_array[0]="0"; # for empty variant list
			while ($row = mysql_fetch_assoc($result)) {
				#print_r($row);
				$variant_ids_array[$row["id"]]=$row["id"];
			};#while
			$variant_ids=join(",",$variant_ids_array);
			#print $variant_ids;

			# Variant Causal Score

			# default variant score
			$query="SELECT flags.score AS flag_score, flags.label AS flag_label, flags.id AS flag_id, flags.definition AS flag_definition
				FROM flags
				WHERE type = 'variant'
				  AND score=0 # default
				";
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			while ($row = mysql_fetch_assoc($result)) {
				$variant_score_default=$row;
			};#foreach
			$query="SELECT min(flags.score) AS flag_score_min, max(flags.score) AS flag_score_max
				FROM flags
				WHERE type = 'variant'
				";
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			while ($row = mysql_fetch_assoc($result)) {
				$variant_score_min=$row["flag_score_min"];
				$variant_score_max=$row["flag_score_max"];
			};#foreach

			$variant_causal_score_table=array();
			$timestart=microtime(true);
			$query="SELECT variant.id AS variant_id, flags.score AS flag_score, flags.label AS flag_label, flags.id AS flag_id, flags.definition AS flag_definition
				FROM variant
				LEFT OUTER JOIN comment ON (comment.id=variant.id AND comment.type='variant')
				LEFT OUTER JOIN flags ON (flags.id=comment.flag_id AND flags.type=comment.type)
				WHERE variant.id IN ($variant_ids)
				GROUP BY variant_id
				ORDER BY abs(flags.score) ASC";
			#print "<pre>$query</pre><BR>"; # AND sample.manifest_id=$manifest_id
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			$timestop=microtime(true);
			$timeDBexportVariantSampleALL=$timestop-$timestart;
			while ($row = mysql_fetch_assoc($result)) {
				$variant_id=$row["variant_id"];
				#$variant_flag_score=$row["flag_score"];
				#$variant_flag_label=$row["flag_label"];
				#$variant_flag_id=$row["flag_id"];
				#$variant_flag_definition=$row["flag_definition"];
				$variant_causal_score_table[$variant_id]=$row;
			};#foreach
			#print "<pre>";print_r($variant_causal_score_table);print "</pre>";


			# Validation SCORE
			$timestart=microtime(true);

			$sample_where_clause=(trim("$sample_id")!="")?"  AND variant_sample.sample_id!=$sample_id ":"";
			$manifest_where_clause=(trim("$manifest_id")!="")?" AND sample.manifest_id=$manifest_id ":"";

			$query="SELECT variant_sample.variant_id AS variant_id, variant_sample.sample_id AS sample_id, flags.label, flags.id
				FROM variant_sample
				LEFT OUTER JOIN comment ON (comment.id=variant_sample.id AND comment.type='variant_sample')
				LEFT OUTER JOIN variant ON (variant.id=variant_sample.variant_id)
				LEFT OUTER JOIN flags ON (flags.id=comment.flag_id AND flags.type=comment.type)
				LEFT OUTER JOIN sample ON (sample.id=variant_sample.sample_id)
				WHERE variant_sample.variant_id IN ($variant_ids)
				  $sample_where_clause
				  $manifest_where_clause
				GROUP BY variant_id, sample_id
				ORDER BY abs(flags.score) ASC";

			#print "<pre>$query</pre><BR>"; # AND sample.manifest_id=$manifest_id
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			$timestop=microtime(true);
			$timeDBexportVariantSampleALL=$timestop-$timestart;
			#echo "Time DBexport VariantSample ALL: $timeDBexportVariantSampleALL<BR>";

		   	$summary= "<table class=' '>\n";
			$validation_table=array();
			$flag_0="inconnue";
			$nb=0;
			while ($row = mysql_fetch_assoc($result)) {
				$nb++;
				#print "$nb "; #print "<BR>";
				#print_r($row); print "<BR>";
				$variant_id=$row['variant_id'];
				$sample_id2=$row['sample_id'];
				$label=($row['label']=="")?"$flag_0":$row['label'];
				$flag_score=$row['score']+0;
				$flag_label[$flag]=$label;
				#$value=str_replace(",","<BR>",$row['value']);
				#$comment=$row['comment'];
				$validation_table_all[$variant_id][$flag_score][$sample_id2]=$row;
				#$validation_table_variant_id[$variant_id]=$value;
			};#while
			#print "<BR>total: $nb <BR>"; #print "<BR>";
			#die();
			#ksort($validation_table);
			#echo "<PRE>"; print_r($validation_table_all); echo "</PRE>";
			/*
			$variant_checked_summary_tmp="";
			$VDScore=0;
			$validation_summary="<TABLE><TR  valign='top' halign='left' style=''>
						<TD width=''>Validation</TD>
						<TD width='10px'>&nbsp;</TD>
						<TD width=''>Number of variants</TD>
						<TD width='10px'>&nbsp;</TD>
						<TD width=''>Score added</TD>
						</TR>";
			foreach ($validation_table as $flag=>$sample_ids) {
				$label=$flag_label[$flag];

				$variant_checked_summary_tmp.=$summary;
				$sign= $flag < 0 ? -1 : ( $flag > 0 ? 1 : 0 ) ;
				$VDScore_label=($sign*(count($sample_ids)/(pow(10,abs($flag)-1))));
				$VDScore+=$VDScore_label;
				$validation_summary.="<TR><TD>$label</TD><TD></TD><TD>".count($sample_ids)." samples</TD><TD></TD><TD>$VScore_label</TD></TR>";
				#echo "$VScore+=($sign*(count($sample_ids)/(pow(10,abs($flag)-1))));<BR>";
			};#foreach
			#$VScoreLog=-log(2,10);
			#$VScoreLog=(1/$VScore);
			$variant_array["VDSCORE"]="$VDScore";
			#echo "<BR><B>".$variant_array["VARIANT_ID"].".$sample_id: VDScore=$VDScore (".sprintf('%1.2E',$VDScore).")</B><TABLE valign='left' halign='top'>$variant_checked_summary_header$variant_checked_summary_tmp</TABLE><BR>";
			$validation_summary.="</TABLE>";
			*/

		};#if




	$header_main_array=array();
	$header_array=array();
	#$Variants_VCFContent_split=split("\n",$Variants_VCFContent);
	foreach (explode("\n",$Variants_VCFContent) as $line_nb=>$line) {
		#$Variants_HTMLContent.= "$line<BR>\n";
		$line=trim($line);
		if ($line=="") {
			continue;
		};#if
		#$line = "#qeydtsdfth";

		$pattern0 = '/^##(.*)=(.*)/';
		preg_match($pattern0, $line, $matches0, PREG_OFFSET_CAPTURE);

		$pattern = '/^##(.*)=<(.*)>/';
		preg_match($pattern, $line, $matches, PREG_OFFSET_CAPTURE);
		#print_r($matches); print "<BR>";

		$pattern2 = '/^#/';
		preg_match($pattern2, $line, $matches2, PREG_OFFSET_CAPTURE);



		#print "<pre>LINE"; print $line; print "</pre>";

		if (!empty($matches)) {

			# VCF header


			#print "<pre>MATCHES"; print_r($matches); print "</pre>";

			$type=$matches[1][0];
			$def=$matches[2][0];
			#print "<pre>DEF"; print_r($def); print "</pre>";

			$def_split=explode(",",$def);
			$def_array=array();
			$varID="";
			foreach ($def_split as $k=>$varval) {
				$varval_split=explode("=",$varval,2);
				#print_r($varval_split);
				$var=$varval_split[0];
				$val=$varval_split[1];
				if ($var=="ID") {
					$varID=$val;
				} else {
					$def_array[$var]=$val;

					# Special annotation for database
					if ($var=="Description") {
						$patternS = '/^.*\[(.*)\]/';
						preg_match($patternS, $line, $matchesS, PREG_OFFSET_CAPTURE);
						#print "<pre>def_array"; print_r($matchesS); print "</pre>";
						if (!empty($matchesS)) {
							$defS=$matchesS[1][0];
							$defS_split=explode(";",$defS);
							#print "<pre>def_array"; print_r($defS_split); print "</pre>";
							foreach ($defS_split as $k=>$varvalS ) {
								$varvalS_split=explode("=",$varvalS,2);
								$varS=$varvalS_split[0];
								$valS=$varvalS_split[1];
								#print "<pre>"; print("$varS=>$valS"); print "</pre>";
								$def_array[$varS]=$valS;
							};#foreach
						};#if
					};#if
				};#if
			};#foreach
			#print "<pre>def_array"; print_r($def_array); print "</pre>";
			$header_array[$type][$varID]=$def_array;

			/*
			$pattern3 = '/ID=(.*),/';
			preg_match($pattern3, $matches[2][0], $matches3, PREG_OFFSET_CAPTURE);
			if (!empty($matches3)) {
				print "<pre>"; print_r($matches3); print "</pre>";
			};#if
			*/


		} elseif (!empty($matches0)) {
			#print "<pre>MAIN INFOS DEF HEADER</pre>";
			#print "<pre>MATCHES0"; print_r($matches0); print "</pre>";
			$header_main_array[$matches0[1][0]]=$matches0[2][0];

		} elseif (!empty($matches2)) {

			#print "<pre>HEADER_ARRAY1"; print_r($header_main_array); print "</pre>";
			#print "<pre>HEADER_ARRAY2"; print_r($header_array); print "</pre>";

			# VCF Variants Header
			#$Variants_HTMLContent.="HEADER $line";
			$variant_header=explode("\t",$line);

			#$Variants_HTMLContent.= "$line<BR>\n";

			if (!empty($matches)) {
				# VCF Header

			} else {



			};#if

		} else {

			# Find Variant Information
			$variant_array=VCFLinetoArray($line,$variant_header);
			#print "<pre>"; print_r($variant_header); print "</pre>";
			#print "<pre>"; print_r($variant_array); print "</pre>";
			#$Variants_HTMLContent.="".join("\t",$variant_array)."<BR>";
			$nb_variant++;

			# variant in database
			if ($sample_id!="" && $variant_array["VARIANT_ID"]!="" && 1) {

				# variant flag
				$flag="0";
				$comment="";
				$variant_id=$variant_array["VARIANT_ID"];
				$timestart=microtime(true);
				$query="SELECT comment.id, comment.type, comment.flag_id, comment.comment, comment.history
					FROM comment
					WHERE comment.id='$variant_id.$sample_id'
					  AND comment.type='variant_sample'
					";
				#print "$query<BR>";
				$result = mysql_query($query,$mydbh) or die ("OUPS! ".mysql_error($mydbh)." <BR><PRE>$query</PRE>");
				$timestop=microtime(true);
				$timeDBexportVariantSampleComment=$timestop-$timestart;
				#echo "Time DBexport VariantSample Comment: $timeDBexportVariantSampleComment<BR>";
				$summary= "<table class=' '>\n";
				$comment_id="";
				$comment_type="";
				$flag="0";
				$rank="0";
				$comment="";
				$history="";
				$history_html="";
				while ($row = mysql_fetch_assoc($result)) {
					$comment_id=$row['id'];
					$comment_type=$row['type'];
					$comment_flag_id=$row['flag_id'];

					$comment=$row['comment'];
					$history=$row['history'];
					$history_html="<pre><small>$history</small></pre>";
					$summary.= "<tr>"
						."<td>".$row['flag_id']."</td>"
						."<td>".$row['comment']."</td>"
						."<td>".$row['history']."</td>"
					."</tr>";
					#echo "<PRE>$comment_id\n$comment_type\n$history</PRE><BR>";
				};#while
				#print "flag=$flag<BR>";
				$summary.= "</table>\n";
				#echo "<PRE>$query<BR>$summary</PRE>";

				$timestart=microtime(true);
				$query="SELECT id, type, score, rank, label, code, definition
					FROM flags
					WHERE type='variant_sample'
					ORDER BY rank ASC
					";
				#print "$query<BR>";
				$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
				$timestop=microtime(true);
				$timeDBexportVariantSampleFlag=$timestop-$timestart;
				#echo "Time DBexport VariantSample Flag: $timeDBexportVariantSampleFlag<BR>";
			   	$summary= "";
				$flags_options="";
				while ($row = mysql_fetch_assoc($result)) {
					$flag_score=$row['score'];
					$flag_rank=$row['rank'];
					$flag_id=$row['id'];
					#$flag_description=$row['label']." [".$row['level']."]";
					$flag_description=$row['label']; #." [".$row['level']."]";
					$summary.= "<tr>";
					$summary_head="<tr>";
					foreach ($row as $col=>$val) {
						$summary.= "<td>$val</td>";
						$summary_head.= "<td>$col</td>";
					};#foreach
					$summary.= "</tr>";
					$summary_head.= "</tr>";
					#$flags_array[$row['flag_value']]=$row['flag_label'];
					$flags_options.="<OPTION value='$flag_id' ".(($flag_id==$comment_flag_id)?"selected":"").">$flag_description [$flag_score]</OPTION>\n";
				};#while
				$summary= "<table class=' '>$summary_head $summary</table>\n";
				#echo $summary."<SELECT>$flags_options</SELECT>";

				# Validation SCORE
				/*
				$timestart=microtime(true);
				$query="SELECT variant_sample.variant_id AS variant_id, variant_sample.sample_id AS sample_id, variant_annotation.value, flags.label, flags.flag, comment.comment, sample.manifest_id
				FROM variant_sample
				LEFT OUTER JOIN comment ON (comment.id=variant_sample.id AND comment.type='variant_sample')
				LEFT OUTER JOIN variant ON (variant.id=variant_sample.variant_id)
				LEFT OUTER JOIN variant_annotation ON (variant_annotation.variant_id=variant.id)
				LEFT OUTER JOIN annotation ON (annotation.id=variant_annotation.annotation_id AND annotation.source='hgvs')
				LEFT OUTER JOIN flags ON (flags.flag=comment.flag AND flags.type=comment.type)
				LEFT OUTER JOIN sample ON (sample.id=variant_sample.sample_id)
				WHERE variant_sample.variant_id=$variant_id
				  AND variant_sample.sample_id!=$sample_id
				  AND sample.manifest_id=$manifest_id
				ORDER BY abs(comment.flag) ASC
				";

				#print "<pre>$query</pre><BR>";
				$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			   	$timestop=microtime(true);
				$timeDBexportVariantSample=$timestop-$timestart;
				#echo "Time DBexport VariantSample: $timeDBexportVariantSample<BR>";
				$summary= "<table class=' '>\n";
				$validation_table=array();
				$flag_0="inconnue";
				while ($row = mysql_fetch_assoc($result)) {
					#print_r($row);
					#print_r($row["manifest_id"]);
					$variant_id=$row['variant_id'];
					$sample_id2=$row['sample_id'];
					$label=($row['label']=="")?"$flag_0":$row['label'];
					$flag=$row['flag']+0;
					$flag_label[$flag]=$label;
					#$value=str_replace(",","<BR>",$row['value']);
					#$comment=$row['comment'];
					$validation_table[$flag][$sample_id2]=$row;
					#$validation_table_variant_id[$variant_id]=$value;
				};#while
				*/


				$validation_table=$validation_table_all[$variant_id];

				ksort($validation_table);

				#echo "$variant_id<PRE>"; print_r($validation_table); echo "</PRE>";

				$variant_checked_summary_tmp="";
				$VDScore=0;
				$VDInternalFrequency=0;
				$validation_summary="";
				$variant_array["VDSCORE"]="0";
				$variant_array["VDIF"]=0;
				#echo "<pre>"; print_r($validation_table); echo "</pre>";
				if ($validation_table!="") {

					$validation_summary_flag="";
					$nb_sample_ids=0;
					foreach ($validation_table as $flag=>$sample_ids) {
						$label=$flag_label[$flag];
						$nb_sample_ids+=count($sample_ids);
						$variant_checked_summary_tmp.=$summary;
						$sign= $flag < 0 ? -1 : ( $flag > 0 ? 1 : 0 ) ;
						$VDScore_label=($sign*(count($sample_ids)/(pow(10,abs($flag)-1))));
						$VDScore+=$VDScore_label;
						#$validation_summary_flag.="<TR><TD>$label</TD><TD></TD><TD>".count($sample_ids)." samples</TD><TD></TD><TD>$VDScore_label</TD></TR>";
						#$validation_summary_flag.="<TR><TD>".count($sample_ids)." </TD><TD>$label samples</TD><TD></TD><!--<TD></TD><TD>$VDScore_label</TD>--></TR>";
						$validation_summary_flag.=(($validation_summary_flag=="")?"":", ").count($sample_ids)." $label";
						#echo "$VScore+=($sign*(count($sample_ids)/(pow(10,abs($flag)-1))));<BR>";
					};#foreach
					#$VScoreLog=-log(2,10);
					#$VScoreLog=(1/$VScore);
					$variant_array["VDSCORE"]="$VDScore";
					$VDInternalFrequency=$nb_sample_ids/$nb_sample_manifest;
					$variant_array["VDIF"]=$VDInternalFrequency;
					#echo "<BR><B>".$variant_array["VARIANT_ID"].".$sample_id: VDScore=$VDScore (".sprintf('%1.2E',$VDScore).")</B><TABLE valign='left' halign='top'>$variant_checked_summary_header$variant_checked_summary_tmp</TABLE><BR>";
					$validation_summary="
						This variant was identified in $nb_sample_ids/$nb_sample_manifest other sample".(($nb_sample_manifest>1)?"s":"")."
						(Internal Frequency IF=".round($variant_array["VDIF"]*100,2)."%)<BR>

						$validation_summary_flag".(($validation_summary_flag=="")?"":"<BR>");
					if ($nb_sample_ids>0) {
					/*$validation_summary.="
						<!--<TABLE>
							<TR  valign='top' halign='left' style=''>
								<TD width=''>Validation</TD>
								<TD width='10px'>&nbsp;</TD>
								<TD width=''>Number of variants</TD>
								</TR>-->
							$validation_summary_flag2
						</TABLE>-->";*/
					};#if
				} else {

					$validation_summary="
						This variant was identified in 0/$nb_sample_manifest other sample".(($nb_sample_manifest>1)?"s":"")."
						(Internal Frequency IF=0%)<BR>
						";

				};#if
			};#if


			# SNAPSHOT
			if ($snapshot_array[$variant_array["CHROM"]][$variant_array["POS"]]!="") {
				$snapshot_table="<BR><TABLE><TR><TD>Snapshots:</TD><TD width='20'></TD>";
				#print "CHROM:".$variant_array["CHROM"]."<BR>";
				#print "POS:".$variant_array["POS"]."<BR>";
				foreach ($snapshot_array[$variant_array["CHROM"]][$variant_array["POS"]] as $aligner=>$snapshot_id) {
					#$snapshot_array[$snapshot_chrom][$snapshot_pos][$vcf_aligner]=$snapshot_id
					#print "<pre>"; print_r($row); print "</pre>";
					$URL="blob_image.php?id=$snapshot_id";
					$LINK="<span class=' ' onclick=\"javascript: var newWin = window.open('$URL','SNAPSHOT','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=650, height=800, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>$aligner</span>";
					$snapshot_table.="<TD>$LINK</TD><TD width='20'></TD>";

				};#foreach
				$snapshot_table.="</TR></TABLE>";
			};#if

			# Variant Output

			# Annotation
			$annotation_table="<TABLE>";
			foreach ($variant_array["INFO"] as $annotation_source=>$annotation_value) {
				#if (in_array($annotation_source,$annotation_list)) {
					$annotation_table.="<TR><TD valign='top' title=\"".str_replace("\"","",$header_array["INFO"][$annotation_source]["Description"])."\">$annotation_source</TD><TD width='20'></TD><TD valign='top'>".str_replace(",",", ",$annotation_value)."</TD></TR>\n";
				#};#if
			};#foreach
			$annotation_table.="</TABLE>";

			# Annotation classified
			$annotation_array=array();
			foreach ($variant_array["INFO"] as $annotation_source=>$annotation_value) {
				$annotationType=$header_array["INFO"][$annotation_source]["AnnotationType"];
				$annotationType=($annotationType=="")?"unknown":$annotationType;
				$annotation_array[$annotationType][$annotation_source]="<TD valign='top' title=\"".str_replace("\"","",$header_array["INFO"][$annotation_source]["Description"])."\">$annotation_source</TD><TD width='20'></TD><TD valign='top'>".str_replace(",",", ",$annotation_value)."</TD>";
				#if (in_array($annotation_source,$annotation_list)) {
				#$annotation_table.="<TR><TD valign='top' title=\"".str_replace("\"","",$header_array["INFO"][$annotation_source]["Description"])."\">$annotation_source</TD><TD width='20'></TD><TD valign='top'>".str_replace(",",", ",$annotation_value)."</TD></TR>\n";
				#};#if
			};#foreach
			$annotation_table_known="";
			$annotation_table_unknown="";
			foreach ($annotation_array as $annotationType=>$annotationType_array) {
				if ($annotationType!="unknown") {
					$annotation_table_known.="<TR><TD colspan='10'><B>".strtoupper($annotationType)."</B></TD></TR>\n";
					foreach ($annotationType_array as $annotation_source=>$annotation_source_content) {
						$annotation_table_known.="<TR>$annotation_source_content</TR>\n";
					};#foreach
				} else {
					$annotation_table_unknown.="<TR><TD colspan='10'><B>".strtoupper($annotationType)."</B></TD></TR>\n";
					foreach ($annotationType_array as $annotation_source=>$annotation_source_content) {
						$annotation_table_unknown.="<TR>$annotation_source_content</TR>\n";
					};#foreach
				};#if
			};#foreach

			#if () {

			#};#if
			$annotation_table="<TABLE>$annotation_table_known$annotation_table_unknown</TABLE>";
			#$annotation_table_unknown="<TABLE>$annotation_table_unknown</TABLE>";

			/*
			$annotation_table="<TABLE>
					<TR>
						<TD><B>Frequency</B></TD>
						<TD width='20px'></TD>
						<TD><B>Score</B></TD>
						<TD width='20px'></TD>
						<TD><B>Prediction</B></TD>
					</TR>
					<TR>
						<TD><TABLE><TR>".join("</TR><TR>",$annotation_array["frequency"])."</TR></TABLE></TD>
						<TD width='20px'></TD>
						<TD><TABLE><TR>".join("</TR><TR>",$annotation_array["score"])."</TR></TABLE></TD>
						<TD width='20px'></TD>
						<TD><TABLE><TR>".join("</TR><TR>",$annotation_array["prediction"])."</TR></TABLE></TD>
					</TR>

				</TABLE>";
			*/
			# Samples
			if (count($variant_array["SAMPLES"])>0) {
				$sample_table="";
				$sample_table_header="";
				$sample_table_lines="";
				# create header
				$sample_table_header_array=array();
				foreach ($variant_array["SAMPLES"] as $sample=>$qual_array) {
					foreach ($qual_array as $qual_name=>$qual_value) {
						$sample_table_header_array[$qual_name]=1;
					};#foreach
				};#foreach
				#print "<PRE>"; print_r($sample_table_header_array); print "</PRE>";

				foreach ($variant_array["SAMPLES"] as $sample=>$qual_array) {
					$sample_table_header_tmp="<TR><TD></TD>";
					$sample_table_lines.="<TR><TD>$sample</TD>";
					#foreach ($qual_array as $qual_name=>$qual_value) {
					foreach ($sample_table_header_array as $qual_name=>$qual_value_header) {
						$qual_value=$qual_array[$qual_name];
						$sample_table_header_tmp.="<TD width='20'></TD><TD>$qual_name</TD>";
						$sample_table_lines.="<TD width='20'></TD><TD>$qual_value</TD>";
					};#foreach
					$sample_table_header_tmp.="</TR>";
					if (count($qual_array)>1) {
						$sample_table_header=$sample_table_header_tmp;
					};#if
					$sample_table_lines.="</TR>";
				};#foreach
				$sample_table="<TABLE>$sample_table_header$sample_table_lines</TABLE>";
			};#if


			# Validation
			# <!--<option value=''>Flag the variant...</option>-->
			$validation="";
			$validation_width="300px";
			#print date(DATE_ATOM);
			#print date_default_timezone_set('UTC');

			$validation_form=0;
			if ($sample_id!="" && $variant_array["VARIANT_ID"]!="") {
				$validation_form=1;
			};#if

			if ($validation_form) {
				$variant_id=$variant_array["VARIANT_ID"];
				$validation="	<iframe id='ponyo_frame' name='iframe_flags_$nb_variant' border='0' style='visibility:hidden;float:right;display:none;'></iframe>
						<form id='form_flags_$nb_variant' method='post' action='validation.php' target='iframe_flags_$nb_variant' style='height:0px;visibility:;display:;'>
							<TABLE width='$validation_width' border=0 class='' style='padding:0px;margin-top:0px;'>
								<TR>

									<TD valign='top' width='99%'>
										<SELECT name='flag_id' style='margin-bottom:-5px;' width='99%'>$flags_options</SELECT>
									</TD>
									<TD valign='top'>
										<input type='submit' id='validation' name='validation' value='Save' class='btn button search' onclick=\"javascript:;\">
										<input data-original-title='' value='sample_id' type='hidden'>
										<input type='hidden' id='type' name='type' style='visibility:hidden;float:right;display:none;' value='variant_sample'></input>
										<input type='hidden' id='project_id' name='project_id' style='visibility:hidden;float:right;display:none;' value='$project_id'></input>
										<input type='hidden' id='run' name='run' style='visibility:hidden;float:right;display:none;' value='$input_run'></input>
										<input type='hidden' id='id' name='id' style='visibility:hidden;float:right;display:none;' value='".$variant_array["VARIANT_ID"].".".$sample_id."'></input>
										<input type='hidden' id='author' name='author' style='visibility:hidden;float:right;display:none;' value='".USERNAME."'></input>
										<input data-original-title='' onclick=\"document.getElementById('form_flags_$nb_variant').submit();\" type='hidden'>
									</TD>
									<TD valign='top'>
										<input type='button' id='history' name='history' value='H' class='btn button search' onclick=\"javascript:;document.getElementById('variant_validation_history_$nb_variant').src='validation_history.php?id=$variant_id.$sample_id&type=variant_sample&project_id=$project_id&rand='+Math.floor((Math.random()*10000000)+1);show_hide('variant_validation_history_$nb_variant');\"></TD>

								</TR>
								<TR>

									<TD colspan=3 align='center' valign='top'>
										<TEXTAREA id='comment' name='comment' class='' style='border: none; width: 99%; height: 100%; -webkit-box-sizing: border-box; -moz-box-sizing: border-box; box-sizing: border-box; padding:0px;' rows=3>$comment</TEXTAREA>
									</TD>

								</TR>
								<TR>
									<TD colspan=3 align='left' valign='right' zalign='top' >
										<iframe id='variant_validation_history_$nb_variant' name='variant_validation_history_$nb_variant' border='1' style='float:left;text-align:right;display:none;z-index:2;position:absolute;' marginwidth='0' marginheight='0' frameborder='0' vspace='0' hspace='0' width='$validation_width'>$history_html</iframe>
									</TD>

								</TR>
							</TABLE>
						</form>

				";
			};#if

			# Output
			#<pre>".join("\t",$variant_array)."<BR>".implode_with_key($variant_array["INFO"],"\t","<BR>")."<BR>".implode_with_key($variant_array["SAMPLES"],"\t","<BR>")."<BR><BR><BR></pre> ".implode_with_key($annotation_list,"\t","<BR>")."
			#$Variants_HTMLContent.="

			# Ensembl Link
			$ensembl_prefix="grch37"; # www by default
			$ensembl_links="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Location/View?db=core;r=".$variant_array["REGION"]."' target='ensembl'>Location</a></TD><TD width='20'></TD>";
			if ($variant_array["ENSEMBL"]!="") {
				$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=".$variant_array["ENSEMBL"]."' target='ensembl'>Gene</a></TD><TD width='20'></TD>";
			};#if
			if ($variant_array["ENSEMBL"]!="" && $variant_array["REGION"]!=""
				&& ( $variant_array["DBSNP"]!="" || $variant_array["COSMIC"]!="" || $variant_array["HGMD"]!="" )) {
				#print "$v_id=(".$variant_array["DBSNP"]." ".$variant_array["HGMD"]." ".$variant_array["COSMIC"]."";
				if ($variant_array["DBSNP"]!="") {
					$v_id=$variant_array["DBSNP"];
				} elseif ($variant_array["HGMD"]!="") {
					$v_id=$variant_array["HGMD"];
				} elseif ($variant_array["COSMIC"]!="") {
					$v_id=$variant_array["COSMIC"];
				};#if
				#$v_id=($variant_array["DBSNP"]!="")?$variant_array["DBSNP"]:($variant_array["HGMD"]!="")?$variant_array["HGMD"]:($variant_array["COSMIC"]!="")?$variant_array["COSMIC"]:"";
				#print "v_id: $v_id<BR>";
				$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Variation/Mappings?db=core;gene=".$variant_array["ENSEMBL"].";v=$v_id;r=".$variant_array["REGION"].";vdb=variation' target='ensembl'>Variation</a></TD><TD width='20'></TD>";
				$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Variation/Context?db=core;gene=".$variant_array["ENSEMBL"].";v=$v_id;r=".$variant_array["REGION"].";vdb=variation' target='ensembl'>Context</a></TD><TD width='20'></TD>";
				$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Variation/Phenotype?db=core;gene=".$variant_array["ENSEMBL"].";v=$v_id;r=".$variant_array["REGION"].";vdb=variation' target='ensembl'>Phenotype</a></TD><TD width='20'></TD>";
				if ($variant_array["DBSNP"]!="") {
					$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Variation/Population?db=core;gene=".$variant_array["ENSEMBL"].";v=$v_id;r=".$variant_array["REGION"].";vdb=variation' target='ensembl'>Population</a></TD><TD width='20'></TD>";
				};#if

			};#if
			if ($variant_array["ENSEMBL"]!="") {
				$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Gene/Variation_Gene/Table?db=core;gene=".$variant_array["ENSEMBL"].";v=$v_id;r=".$variant_array["REGION"].";vdb=variation' target='ensembl'>Variation_Gene</a></TD><TD width='20'></TD>";
			};#if

			# Global filter
			$global_filter=0;
			if ($options["global_filter"]=="") {
				$global_filter=1;
			} else {
				foreach ($variant_array as $info=>$info_value) {
					if (!is_array($info_value)) {
						#print "$info=>$info_value $global_filter ";
						$global_filter=($global_filter || (strpos($info_value,$options["global_filter"])!== false));
						#print " => $global_filter <BR>";
					};#if
				};#foreach
			};#if
			#print "global_filter=$global_filter<BR>";

			if (	(!$options["hardfiltering"] || $variant_array["FILTER_FLAG"]=="PASS")
				&& $global_filter
				&& ($options["gene_filter"]==""
					|| strpos($variant_array["SYMBOL"],$options["gene_filter"])!== false
					|| strpos($variant_array["HGVS"],$options["gene_filter"])!== false
					)
			) {

			$variant_nb++;
			$variant_list["SYMBOL"][$variant_array["SYMBOL"]]++;
			$variant_list["FILTER_FLAG"][$variant_array["FILTER_FLAG"]]++;

			if ($options["orderby"]=="") {
				$orderby="FILTER_SCORE";
			} else {
				$orderby=strtoupper($options["orderby"]);
			};#if

			#echo "<pre>"; print_r($variant_list["SYMBOL"]); echo "</pre>";
			$variant_id=$variant_array["VARIANT_ID"];
			#.sprintf('%1.2E',$VScore).
			#$order_by="SYMBOL";
			#echo "$orderby ".$variant_array[$orderby]." ".($variant_array[$orderby]+1)."<BR>";
			#print_r($variant_array); echo "<BR>";
			#$Variants_HTMLContent[$variant_array["FILTER_SCORE"]].="
			#".$variant_array["DBSNP"]."



			# GROUP DATABASE
			#echo "<PRE>";
			$group_database_annotation="";
			foreach ($user_groups as $group_id=>$group_ref) {
				#echo "# GROUP ID: $group_id\n";
				#echo "# GROUP REF: $group_ref\n";
				#print_r($variant_array["INFO"]);
				if (trim($variant_array["INFO"]["database_$group_ref"])!="") {
					#echo "# ".$variant_array["INFO"]["database_$group_ref"]."\n";
					$sep=(($group_database_annotation)=="")?"":". ";
					$group_database_annotation.=$sep.trim($variant_array["INFO"]["database_$group_ref"]);
				};#if
			};#foreach
			#echo "</PRE>";
			#&nbsp;&nbsp;&nbsp;<span title='Variant Depth'>DP=".$variant_array["DP"]."</span>
			#&nbsp;&nbsp;&nbsp;<span title='Variant Depth'>DP=".$depthbed_stats_array[str_replace("chr","",$variant_array["CHROM"])][$variant_array["POS"]]["min"]."</span>




			$customHGVS=find_hgvs($variant_array["HGVS"],$variant_array["SYMBOL"],$user_groups,$customNM_infos);
			#print "<PRE>!$customHGVS</PRE>";
			$customHGVS=str_replace(",","<BR>",$customHGVS);

			$variant_color="gray";
			$variant_validation="";
			#$variant_score_max=5;
			#$variant_score_min=0;
			$variant_score=($variant_causal_score_table[$variant_id]["flag_score"]!="")?$variant_causal_score_table[$variant_id]["flag_score"]:"?"; #print_r($variant_causal_score_table[$variant_id]);
			$variant_score_definition=($variant_causal_score_table[$variant_id]["flag_definition"]!="")?$variant_causal_score_table[$variant_id]["flag_definition"]:$variant_score_default["flag_definition"]; #print_r($variant_causal_score_table[$variant_id]);
			#print_r($variant_score);
			$variant_color_scale_max=512;
			$variant_color_scale_min=1;
			$variant_color_code=($variant_score*$variant_color_scale_max)/$variant_score_max;
			$variant_color_code_inv=round($variant_color_scale_max-$variant_color_code);
			$variant_color_codeR=(($variant_color_code_inv-255)>255)?255:$variant_color_code_inv-255;
			$variant_color_codeG=(($variant_color_code_inv-255)>255)?255:$variant_color_code_inv-255;
			$variant_color_codeB=($variant_color_code_inv>255)?255:$variant_color_code_inv;
			$variant_color="rgb($variant_color_codeR,$variant_color_codeG,$variant_color_codeB)";
			$variant_validation="<DIV  onclick=\"javascript:show_hide('variant_content_$nb_variant');show_hide('variant_links_$nb_variant');\"
										style='halign=middle;float:right;cursor:pointer;background-color:$variant_color;height:15px;width:15px;'
										title='Variant Causative Score'></DIV>";
			$variant_validation="";



			$variant_frequency_tag="Unknown";
			$variant_rare_frequency=0.05;
			$variant_error_frequency=0.1;
			$COV_IF_EF=(abs($variant_array["POPFREQ"]-$variant_array["VDIF"])/abs($variant_array["POPFREQ"]+$variant_array["VDIF"]))+0;
			if ($variant_array["DBSNP_NONFLAGGED"] != "") {
				$variant_frequency_tag="Polymorphism";
				if (	($variant_array["VDIF"] > $variant_error_frequency && $variant_array["POPFREQ"] < $variant_error_frequency)
					|| ($variant_array["VDIF"] > $variant_error_frequency && $variant_array["POPFREQ"] > $variant_error_frequency && $COV_IF_EF>0.5)
				) {
					$variant_frequency_tag.="/Error?".$variant_array["EF"];
				#} else {
				#	$variant_frequency_tag="Polymorphism";
				};#if
			} elseif ($variant_array["DBSNP"] != "" && $variant_array["VDIF"] <= $variant_error_frequency) {
				$variant_frequency_tag="Rare/Pathogenic";
			} elseif ($variant_array["VDIF"] > $variant_error_frequency) {
				$variant_frequency_tag="Error?";
			};#if
			#$variant_frequency_tag.=$COV_IF_EF;
			$variant_frequency_tag_title="Variant type: Based on dbSNP annotation and Internal and External frequencies (Briefly, polymorphism of annotated with dbSNP NonFlagged, Rare/pathogenic is annotated in dbSNP but not as a polymorphism, Error if IF and EF inconsistante. Rare max frequency = $variant_rare_frequency, Error min frequency = $variant_error_frequency)";

			$Variants_HTMLContent[$variant_array[$orderby]][$variant_nb]="

				<SCRIPT>
				function show_hide(id) {
					if(document.getElementById(id).style.display=='none') {
						document.getElementById(id).style.display='inline';
					} else {
						document.getElementById(id).style.display='none';
					}
					return true;
				}
				</SCRIPT>
				<table border=0 width='100%' class='' ><TR><TD>


					<table border=0 width='100%' class='variant_mini_header' >
						<TR>
							<TD class=' variant_header' valign='top' halign='left' height='' width='350'>
								<div  style='float:;' height='100%'>

									<DIV  onclick=\"javascript:show_hide('variant_content_$nb_variant');show_hide('variant_links_$nb_variant');\"
										style='float:right;cursor:pointer;'
										title='show or hide variant content'>[+]</DIV>

									<b>".(($variant_array["SYMBOL"]=="")?"unknown":$variant_array["SYMBOL"]."<BR>")."</b><br>
									$variant_validation
									<!--[#id:$variant_id]--> ".$variant_array["CHROMPOS"]."<br>
									".(($customHGVS=="")?"":"$customHGVS<BR>")."
									<DIV>
									<!--".($variant_array["DBSNP"]==""?"":" &nbsp;&nbsp;&nbsp;<span title='Polymorphism'>".(($variant_array["DBSNP"] != "" && $variant_array["DBSNP_NONFLAGGED"] == "") ? " [Rare/Pathogenic]" : "[SNP]")." </span>")."-->
									</DIV>


								</div>
							</TD>
							<TD valign='top' align='left'  class=' variant_mini_header' width=''>
								<div class=''style='float:middle;' height='100%'>
								<b>".$variant_array["FILTER_FLAG"]."</b>
								&nbsp;&nbsp;&nbsp;<span title='$variant_frequency_tag_title'>[$variant_frequency_tag]</span>
								&nbsp;&nbsp;&nbsp;<span title='PrioritiZation Score'>PZ=".$variant_array["FILTER_SCORE"]."</span>
								&nbsp;&nbsp;&nbsp;<span title='Variant Depth'>DP=".$variant_array["DP"]."</span>
								&nbsp;&nbsp;&nbsp;<span title='Genotype Quality'>GQ=".$variant_array["GQ"]." </span>
								&nbsp;&nbsp;&nbsp;<span title='Allele Frequency'>VAF=".round($variant_array["AlleleFrequency"]*100,2)."%</span>
								&nbsp;&nbsp;&nbsp;<span title='Variant Internal Frequency'>IF=".round($variant_array["VDIF"]*100,2)."%</span>
								&nbsp;&nbsp;&nbsp;<span title='Variant External Frequency'>EF=".round($variant_array["POPFREQ"]*100,2)."%</span>
								&nbsp;&nbsp;&nbsp;<span title='Variant Causative Score: $variant_score_definition [$variant_score]'>CS=$variant_score</span>

								<!--&nbsp;&nbsp;&nbsp;<span title='ValiDation Score Score'>VD=".((abs($variant_array["VDSCORE"])>0.001||$variant_array["VDSCORE"]==0)?$variant_array["VDSCORE"]:sprintf('%1.2E',$variant_array["VDSCORE"]))."</span>-->
								<!--&nbsp;&nbsp;&nbsp;<span title=''>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;</span>-->



								".($group_database_annotation==""?"":((strlen($group_database_annotation)>4)?"<BR><span title='Group database annotation'>$group_database_annotation</span>":"&nbsp;&nbsp;&nbsp;<span title='Group database annotation'>[$group_database_annotation]</span>"))."
								<BR>
								".$variant_array["FILTER_COMMENT"]."
								<BR>
								</div>

							</TD>".((!$validation_form)?"":"
							<TD valign='top' align='right' width='$validation_width' class=' variant_validation' style='backgroud-color:white;'>
								<div style='float:right;text-align:right;display:;margin-left:0px;margin-top:;' id='variant_validation_$nb_variant'>
									$validation
								</div>
							</TD>")."
						</TR>

					</TABLE>
				</TD></TR>
				</TABLE>
				<table border=0 width='100%' class=' ' >
				<TR id='variant_content_$nb_variant' style='display:none;'>
				<!--<TD width='20px'></TD>-->
				<TD width='99%' class=''  >
				<TABLE border=0 width='100%' height='100%' class='variant_infos'><TR>
					<TD width='20px'></TD>
					<TD >
					$snapshot_table
					<BR>
					$sample_table
					<!--<A HREF='".str_replace(">","|","search.php?s=Search&q=".$variant_array["CHROMPOS"])."'>Variant in the database</A><BR>-->
					<!--<BR>$validation_summary-->
					<!--<BR>-->
					<!--".$depthbed_stats_array[str_replace("chr","",$variant_array["CHROM"])][$variant_array["POS"]]["comment"]."<BR>
					<BR>-->
					<TABLE>
						<TR>
							<TD>Ensembl</TD>
							<TD width='20'></TD>
							$ensembl_links

						</TR>
						<TR>
							<TD>Links</TD>
							<TD width='20'></TD>
							<TD><a href='http://exac.broadinstitute.org/variant/".$variant_array["CHROM"]."-".$variant_array["POS"]."-".$variant_array["REF"]."-".$variant_array["ALT"]."' target='exac'>ExAC</a></TD><TD width='20'></TD>
							<TD><a href='http://www.ncbi.nlm.nih.gov/clinvar?term=(".str_replace("chr","",$variant_array["CHROM"])."[Chromosome]) AND ".$variant_array["POS"]."[Base Position]' target='clinvar'>ClinVar</a> </TD><TD width='20'></TD>
							".(($variant_array["DBSNP"]!="")?"<TD><a href='http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=".$variant_array["DBSNP"]."' target='dbsnp'>dbSNP</a></TD><TD width='20'></TD>":"")."
							".(($variant_array["DBSNP"]!="")?"<TD><a href='https://www.snpedia.com/index.php/".$variant_array["DBSNP"]."' target='SNPedia'>SNPedia</a></TD><TD width='20'></TD>":"")."
							".(($variant_array["DBSNP"]!="")?"<TD><a href='http://www.ncbi.nlm.nih.gov/pubmed?term=".$variant_array["SYMBOL"]." AND ".$variant_array["PNOMEN"]."' target='PubMed'>PubMed</a></TD><TD width='20'></TD>":"")."
							".(($variant_array["DBSNP"]!="")?"<TD><a href='http://www.ncbi.nlm.nih.gov/variation/view/?filters=source:dbsnp&q=".$variant_array["DBSNP"]."' target='VariationViewer'>Variation</a></TD><TD width='20'></TD>":"")."
							".(($variant_array["GNOMEN"]!="" && $variant_array["PNOMEN"]!="")?"<TD><a href='https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/LitVar/#!?query=".$variant_array["GNOMEN"]."%20".$variant_array["PNOMEN"]."' target='LitVar'>LitVar</a></TD><TD width='20'></TD>":"")."
							".(($variant_array["GNOMEN"]!="" && $variant_array["PNOMEN"]!="")?"<TD><a href='https://cancer.sanger.ac.uk/cosmic/search?q=".$variant_array["GNOMEN"]."+".$variant_array["PNOMEN"]."' target='COSMIC'>COSMIC</a></TD><TD width='20'></TD>":"")."
							".(($variant_array["GNOMEN"]!="" && $variant_array["PNOMEN"]!="")?"<TD><a href='https://varsome.com/variant/hg19/".$variant_array["GNOMEN"]."%20".$variant_array["PNOMEN"]."' target='VarSome'>VarSome</a></TD><TD width='20'></TD>":"")."
							<TD width='20'></TD>
						</TR>
						<TR><TD></TD>
							<TD width='20'></TD>
							".(($variant_array["GNOMEN"]!="")?"<TD><a href='https://www.ncbi.nlm.nih.gov/gap/advanced_search/?TERM=".$variant_array["GNOMEN"]."' target='dbGaP'>dbGaP</a></TD><TD width='20'></TD>":"")."
							".(($variant_array["SYMBOL"]!="")?"<TD><a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=".$variant_array["SYMBOL"]."' target='genecards'>GeneCards</a></TD>
							<TD width='20'></TD>":"")."
							".(($variant_array["CHROM"]!="" && $variant_array["POS"]!="" && $variant_array["ALT"]!="" && $variant_array["REF"]!="")?"<TD><a href='http://www.fudan-pgx.org/premedkb/index.html#/search/result?queryType=3&step=1&num=1&term=%27".$variant_array["CHROMNUM"]."-".$variant_array["POS"]."-".$variant_array["REF"]."-".$variant_array["ALT"]."%27%5Bvariant%5D' target='_blank'>PreMedKB</a></TD>
							<TD width='20'></TD>":"")."
							".(($variant_array["DBSNP"]!="")?"<TD><a href='https://www.pharmgkb.org/rsid/".$variant_array["DBSNP"]."' target='PharmGKB'>PharmGKB</a></TD><TD width='20'></TD>":"")."
							<TD width='20'></TD>

						</TR>
					</TABLE>
					<BR>
					$annotation_table

					</TD></TR></TABLE>

				</TD></TR></TABLE>
			";
			};#if


		};#if
	};#foreach

	# Summary
	$summary="<!--<div class=' variant_header' style='clear:left;width:100%'>-->
			<table border=0 width='100%' class='  variant_header block' >
			<TR>
			<TD width='20px'></TD>
			<TD>
			<BR>
			".(count($variant_list["SYMBOL"])+0)." gene".(((count($variant_list["SYMBOL"])+0)>1)?"s":"")."&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
			".($variant_nb+0)." variant".((($variant_nb+0)>1)?"s":"")."&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
			".($variant_list["FILTER_FLAG"]["PASS"]+0)." PASS &nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp; ".(array_sum($variant_list["FILTER_FLAG"])-$variant_list["FILTER_FLAG"]["PASS"]+0)." FILTERED
			<BR><BR></TD></TR>
			</TABLE>

		";
		# $variant_list["FILTER_FLAG"]["FILTERED"]+$variant_list["FILTER_FLAG"]["FILTERED,PASS"]
	# ASC/DESC sorting
	if ($options["ascdesc"]=="DESC") {
		krsort($Variants_HTMLContent);
	} elseif ($options["ascdesc"]=="ASC") {
		ksort($Variants_HTMLContent);
	};#if

	$v_nb=0;
	$return=$summary."";
	$limited=0;
	foreach ($Variants_HTMLContent as $korder=>$Variants_list) {
		foreach ($Variants_list as $Variant_id=>$Variant_profile) {
			$v_nb++;
			if ($v_nb<=$options["limit"] || $options["limit"]=="") {
				$return.=$Variant_profile;
			} elseif (!$limited) {
				$return.="
						<table border=0 width=100%' class='  variant_header block' >
						<TR>
						<TD width='20px'></TD>
						<TD>

						... List of variants limited to ".$options["limit"]." (out of ".($variant_nb+0).")
						</TD></TR>
						</TABLE>


					";
				$limited=1;
			};#if
		};#foreach
	};#foreach

	#$return=$summary.implode("",$Variants_HTMLContent);

	return $return;

	#return $Variants_HTMLContent;

}#function VCFtoHTML


function VCFFiletoHTML($Variants_VCFContent,$sample_id=0,$annotation_list=array("ALL"),$options=array()) {

	#require "connect.php";
	require "config.php";

	if ($project_id=="") {
		$project_id=1;
	};#if

	if (0) {
	$handle = fopen($Variants_VCFContent, "r");
	if ($handle) {
	    while (($line = fgets($handle)) !== false) {
		// process the line read.
		#print "<plaintext>$line</plaintext>";
		$line_split=explode("\t",$line);
		if ($line_split[0][0] != "#") {
			#echo $line_split[0].":".$line_split[1]."<BR>";
			$all_positions[$line_split[0]][$line_split[1]]=1;
		}#if
	    }
	    fclose($handle);
	} else {
	    // error opening the file.
		print "error opening file";
	}
	#print "<pre>"; print_r($all_positions); print "</pre>";
	#return 0;
	};#if

	# found positions
	if (0) {
	$all_positions=array();
	foreach (explode("\n",$Variants_VCFContent) as $line_nb=>$line) {
		$line_split=explode("\t",$line);
		if ($line_split[0][0] != "#") {
			#echo $line_split[0].":".$line_split[1]."<BR>";
			$all_positions[$line_split[0]][$line_split[1]]=1;
		}#if
	};#foreach
	};#fi





	# SNAPSHOT
	if ($options["vcf_ids"]!="" || $options["report"]!="") {
		$variant_files=$options["vcf_ids"];
		#print_r($options["vcf_ids"]);
		if ($options["vcf_ids"]!="") {
			$where_clause="vcf.id IN ($variant_files)";
		} elseif ($options["report"]!="") {
			$where_clause="vcf.sample_id=$sample_id";
		};#if

		# SNAPSHOT
		$query="SELECT snapshot.id AS snapshot_id, snapshot.ref AS snapshot_ref, snapshot.chrom AS snapshot_chrom, snapshot.pos AS snapshot_pos, vcf.aligner AS vcf_aligner, vcf.bam, vcf.ref
			FROM vcf
			INNER JOIN snapshot_vcf ON (snapshot_vcf.vcf_id=vcf.id)
			INNER JOIN snapshot ON (snapshot.id=snapshot_vcf.snapshot_id)
			WHERE $where_clause
			";
		#print "$query<BR>";

		$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));

		while ($row = mysql_fetch_assoc($result)) {
			$snapshot_id=$row["snapshot_id"];
			$snapshot_ref=$row["snapshot_ref"];
			$snapshot_chrom=$row["snapshot_chrom"];
			$snapshot_pos=$row["snapshot_pos"];
			$vcf_aligner=$row["vcf_aligner"];
			$snapshot_array[$snapshot_chrom][$snapshot_pos][$vcf_aligner]=$snapshot_id;
			#print "<pre>"; print_r($row); print "</pre>";
			#$URL="blob_image.php?id=$snapshot_id";
			#echo "<span class=' ' onclick=\"javascript: var newWin = window.open('$URL','SNAPSHOT','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=650, height=800, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>$snapshot_ref</span><BR>";
		};#while

	};#if


	# QC
	if (0 && $options["vcf_ids"]!="" || $options["report"]!="") {
		#echo "QC_depthbed";
		$variant_files=$options["vcf_ids"];
		#print_r($options["vcf_ids"]);
		if ($options["vcf_ids"]!="") {
			$where_clause="vcf.id IN ($variant_files)";
		} elseif ($options["report"]!="") {
			$where_clause="vcf.sample_id=$sample_id";
		};#if

		# QC
		$query="SELECT distinct QC.id, QC.ref, QC.QC
			FROM vcf
			INNER JOIN QC_association ON (QC_association.object_id=vcf.id AND QC_association.object_type='vcf')
			INNER JOIN QC ON (QC.id=QC_association.QC_id)
			WHERE QC.metrics='depthbed'
			  AND QC.type='AlignedBAM'
			  AND $where_clause";
		#print "$query<BR>";

		$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		while ($row = mysql_fetch_assoc($result)) {
			#print "<PRE>"; print_r($row);  print "</PRE>";
			$QC_id=$row["id"];
			$QC_ref=$row["ref"];
			$QC_QC=$row["QC"];
			$QC_QC_split=explode("\n",$QC_QC);
			#print "<PRE>"; print_r($row);  print "</PRE>";

			foreach ($QC_QC_split as $line_nb=>$QC_QC_line) {
				#print "$line_nb=>$QC_QC_line<BR>\n";
				$QC_QC_line_split=explode("\t",$QC_QC_line);
				$QC_QC_line_chrom=str_replace("chr","",$QC_QC_line_split[0]);
				$QC_QC_line_pos=$QC_QC_line_split[1];
				if ($all_positions[$QC_QC_line_split[0]][$QC_QC_line_pos] > 0 ) {
					print "# Depth position ".$QC_QC_line_split[0].$QC_QC_line_pos."<BR>";
					$QC_QC_line_depth=$QC_QC_line_split[2];
					$depthbed_array[$QC_QC_line_chrom][$QC_QC_line_pos][$QC_id]=$QC_QC_line_depth;
					$depthbed_QCidref_array[$QC_id]=$QC_ref;
				};#if
			};#foreach
			#print "<pre>"; print_r($row); print "</pre>";
			#$URL="blob_image.php?id=$snapshot_id";
			#echo "<span class=' ' onclick=\"javascript: var newWin = window.open('$URL','SNAPSHOT','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=650, height=800, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>$snapshot_ref</span><BR>";
		};#while
		#die();
		foreach ($depthbed_array as $chrom=>$poss) {
			foreach ($poss as $pos=>$QCs) {
				$depthbed_stats_array[$chrom][$pos]["min"]=min($QCs);
				$depthbed_stats_array[$chrom][$pos]["max"]=max($QCs);
				$depthbed_stats_array[$chrom][$pos]["average"]=array_sum($QCs)/count($QCs);
				$depthbed_stats_array[$chrom][$pos]["stddev"]=standard_deviation($QCs);
				#$depthbed_stats_array[$chrom][$pos]["list"]=join(",",$QCs);
				$depthbed_stats_array[$chrom][$pos]["comment"]="Depth statistics:";
				$depthbed_stats_array[$chrom][$pos]["comment"].="&nbsp;&nbsp;&nbsp;Min=".$depthbed_stats_array[$chrom][$pos]["min"];
				$depthbed_stats_array[$chrom][$pos]["comment"].="&nbsp;&nbsp;&nbsp;Max=".$depthbed_stats_array[$chrom][$pos]["max"];
				$depthbed_stats_array[$chrom][$pos]["comment"].="&nbsp;&nbsp;&nbsp;Average=".$depthbed_stats_array[$chrom][$pos]["average"];
				$depthbed_stats_array[$chrom][$pos]["comment"].="&nbsp;&nbsp;&nbsp;Stadard Deviation=".$depthbed_stats_array[$chrom][$pos]["stddev"];

				foreach ($QCs as $QC_id=>$depth) {
					$depthbed_stats_array[$chrom][$pos]["QCrefs"].=$depthbed_QCidref_array[$QC_id].",";
				};#foreach
			};#foreach
		};#foreach

		#print "<PRE>"; print_r($depthbed_stats_array);  print "</PRE>";

	};#if


	$timestart=microtime(true);
	if (trim($sample_id) != "") {
		$query="SELECT sample.manifest_id, count(sample.id) AS nb_sample
			FROM sample
			WHERE sample.manifest_id IN
				(SELECT manifest.id
				FROM manifest
				INNER JOIN sample AS sample ON (manifest.id=sample.manifest_id)
				WHERE sample.id=$sample_id)
			GROUP BY sample.manifest_id
			";
		#print "<pre>$query</pre><BR>";

		$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));

		$manifest_id=0;
		while ($row = mysql_fetch_assoc($result)) {
			#echo "<pre>"; print_r($row); echo "</pre>";
			$manifest_id=$row["manifest_id"];
			$nb_sample_manifest=$row["nb_sample"]-1;
		};#while
	};#if
	$timestop=microtime(true);
	$timeDBexportManifestNbSample=$timestop-$timestart;
	#echo "Time DBexport Manifest NbSample: $timeDBexportManifestNbSample<BR>";


	$Variants_HTMLContent=array();
	$nb_variant=0;

	# First read of the VCF
		if ($options["vcf_ids"]!="" && 1) {

			$variant_files=$options["vcf_ids"];

			# List of variants
			$query="SELECT variant.id
			FROM variant
			INNER JOIN variant_vcf ON (variant_vcf.variant_id=variant.id)
			WHERE variant_vcf.vcf_id IN ($variant_files)";
			#print "<pre>$query</pre><BR>";
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			$variant_ids_array[0]="0"; # for empty variant list
			while ($row = mysql_fetch_assoc($result)) {
				#print_r($row);
				$variant_ids_array[$row["id"]]=$row["id"];
			};#while
			$variant_ids=join(",",$variant_ids_array);
			#print $variant_ids;

			# Variant Causal Score

			# default variant score
			$query="SELECT flags.score AS flag_score, flags.label AS flag_label, flags.id AS flag_id, flags.definition AS flag_definition
				FROM flags
				WHERE type = 'variant'
				  AND score=0 # default
				";
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			while ($row = mysql_fetch_assoc($result)) {
				$variant_score_default=$row;
			};#foreach
			$query="SELECT min(flags.score) AS flag_score_min, max(flags.score) AS flag_score_max
				FROM flags
				WHERE type = 'variant'
				";
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			while ($row = mysql_fetch_assoc($result)) {
				$variant_score_min=$row["flag_score_min"];
				$variant_score_max=$row["flag_score_max"];
			};#foreach

			$variant_causal_score_table=array();
			$timestart=microtime(true);
			$query="SELECT variant.id AS variant_id, flags.score AS flag_score, flags.label AS flag_label, flags.id AS flag_id, flags.definition AS flag_definition
				FROM variant
				LEFT OUTER JOIN comment ON (comment.id=variant.id AND comment.type='variant')
				LEFT OUTER JOIN flags ON (flags.id=comment.flag_id AND flags.type=comment.type)
				WHERE variant.id IN ($variant_ids)
				GROUP BY variant_id
				ORDER BY abs(flags.score) ASC";
			#print "<pre>$query</pre><BR>"; # AND sample.manifest_id=$manifest_id
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			$timestop=microtime(true);
			$timeDBexportVariantSampleALL=$timestop-$timestart;
			while ($row = mysql_fetch_assoc($result)) {
				$variant_id=$row["variant_id"];
				#$variant_flag_score=$row["flag_score"];
				#$variant_flag_label=$row["flag_label"];
				#$variant_flag_id=$row["flag_id"];
				#$variant_flag_definition=$row["flag_definition"];
				$variant_causal_score_table[$variant_id]=$row;
			};#foreach
			#print "<pre>";print_r($variant_causal_score_table);print "</pre>";


			# Validation SCORE
			$timestart=microtime(true);

			$sample_where_clause=(trim("$sample_id")!="")?"  AND variant_sample.sample_id!=$sample_id ":"";
			$manifest_where_clause=(trim("$manifest_id")!="")?" AND sample.manifest_id=$manifest_id ":"";

			$query="SELECT variant_sample.variant_id AS variant_id, variant_sample.sample_id AS sample_id, flags.label, flags.id
				FROM variant_sample
				LEFT OUTER JOIN comment ON (comment.id=variant_sample.id AND comment.type='variant_sample')
				LEFT OUTER JOIN variant ON (variant.id=variant_sample.variant_id)
				LEFT OUTER JOIN flags ON (flags.id=comment.flag_id AND flags.type=comment.type)
				LEFT OUTER JOIN sample ON (sample.id=variant_sample.sample_id)
				WHERE variant_sample.variant_id IN ($variant_ids)
				  $sample_where_clause
				  $manifest_where_clause
				GROUP BY variant_id, sample_id
				ORDER BY abs(flags.score) ASC";

			#print "<pre>$query</pre><BR>"; # AND sample.manifest_id=$manifest_id
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			$timestop=microtime(true);
			$timeDBexportVariantSampleALL=$timestop-$timestart;
			#echo "Time DBexport VariantSample ALL: $timeDBexportVariantSampleALL<BR>";

		   	$summary= "<table class=' '>\n";
			$validation_table=array();
			$flag_0="inconnue";
			$nb=0;
			while ($row = mysql_fetch_assoc($result)) {
				$nb++;
				#print "$nb "; #print "<BR>";
				#print_r($row); print "<BR>";
				$variant_id=$row['variant_id'];
				$sample_id2=$row['sample_id'];
				$label=($row['label']=="")?"$flag_0":$row['label'];
				$flag_score=$row['score']+0;
				$flag_label[$flag]=$label;
				#$value=str_replace(",","<BR>",$row['value']);
				#$comment=$row['comment'];
				$validation_table_all[$variant_id][$flag_score][$sample_id2]=$row;
				#$validation_table_variant_id[$variant_id]=$value;
			};#while
			#print "<BR>total: $nb <BR>"; #print "<BR>";
			#die();
			#ksort($validation_table);
			#echo "<PRE>"; print_r($validation_table_all); echo "</PRE>";
			/*
			$variant_checked_summary_tmp="";
			$VDScore=0;
			$validation_summary="<TABLE><TR  valign='top' halign='left' style=''>
						<TD width=''>Validation</TD>
						<TD width='10px'>&nbsp;</TD>
						<TD width=''>Number of variants</TD>
						<TD width='10px'>&nbsp;</TD>
						<TD width=''>Score added</TD>
						</TR>";
			foreach ($validation_table as $flag=>$sample_ids) {
				$label=$flag_label[$flag];

				$variant_checked_summary_tmp.=$summary;
				$sign= $flag < 0 ? -1 : ( $flag > 0 ? 1 : 0 ) ;
				$VDScore_label=($sign*(count($sample_ids)/(pow(10,abs($flag)-1))));
				$VDScore+=$VDScore_label;
				$validation_summary.="<TR><TD>$label</TD><TD></TD><TD>".count($sample_ids)." samples</TD><TD></TD><TD>$VScore_label</TD></TR>";
				#echo "$VScore+=($sign*(count($sample_ids)/(pow(10,abs($flag)-1))));<BR>";
			};#foreach
			#$VScoreLog=-log(2,10);
			#$VScoreLog=(1/$VScore);
			$variant_array["VDSCORE"]="$VDScore";
			#echo "<BR><B>".$variant_array["VARIANT_ID"].".$sample_id: VDScore=$VDScore (".sprintf('%1.2E',$VDScore).")</B><TABLE valign='left' halign='top'>$variant_checked_summary_header$variant_checked_summary_tmp</TABLE><BR>";
			$validation_summary.="</TABLE>";
			*/

		};#if


	# NB_GENE
	$command="grep -oiP 'symbol=\w+' $Variants_VCFContent | sort | uniq -c | wc -l ";
	$output_exec = shell_exec($command);
	#print "<pre>$output_exec</pre><BR>";
	$NB_GENE=$output_exec;

	# NB_VARIANT
	$command="grep -cv ^# $Variants_VCFContent";
	$output_exec = shell_exec($command);
	#print "<pre>$output_exec</pre><BR>";
	$NB_VARIANT=$output_exec;

	# NB_PASS
	$command="grep -oiP 'PZFlag=PASS' $Variants_VCFContent | grep ^# -cv";
	$output_exec = shell_exec($command);
	#print "<pre>$output_exec</pre><BR>";
	$NB_PASS=$output_exec;

	# NB_VARIANT
	$command="grep -oiP 'PZFlag=FILTERED' $Variants_VCFContent | grep ^# -cv";
	$output_exec = shell_exec($command);
	#print "<pre>$output_exec</pre><BR>";
	$NB_FILTERED=$output_exec;


	#grep -oiP 'PZFlag=PASS' validation/annotated.vcf | grep ^# -cv

	#break;

	$header_main_array=array();
	$header_array=array();
	#$Variants_VCFContent_split=split("\n",$Variants_VCFContent);
	#foreach (split("\n",$Variants_VCFContent) as $line_nb=>$line) {
	$nb_variant_in_html=0;
	$line_nb=0;
	$handle = fopen($Variants_VCFContent, "r");
	if ($handle) {
	    while (($line = fgets($handle)) !== false) {

		$line_nb++;



		#echo $line_nb;

		#$Variants_HTMLContent.= "$line<BR>\n";
		$line=trim($line);
		if ($line=="") {
			continue;
		};#if
		#$line = "#qeydtsdfth";

		$pattern0 = '/^##(.*)=(.*)/';
		preg_match($pattern0, $line, $matches0, PREG_OFFSET_CAPTURE);

		$pattern = '/^##(.*)=<(.*)>/';
		preg_match($pattern, $line, $matches, PREG_OFFSET_CAPTURE);
		#print_r($matches); print "<BR>";

		$pattern2 = '/^#/';
		preg_match($pattern2, $line, $matches2, PREG_OFFSET_CAPTURE);



		#print "<pre>LINE"; print $line; print "</pre>";

		if (!empty($matches)) {

			# VCF header


			#print "<pre>MATCHES"; print_r($matches); print "</pre>";

			$type=$matches[1][0];
			$def=$matches[2][0];
			#print "<pre>DEF"; print_r($def); print "</pre>";

			$def_split=explode(",",$def);
			$def_array=array();
			$varID="";
			foreach ($def_split as $k=>$varval) {
				$varval_split=explode("=",$varval,2);
				#print_r($varval_split);
				$var=$varval_split[0];
				$val=$varval_split[1];
				if ($var=="ID") {
					$varID=$val;
				} else {
					$def_array[$var]=$val;

					# Special annotation for database
					if ($var=="Description") {
						$patternS = '/^.*\[(.*)\]/';
						preg_match($patternS, $line, $matchesS, PREG_OFFSET_CAPTURE);
						#print "<pre>def_array"; print_r($matchesS); print "</pre>";
						if (!empty($matchesS)) {
							$defS=$matchesS[1][0];
							$defS_split=explode(";",$defS);
							#print "<pre>def_array"; print_r($defS_split); print "</pre>";
							foreach ($defS_split as $k=>$varvalS ) {
								$varvalS_split=explode("=",$varvalS,2);
								$varS=$varvalS_split[0];
								$valS=$varvalS_split[1];
								#print "<pre>"; print("$varS=>$valS"); print "</pre>";
								$def_array[$varS]=$valS;
							};#foreach
						};#if
					};#if
				};#if
			};#foreach
			#print "<pre>def_array"; print_r($def_array); print "</pre>";
			$header_array[$type][$varID]=$def_array;

			/*
			$pattern3 = '/ID=(.*),/';
			preg_match($pattern3, $matches[2][0], $matches3, PREG_OFFSET_CAPTURE);
			if (!empty($matches3)) {
				print "<pre>"; print_r($matches3); print "</pre>";
			};#if
			*/


		} elseif (!empty($matches0)) {
			#print "<pre>MAIN INFOS DEF HEADER</pre>";
			#print "<pre>MATCHES0"; print_r($matches0); print "</pre>";
			$header_main_array[$matches0[1][0]]=$matches0[2][0];

		} elseif (!empty($matches2)) {

			#print "<pre>HEADER_ARRAY1"; print_r($header_main_array); print "</pre>";
			#print "<pre>HEADER_ARRAY2"; print_r($header_array); print "</pre>";

			# VCF Variants Header
			#$Variants_HTMLContent.="HEADER $line";
			$variant_header=explode("\t",$line);

			#$Variants_HTMLContent.= "$line<BR>\n";

			if (!empty($matches)) {
				# VCF Header

			} else {



			};#if

		} else {






			# Find Variant Information
			$variant_array=VCFLinetoArray($line,$variant_header,$options);
			#print "<pre>"; print_r($variant_header); print "</pre>";
			#print "<pre>"; print_r($variant_array); print "</pre>";
			#$Variants_HTMLContent.="".join("\t",$variant_array)."<BR>";
			$nb_variant++;

			$variant_nb++;
			$variant_list["SYMBOL"][$variant_array["SYMBOL"]]++;
			$variant_list["FILTER_FLAG"][$variant_array["FILTER_FLAG"]]++;


			if ($nb_variant_in_html<=$options["limit"] || $options["limit"]=="") {

				if ($options["orderby"]=="") {
					$orderby="FILTER_SCORE";
				} else {
					$orderby=strtoupper($options["orderby"]);
				};#if

				# Gene filter
				if (0) {
					$gene_filter=0;
					if ($options["gene_filter"]=="") {
						$gene_filter=1;
					} else {
						foreach ($variant_array as $info=>$info_value) {
							if (!is_array($info_value)) {
								$gene_filter_array=explode(" ",$options["gene_filter"]);
								foreach (array($variant_array["SYMBOL"],$variant_array["HGVS"]) as $gene_filter_one) {
									#print "$info=>$info_value $gene_filter ";
									$gene_filter=($gene_filter || (strpos($info_value,$options["gene_filter"])!== false));
									#print " => $gene_filter <BR>";
								};#foreach
							};#if
						};#foreach
					};#if
				} else {
					$gene_filter=1;
				};#if

				# Global filter
				if (0) {
					$global_filter=0;
					if ($options["global_filter"]=="") {
						$global_filter=1;
					} else {
						foreach ($variant_array as $info=>$info_value) {
							if (!is_array($info_value)) {
								$global_filter_array=explode(" ",$options["global_filter"]);
								foreach ($global_filter_array as $global_filter_one) {
									#print "$info=>$info_value $global_filter ";
									$global_filter=($global_filter || (strpos($info_value,$global_filter_one)!== false));
									#print " => $global_filter <BR>";
								};#foreach
							};#if
						};#foreach
					};#if
				} else {
					$global_filter=1;
				};#if

				# hard filter
				if (0) {
					$hard_filter=0;
					if ($options["hardfiltering"]=="") {
						$hard_filter=1;
					} else {
						$hard_filter=$variant_array["FILTER_FLAG"]=="PASS";
					};#if
				} else {
					$hard_filter=1;
				};#if





				#print "global_filter=$global_filter<BR>";

				if (	1
					&& $hard_filter
					&& $global_filter
					&& $gene_filter
					#(!$options["hardfiltering"] || $variant_array["FILTER_FLAG"]=="PASS")
					#&& ($options["gene_filter"]==""
					#	|| strpos($variant_array["SYMBOL"],$options["gene_filter"])!== false
					#	|| strpos($variant_array["HGVS"],$options["gene_filter"])!== false
					#	)
				) {








						# variant in database
						if ($sample_id!="" && $variant_array["VARIANT_ID"]!="" && 1) {

							# variant flag
							$flag="0";
							$comment="";
							$variant_id=$variant_array["VARIANT_ID"];
							$timestart=microtime(true);
							$query="SELECT comment.id, comment.type, comment.flag_id, comment.comment, comment.history
								FROM comment
								WHERE comment.id='$variant_id.$sample_id'
								  AND comment.type='variant_sample'
								";
							#print "$query<BR>";
							$result = mysql_query($query,$mydbh) or die ("OUPS! ".mysql_error($mydbh)." <BR><PRE>$query</PRE>");
							$timestop=microtime(true);
							$timeDBexportVariantSampleComment=$timestop-$timestart;
							#echo "Time DBexport VariantSample Comment: $timeDBexportVariantSampleComment<BR>";
							$summary= "<table class=' '>\n";
							$comment_id="";
							$comment_type="";
							$flag="0";
							$rank="0";
							$comment="";
							$history="";
							$history_html="";
							while ($row = mysql_fetch_assoc($result)) {
								$comment_id=$row['id'];
								$comment_type=$row['type'];
								$comment_flag_id=$row['flag_id'];

								$comment=$row['comment'];
								$history=$row['history'];
								$history_html="<pre><small>$history</small></pre>";
								$summary.= "<tr>"
									."<td>".$row['flag_id']."</td>"
									."<td>".$row['comment']."</td>"
									."<td>".$row['history']."</td>"
								."</tr>";
								#echo "<PRE>$comment_id\n$comment_type\n$history</PRE><BR>";
							};#while
							#print "flag=$flag<BR>";
							$summary.= "</table>\n";
							#echo "<PRE>$query<BR>$summary</PRE>";

							$timestart=microtime(true);
							$query="SELECT id, type, score, rank, label, code, definition
								FROM flags
								WHERE type='variant_sample'
								ORDER BY rank ASC
								";
							#print "$query<BR>";
							$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
							$timestop=microtime(true);
							$timeDBexportVariantSampleFlag=$timestop-$timestart;
							#echo "Time DBexport VariantSample Flag: $timeDBexportVariantSampleFlag<BR>";
						   	$summary= "";
							$flags_options="";
							while ($row = mysql_fetch_assoc($result)) {
								$flag_score=$row['score'];
								$flag_rank=$row['rank'];
								$flag_id=$row['id'];
								#$flag_description=$row['label']." [".$row['level']."]";
								$flag_description=$row['label']; #." [".$row['level']."]";
								$summary.= "<tr>";
								$summary_head="<tr>";
								foreach ($row as $col=>$val) {
									$summary.= "<td>$val</td>";
									$summary_head.= "<td>$col</td>";
								};#foreach
								$summary.= "</tr>";
								$summary_head.= "</tr>";
								#$flags_array[$row['flag_value']]=$row['flag_label'];
								$flags_options.="<OPTION value='$flag_id' ".(($flag_id==$comment_flag_id)?"selected":"").">$flag_description [$flag_score]</OPTION>\n";
							};#while
							$summary= "<table class=' '>$summary_head $summary</table>\n";
							#echo $summary."<SELECT>$flags_options</SELECT>";

							# Validation SCORE
							/*
							$timestart=microtime(true);
							$query="SELECT variant_sample.variant_id AS variant_id, variant_sample.sample_id AS sample_id, variant_annotation.value, flags.label, flags.flag, comment.comment, sample.manifest_id
							FROM variant_sample
							LEFT OUTER JOIN comment ON (comment.id=variant_sample.id AND comment.type='variant_sample')
							LEFT OUTER JOIN variant ON (variant.id=variant_sample.variant_id)
							LEFT OUTER JOIN variant_annotation ON (variant_annotation.variant_id=variant.id)
							LEFT OUTER JOIN annotation ON (annotation.id=variant_annotation.annotation_id AND annotation.source='hgvs')
							LEFT OUTER JOIN flags ON (flags.flag=comment.flag AND flags.type=comment.type)
							LEFT OUTER JOIN sample ON (sample.id=variant_sample.sample_id)
							WHERE variant_sample.variant_id=$variant_id
							  AND variant_sample.sample_id!=$sample_id
							  AND sample.manifest_id=$manifest_id
							ORDER BY abs(comment.flag) ASC
							";

							#print "<pre>$query</pre><BR>";
							$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
						   	$timestop=microtime(true);
							$timeDBexportVariantSample=$timestop-$timestart;
							#echo "Time DBexport VariantSample: $timeDBexportVariantSample<BR>";
							$summary= "<table class=' '>\n";
							$validation_table=array();
							$flag_0="inconnue";
							while ($row = mysql_fetch_assoc($result)) {
								#print_r($row);
								#print_r($row["manifest_id"]);
								$variant_id=$row['variant_id'];
								$sample_id2=$row['sample_id'];
								$label=($row['label']=="")?"$flag_0":$row['label'];
								$flag=$row['flag']+0;
								$flag_label[$flag]=$label;
								#$value=str_replace(",","<BR>",$row['value']);
								#$comment=$row['comment'];
								$validation_table[$flag][$sample_id2]=$row;
								#$validation_table_variant_id[$variant_id]=$value;
							};#while
							*/


							$validation_table=$validation_table_all[$variant_id];

							ksort($validation_table);

							#echo "$variant_id<PRE>"; print_r($validation_table); echo "</PRE>";

							$variant_checked_summary_tmp="";
							$VDScore=0;
							$VDInternalFrequency=0;
							$validation_summary="";
							$variant_array["VDSCORE"]="0";
							$variant_array["VDIF"]=0;
							#echo "<pre>"; print_r($validation_table); echo "</pre>";
							if ($validation_table!="") {

								$validation_summary_flag="";
								$nb_sample_ids=0;
								foreach ($validation_table as $flag=>$sample_ids) {
									$label=$flag_label[$flag];
									$nb_sample_ids+=count($sample_ids);
									$variant_checked_summary_tmp.=$summary;
									$sign= $flag < 0 ? -1 : ( $flag > 0 ? 1 : 0 ) ;
									$VDScore_label=($sign*(count($sample_ids)/(pow(10,abs($flag)-1))));
									$VDScore+=$VDScore_label;
									#$validation_summary_flag.="<TR><TD>$label</TD><TD></TD><TD>".count($sample_ids)." samples</TD><TD></TD><TD>$VDScore_label</TD></TR>";
									#$validation_summary_flag.="<TR><TD>".count($sample_ids)." </TD><TD>$label samples</TD><TD></TD><!--<TD></TD><TD>$VDScore_label</TD>--></TR>";
									$validation_summary_flag.=(($validation_summary_flag=="")?"":", ").count($sample_ids)." $label";
									#echo "$VScore+=($sign*(count($sample_ids)/(pow(10,abs($flag)-1))));<BR>";
								};#foreach
								#$VScoreLog=-log(2,10);
								#$VScoreLog=(1/$VScore);
								$variant_array["VDSCORE"]="$VDScore";
								$VDInternalFrequency=$nb_sample_ids/$nb_sample_manifest;
								$variant_array["VDIF"]=$VDInternalFrequency;
								#echo "<BR><B>".$variant_array["VARIANT_ID"].".$sample_id: VDScore=$VDScore (".sprintf('%1.2E',$VDScore).")</B><TABLE valign='left' halign='top'>$variant_checked_summary_header$variant_checked_summary_tmp</TABLE><BR>";
								$validation_summary="
									This variant was identified in $nb_sample_ids/$nb_sample_manifest other sample".(($nb_sample_manifest>1)?"s":"")."
									(Internal Frequency IF=".round($variant_array["VDIF"]*100,2)."%)<BR>

									$validation_summary_flag".(($validation_summary_flag=="")?"":"<BR>");
								if ($nb_sample_ids>0) {
								/*$validation_summary.="
									<!--<TABLE>
										<TR  valign='top' halign='left' style=''>
											<TD width=''>Validation</TD>
											<TD width='10px'>&nbsp;</TD>
											<TD width=''>Number of variants</TD>
											</TR>-->
										$validation_summary_flag2
									</TABLE>-->";*/
								};#if
							} else {

								$validation_summary="
									This variant was identified in 0/$nb_sample_manifest other sample".(($nb_sample_manifest>1)?"s":"")."
									(Internal Frequency IF=0%)<BR>
									";

							};#if
						};#if







						# SNAPSHOT
						if ($snapshot_array[$variant_array["CHROM"]][$variant_array["POS"]]!="") {
							$snapshot_table="<BR><TABLE><TR><TD>Snapshots:</TD><TD width='20'></TD>";
							#print "CHROM:".$variant_array["CHROM"]."<BR>";
							#print "POS:".$variant_array["POS"]."<BR>";
							foreach ($snapshot_array[$variant_array["CHROM"]][$variant_array["POS"]] as $aligner=>$snapshot_id) {
								#$snapshot_array[$snapshot_chrom][$snapshot_pos][$vcf_aligner]=$snapshot_id
								#print "<pre>"; print_r($row); print "</pre>";
								$URL="blob_image.php?id=$snapshot_id";
								$LINK="<span class=' ' onclick=\"javascript: var newWin = window.open('$URL','SNAPSHOT','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=650, height=800, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>$aligner</span>";
								$snapshot_table.="<TD>$LINK</TD><TD width='20'></TD>";

							};#foreach
							$snapshot_table.="</TR></TABLE>";
						};#if





						# Annotation
						$annotation_table="<TABLE>";
						foreach ($variant_array["INFO"] as $annotation_source=>$annotation_value) {
							#if (in_array($annotation_source,$annotation_list)) {
								$annotation_table.="<TR><TD valign='top' title=\"".str_replace("\"","",$header_array["INFO"][$annotation_source]["Description"])."\">$annotation_source</TD><TD width='20'></TD><TD valign='top'>".str_replace(",",", ",$annotation_value)."</TD></TR>\n";
							#};#if
						};#foreach
						$annotation_table.="</TABLE>";

						# Annotation classified
						$annotation_array=array();
						foreach ($variant_array["INFO"] as $annotation_source=>$annotation_value) {
							$annotationType=$header_array["INFO"][$annotation_source]["AnnotationType"];
							$annotationType=($annotationType=="")?"unknown":$annotationType;
							$annotation_array[$annotationType][$annotation_source]="<TD valign='top' title=\"".str_replace("\"","",$header_array["INFO"][$annotation_source]["Description"])."\">$annotation_source</TD><TD width='20'></TD><TD valign='top'>".str_replace(",",", ",$annotation_value)."</TD>";
							#if (in_array($annotation_source,$annotation_list)) {
							#$annotation_table.="<TR><TD valign='top' title=\"".str_replace("\"","",$header_array["INFO"][$annotation_source]["Description"])."\">$annotation_source</TD><TD width='20'></TD><TD valign='top'>".str_replace(",",", ",$annotation_value)."</TD></TR>\n";
							#};#if
						};#foreach
						$annotation_table_known="";
						$annotation_table_unknown="";
						foreach ($annotation_array as $annotationType=>$annotationType_array) {
							if ($annotationType!="unknown") {
								$annotation_table_known.="<TR><TD colspan='10'><B>".strtoupper($annotationType)."</B></TD></TR>\n";
								foreach ($annotationType_array as $annotation_source=>$annotation_source_content) {
									$annotation_table_known.="<TR>$annotation_source_content</TR>\n";
								};#foreach
							} else {
								$annotation_table_unknown.="<TR><TD colspan='10'><B>".strtoupper($annotationType)."</B></TD></TR>\n";
								foreach ($annotationType_array as $annotation_source=>$annotation_source_content) {
									$annotation_table_unknown.="<TR>$annotation_source_content</TR>\n";
								};#foreach
							};#if
						};#foreach

						#if () {

						#};#if
						$annotation_table="<TABLE>$annotation_table_known$annotation_table_unknown</TABLE>";
						#$annotation_table_unknown="<TABLE>$annotation_table_unknown</TABLE>";







						# Samples
						if (count($variant_array["SAMPLES"])>0) {
							$sample_table="";
							$sample_table_header="";
							$sample_table_lines="";
							# create header
							$sample_table_header_array=array();
							foreach ($variant_array["SAMPLES"] as $sample=>$qual_array) {
								foreach ($qual_array as $qual_name=>$qual_value) {
									$sample_table_header_array[$qual_name]=1;
								};#foreach
							};#foreach
							#print "<PRE>"; print_r($sample_table_header_array); print "</PRE>";

							foreach ($variant_array["SAMPLES"] as $sample=>$qual_array) {
								$sample_table_header_tmp="<TR><TD></TD>";
								$sample_table_lines.="<TR><TD>$sample</TD>";
								#foreach ($qual_array as $qual_name=>$qual_value) {
								foreach ($sample_table_header_array as $qual_name=>$qual_value_header) {
									$qual_value=$qual_array[$qual_name];
									$sample_table_header_tmp.="<TD width='20'></TD><TD>$qual_name</TD>";
									$sample_table_lines.="<TD width='20'></TD><TD>$qual_value</TD>";
								};#foreach
								$sample_table_header_tmp.="</TR>";
								if (count($qual_array)>1) {
									$sample_table_header=$sample_table_header_tmp;
								};#if
								$sample_table_lines.="</TR>";
							};#foreach
							$sample_table="<TABLE>$sample_table_header$sample_table_lines</TABLE>";
						};#if





						## Validation form

						# Validation
						# <!--<option value=''>Flag the variant...</option>-->
						$validation="";
						$validation_width="300px";
						#print date(DATE_ATOM);
						#print date_default_timezone_set('UTC');

						$validation_form=0;
						if ($sample_id!="" && $variant_array["VARIANT_ID"]!="") {
							$validation_form=1;
						};#if

						if ($validation_form) {
							$variant_id=$variant_array["VARIANT_ID"];
							$validation="	<iframe id='ponyo_frame' name='iframe_flags_$nb_variant' border='0' style='visibility:hidden;float:right;display:none;'></iframe>
									<form id='form_flags_$nb_variant' method='post' action='validation.php' target='iframe_flags_$nb_variant' style='height:0px;visibility:;display:;'>
										<TABLE width='$validation_width' border=0 class='' style='padding:0px;margin-top:0px;'>
											<TR>

												<TD valign='top' width='99%'>
													<SELECT name='flag_id' style='margin-bottom:-5px;' width='99%'>$flags_options</SELECT>
												</TD>
												<TD valign='top'>
													<input type='submit' id='validation' name='validation' value='Save' class='btn button search' onclick=\"javascript:;\">
													<input data-original-title='' value='sample_id' type='hidden'>
													<input type='hidden' id='type' name='type' style='visibility:hidden;float:right;display:none;' value='variant_sample'></input>
													<input type='hidden' id='project_id' name='project_id' style='visibility:hidden;float:right;display:none;' value='$project_id'></input>
													<input type='hidden' id='run' name='run' style='visibility:hidden;float:right;display:none;' value='$input_run'></input>
													<input type='hidden' id='id' name='id' style='visibility:hidden;float:right;display:none;' value='".$variant_array["VARIANT_ID"].".".$sample_id."'></input>
													<input type='hidden' id='author' name='author' style='visibility:hidden;float:right;display:none;' value='".USERNAME."'></input>
													<input data-original-title='' onclick=\"document.getElementById('form_flags_$nb_variant').submit();\" type='hidden'>
												</TD>
												<TD valign='top'>
													<input type='button' id='history' name='history' value='H' class='btn button search' onclick=\"javascript:;document.getElementById('variant_validation_history_$nb_variant').src='validation_history.php?id=$variant_id.$sample_id&type=variant_sample&project_id=$project_id&rand='+Math.floor((Math.random()*10000000)+1);show_hide('variant_validation_history_$nb_variant');\"></TD>

											</TR>
											<TR>

												<TD colspan=3 align='center' valign='top'>
													<TEXTAREA id='comment' name='comment' class='' style='border: none; width: 99%; height: 100%; -webkit-box-sizing: border-box; -moz-box-sizing: border-box; box-sizing: border-box; padding:0px;' rows=3>$comment</TEXTAREA>
												</TD>

											</TR>
											<TR>
												<TD colspan=3 align='left' valign='right' zalign='top' >
													<iframe id='variant_validation_history_$nb_variant' name='variant_validation_history_$nb_variant' border='1' style='float:left;text-align:right;display:none;z-index:2;position:absolute;' marginwidth='0' marginheight='0' frameborder='0' vspace='0' hspace='0' width='$validation_width'>$history_html</iframe>
												</TD>

											</TR>
										</TABLE>
									</form>

							";
						};#if










						# Ensembl Link
						$ensembl_prefix="grch37"; # www by default
						$ensembl_links="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Location/View?db=core;r=".$variant_array["REGION"]."' target='ensembl'>Location</a></TD><TD width='20'></TD>";
						if ($variant_array["ENSEMBL"]!="") {
							$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=".$variant_array["ENSEMBL"]."' target='ensembl'>Gene</a></TD><TD width='20'></TD>";
						};#if
						if ($variant_array["ENSEMBL"]!="" && $variant_array["REGION"]!=""
							&& ( $variant_array["DBSNP"]!="" || $variant_array["COSMIC"]!="" || $variant_array["HGMD"]!="" )) {
							#print "$v_id=(".$variant_array["DBSNP"]." ".$variant_array["HGMD"]." ".$variant_array["COSMIC"]."";
							if ($variant_array["DBSNP"]!="") {
								$v_id=$variant_array["DBSNP"];
							} elseif ($variant_array["HGMD"]!="") {
								$v_id=$variant_array["HGMD"];
							} elseif ($variant_array["COSMIC"]!="") {
								$v_id=$variant_array["COSMIC"];
							};#if
							#$v_id=($variant_array["DBSNP"]!="")?$variant_array["DBSNP"]:($variant_array["HGMD"]!="")?$variant_array["HGMD"]:($variant_array["COSMIC"]!="")?$variant_array["COSMIC"]:"";
							#print "v_id: $v_id<BR>";
							$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Variation/Mappings?db=core;gene=".$variant_array["ENSEMBL"].";v=$v_id;r=".$variant_array["REGION"].";vdb=variation' target='ensembl'>Variation</a></TD><TD width='20'></TD>";
							$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Variation/Context?db=core;gene=".$variant_array["ENSEMBL"].";v=$v_id;r=".$variant_array["REGION"].";vdb=variation' target='ensembl'>Context</a></TD><TD width='20'></TD>";
							$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Variation/Phenotype?db=core;gene=".$variant_array["ENSEMBL"].";v=$v_id;r=".$variant_array["REGION"].";vdb=variation' target='ensembl'>Phenotype</a></TD><TD width='20'></TD>";
							if ($variant_array["DBSNP"]!="") {
								$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Variation/Population?db=core;gene=".$variant_array["ENSEMBL"].";v=$v_id;r=".$variant_array["REGION"].";vdb=variation' target='ensembl'>Population</a></TD><TD width='20'></TD>";
							};#if

						};#if
						if ($variant_array["ENSEMBL"]!="") {
							$ensembl_links.="<TD><a href='http://$ensembl_prefix.ensembl.org/Homo_sapiens/Gene/Variation_Gene/Table?db=core;gene=".$variant_array["ENSEMBL"].";v=$v_id;r=".$variant_array["REGION"].";vdb=variation' target='ensembl'>Variation_Gene</a></TD><TD width='20'></TD>";
						};#if





						#echo "<pre>"; print_r($variant_list["SYMBOL"]); echo "</pre>";
						$variant_id=$variant_array["VARIANT_ID"];
						#.sprintf('%1.2E',$VScore).
						#$order_by="SYMBOL";
						#echo "$orderby ".$variant_array[$orderby]." ".($variant_array[$orderby]+1)."<BR>";
						#print_r($variant_array); echo "<BR>";
						#$Variants_HTMLContent[$variant_array["FILTER_SCORE"]].="
						#".$variant_array["DBSNP"]."



						# GROUP DATABASE
						#echo "<PRE>";
						$group_database_annotation="";
						if (0) {
							foreach ($user_groups as $group_id=>$group_ref) {
								#echo "# GROUP ID: $group_id\n";
								#echo "# GROUP REF: $group_ref\n";
								#print_r($variant_array["INFO"]);
								if (trim($variant_array["INFO"]["database_$group_ref"])!="") {
									#echo "# ".$variant_array["INFO"]["database_$group_ref"]."\n";
									$sep=(($group_database_annotation)=="")?"":". ";
									$group_database_annotation.=$sep.trim($variant_array["INFO"]["database_$group_ref"]);
								};#if
							};#foreach
						};#if
						#echo "</PRE>";
						#&nbsp;&nbsp;&nbsp;<span title='Variant Depth'>DP=".$variant_array["DP"]."</span>
						#&nbsp;&nbsp;&nbsp;<span title='Variant Depth'>DP=".$depthbed_stats_array[str_replace("chr","",$variant_array["CHROM"])][$variant_array["POS"]]["min"]."</span>



						#print $variant_array["HGVS"];
						$customHGVS=find_hgvs($variant_array["HGVS"],$variant_array["NOMEN"],$variant_array["SYMBOL"],$user_groups,$customNM_infos);
						if ($customHGVS == "" && $variant_array["NOMEN"] != "") {
							$customHGVS=explode(",",$variant_array["NOMEN"])[0];
						} else if ($customHGVS == "" && $variant_array["HGVS"] != "") { #hgvsEnsembl
							$customHGVS=explode(",",$variant_array["HGVS"])[0];
						} else if ($customHGVS == "" && $variant_array["HGVSENSEMBL"] != "") { #hgvsEnsembl
							$customHGVS=explode(",",$variant_array["HGVSENSEMBL"])[0];
						} else if ($customHGVS == "" && $variant_array["SNPEFF_HGVS"] != "") { #hgvsEnsembl
							$customHGVS=explode(",",$variant_array["SNPEFF_HGVS"])[0];
						};#if
						#print "<PRE>!$customHGVS</PRE>";
						#print "<PRE>!".$variant_array["NOMEN"]."</PRE>";
						$customHGVS=str_replace(",","<BR>",$customHGVS);

						$variant_color="gray";
						$variant_validation="";
						#$variant_score_max=5;
						#$variant_score_min=0;
						$variant_score=($variant_causal_score_table[$variant_id]["flag_score"]!="")?$variant_causal_score_table[$variant_id]["flag_score"]:"?"; #print_r($variant_causal_score_table[$variant_id]);
						$variant_score_definition=($variant_causal_score_table[$variant_id]["flag_definition"]!="")?$variant_causal_score_table[$variant_id]["flag_definition"]:$variant_score_default["flag_definition"]; #print_r($variant_causal_score_table[$variant_id]);
						#print_r($variant_score);
						$variant_color_scale_max=512;
						$variant_color_scale_min=1;

						if (!isset($variant_score_max) ) { $variant_score_max=521; };
						if (!isset($variant_score) || $variant_score=="?") { $variant_score=0; };
						$variant_color_code=($variant_score*$variant_color_scale_max)/$variant_score_max;
						$variant_color_code_inv=round($variant_color_scale_max-$variant_color_code);
						$variant_color_codeR=(($variant_color_code_inv-255)>255)?255:$variant_color_code_inv-255;
						$variant_color_codeG=(($variant_color_code_inv-255)>255)?255:$variant_color_code_inv-255;
						$variant_color_codeB=($variant_color_code_inv>255)?255:$variant_color_code_inv;
						$variant_color="rgb($variant_color_codeR,$variant_color_codeG,$variant_color_codeB)";
						$variant_validation="<DIV  onclick=\"javascript:show_hide('variant_content_$nb_variant');show_hide('variant_links_$nb_variant');\"
													style='halign=middle;float:right;cursor:pointer;background-color:$variant_color;height:15px;width:15px;'
													title='Variant Causative Score'></DIV>";
						$variant_validation="";



						$variant_frequency_tag="Unknown";
						$variant_rare_frequency=0.05;
						$variant_error_frequency=0.1;
						if ( abs(($variant_array["POPFREQ"]+0)+($variant_array["VDIF"]+0))+0 != 0 ) {
							$COV_IF_EF=(abs($variant_array["POPFREQ"]-$variant_array["VDIF"])/abs($variant_array["POPFREQ"]+$variant_array["VDIF"]))+0;
						} else {
							$COV_IF_EF="?";
						};#if
						if ($variant_array["DBSNP_NONFLAGGED"] != "") {
							$variant_frequency_tag="Polymorphism";
							if (	($variant_array["VDIF"] > $variant_error_frequency && $variant_array["POPFREQ"] < $variant_error_frequency)
								|| ($variant_array["VDIF"] > $variant_error_frequency && $variant_array["POPFREQ"] > $variant_error_frequency && $COV_IF_EF>0.5)
							) {
								$variant_frequency_tag.="/Error?".$variant_array["EF"];
							#} else {
							#	$variant_frequency_tag="Polymorphism";
							};#if
						} elseif ($variant_array["DBSNP"] != "" && $variant_array["VDIF"] <= $variant_error_frequency) {
							$variant_frequency_tag="Rare/Pathogenic";
						} elseif ($variant_array["VDIF"] > $variant_error_frequency) {
							$variant_frequency_tag="Error?";
						};#if
						#$variant_frequency_tag.=$COV_IF_EF;
						$variant_frequency_tag_title="Variant type: Based on dbSNP annotation and Internal and External frequencies (Briefly, polymorphism of annotated with dbSNP NonFlagged, Rare/pathogenic is annotated in dbSNP but not as a polymorphism, Error if IF and EF inconsistante. Rare max frequency = $variant_rare_frequency, Error min frequency = $variant_error_frequency)";

						$nb_variant_in_html++;


						$variant_array["CHROMNUM"]=str_replace("chr", "", $variant_array["CHROM"]);


						#print "$orderby";
						#if ($variant_nb == "") { $variant_nb=1; };
						$Variants_HTMLContent[$variant_array[$orderby]][$variant_nb]="

							<SCRIPT>
							function show_hide(id) {
								if(document.getElementById(id).style.display=='none') {
									document.getElementById(id).style.display='inline';
								} else {
									document.getElementById(id).style.display='none';
								}
								return true;
							}
							</SCRIPT>
							<table border=0 width='100%' class='' ><TR><TD>


								<table border=0 width='100%' class='variant_mini_header' >
									<TR>
										<TD class=' variant_header' valign='top' halign='left' height='' width='350'>
											<div  style='float:;' height='100%'>

												<DIV  onclick=\"javascript:show_hide('variant_content_$nb_variant');show_hide('variant_links_$nb_variant');\"
													style='float:right;cursor:pointer;'
													title='show or hide variant content'>[+]</DIV>

												<b>".(($variant_array["SYMBOL"]=="")?"unknown":$variant_array["SYMBOL"]."<BR>")."</b>
												$variant_validation
												<!--[#id:$variant_id]--> ".$variant_array["CHROMPOS"]."<br>
												".(($customHGVS=="")?"":"$customHGVS<BR>")."
												<DIV>
												<!--".($variant_array["DBSNP"]==""?"":" &nbsp;&nbsp;&nbsp;<span title='Polymorphism'>".(($variant_array["DBSNP"] != "" && $variant_array["DBSNP_NONFLAGGED"] == "") ? " [Rare/Pathogenic]" : "[SNP]")." </span>")."-->
												</DIV>


											</div>
										</TD>
										<TD valign='top' align='left'  class=' variant_mini_header' width=''>
											<div class=''style='float:middle;' height='100%'>
											<b>".$variant_array["FILTER_FLAG"]."</b>
											&nbsp;&nbsp;&nbsp;<span title='$variant_frequency_tag_title'>[$variant_frequency_tag]</span>
											&nbsp;&nbsp;&nbsp;<span title='PrioritiZation Score'>PZ=".$variant_array["FILTER_SCORE"]."</span>
											&nbsp;&nbsp;&nbsp;<span title='Variant Depth'>DP=".$variant_array["DP"]."</span>
											&nbsp;&nbsp;&nbsp;<span title='Genotype Quality'>GQ=".$variant_array["GQ"]." </span>
											&nbsp;&nbsp;&nbsp;<span title='Allele Frequency'>VAF=".round($variant_array["AlleleFrequency"]*100,2)."%</span>
											&nbsp;&nbsp;&nbsp;<span title='Variant Internal Frequency'>IF=".round($variant_array["VDIF"]*100,2)."%</span>
											&nbsp;&nbsp;&nbsp;<span title='Variant External Frequency'>EF=".round($variant_array["POPFREQ"]*100,2)."%</span>
											&nbsp;&nbsp;&nbsp;<span title='Variant Causative Score: $variant_score_definition [$variant_score]'>CS=$variant_score</span>

											<!--&nbsp;&nbsp;&nbsp;<span title='ValiDation Score Score'>VD=".((abs($variant_array["VDSCORE"])>0.001||$variant_array["VDSCORE"]==0)?$variant_array["VDSCORE"]:sprintf('%1.2E',$variant_array["VDSCORE"]))."</span>-->
											<!--&nbsp;&nbsp;&nbsp;<span title=''>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;</span>-->



											".($group_database_annotation==""?"":((strlen($group_database_annotation)>4)?"<BR><span title='Group database annotation'>$group_database_annotation</span>":"&nbsp;&nbsp;&nbsp;<span title='Group database annotation'>[$group_database_annotation]</span>"))."
											<BR>
											".$variant_array["FILTER_COMMENT"]."
											<BR>
											</div>

										</TD>".((!$validation_form)?"":"
										<TD valign='top' align='right' width='$validation_width' class=' variant_validation' style='backgroud-color:white;'>
											<div style='float:right;text-align:right;display:;margin-left:0px;margin-top:;' id='variant_validation_$nb_variant'>
												$validation
											</div>
										</TD>")."
									</TR>

								</TABLE>
							</TD></TR>
							</TABLE>
							<table border=0 width='100%' class=' ' >
							<TR id='variant_content_$nb_variant' style='display:none;'>
							<!--<TD width='20px'></TD>-->
							<TD width='99%' class=''  >
							<TABLE border=0 width='100%' height='100%' class='variant_infos'><TR>
								<TD width='20px'></TD>
								<TD >
								$snapshot_table
								<BR>
								$sample_table
								<!--<A HREF='".str_replace(">","|","search.php?s=Search&q=".$variant_array["CHROMPOS"])."'>Variant in the database</A><BR>-->
								<!--<BR>$validation_summary-->
								<!--<BR>-->
								<!--".$depthbed_stats_array[str_replace("chr","",$variant_array["CHROM"])][$variant_array["POS"]]["comment"]."<BR>-->
								<BR>
								<TABLE>
									<TR>
										<TD>Ensembl</TD>
										<TD width='20'></TD>
										$ensembl_links

									</TR>
									<TR>
										<TD>Links</TD>
										<TD width='20'></TD>
										<TD><a href='http://exac.broadinstitute.org/variant/".$variant_array["CHROM"]."-".$variant_array["POS"]."-".$variant_array["REF"]."-".$variant_array["ALT"]."' target='exac'>ExAC</a></TD><TD width='20'></TD>
										<TD><a href='http://www.ncbi.nlm.nih.gov/clinvar?term=(".str_replace("chr","",$variant_array["CHROM"])."[Chromosome]) AND ".$variant_array["POS"]."[Base Position]' target='clinvar'>ClinVar</a> </TD><TD width='20'></TD>
										".(($variant_array["DBSNP"]!="")?"<TD><a href='http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=".$variant_array["DBSNP"]."' target='dbsnp'>dbSNP</a></TD><TD width='20'></TD>":"")."
										".(($variant_array["DBSNP"]!="")?"<TD><a href='https://www.snpedia.com/index.php/".$variant_array["DBSNP"]."' target='SNPedia'>SNPedia</a></TD><TD width='20'></TD>":"")."
										".(($variant_array["DBSNP"]!="")?"<TD><a href='http://www.ncbi.nlm.nih.gov/pubmed?term=".$variant_array["SYMBOL"]." AND ".$variant_array["PNOMEN"]."' target='PubMed'>PubMed</a></TD><TD width='20'></TD>":"")."
										".(($variant_array["DBSNP"]!="")?"<TD><a href='http://www.ncbi.nlm.nih.gov/variation/view/?filters=source:dbsnp&q=".$variant_array["DBSNP"]."' target='VariationViewer'>Variation</a></TD><TD width='20'></TD>":"")."
										".(($variant_array["GNOMEN"]!="" && $variant_array["PNOMEN"]!="")?"<TD><a href='https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/LitVar/#!?query=".$variant_array["GNOMEN"]."%20".$variant_array["PNOMEN"]."' target='LitVar'>LitVar</a></TD><TD width='20'></TD>":"")."
										".(($variant_array["GNOMEN"]!="" && $variant_array["PNOMEN"]!="")?"<TD><a href='https://cancer.sanger.ac.uk/cosmic/search?q=".$variant_array["GNOMEN"]."+".$variant_array["PNOMEN"]."' target='COSMIC'>COSMIC</a></TD><TD width='20'></TD>":"")."
										".(($variant_array["GNOMEN"]!="" && $variant_array["PNOMEN"]!="")?"<TD><a href='https://varsome.com/variant/hg19/".$variant_array["GNOMEN"]."%20".$variant_array["PNOMEN"]."' target='VarSome'>VarSome</a></TD><TD width='20'></TD>":"")."
										<TD width='20'></TD>
									</TR>
									<TR><TD></TD>
										<TD width='20'></TD>
										".(($variant_array["GNOMEN"]!="")?"<TD><a href='https://www.ncbi.nlm.nih.gov/gap/advanced_search/?TERM=".$variant_array["GNOMEN"]."' target='dbGaP'>dbGaP</a></TD><TD width='20'></TD>":"")."
										".(($variant_array["SYMBOL"]!="")?"<TD><a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=".$variant_array["SYMBOL"]."' target='genecards'>GeneCards</a></TD>
										<TD width='20'></TD>":"")."
										".(($variant_array["CHROM"]!="" && $variant_array["POS"]!="" && $variant_array["ALT"]!="" && $variant_array["REF"]!="")?"<TD><a href='http://www.fudan-pgx.org/premedkb/index.html#/search/result?queryType=3&step=1&num=1&term=%27".$variant_array["CHROMNUM"]."-".$variant_array["POS"]."-".$variant_array["REF"]."-".$variant_array["ALT"]."%27%5Bvariant%5D' target='PreMedKB'>PreMedKB</a></TD>
										<TD width='20'></TD>":"")."
										".(($variant_array["DBSNP"]!="")?"<TD><a href='https://www.pharmgkb.org/rsid/".$variant_array["DBSNP"]."' target='PharmGKB'>PharmGKB</a></TD><TD width='20'></TD>":"")."
										<TD width='20'></TD>

									</TR>
								</TABLE>
								<BR>
								$annotation_table

								</TD></TR></TABLE>

							</TD></TR></TABLE>
						";



				};#if global filter


			} else {

				#echo "nb variant $nb_variant_in_html sort de la limit ".$options["limit"];
				break;

			};#if


		};#if
	#};#foreach

	    }; # while
	    fclose($handle);
	} else {
	    // error opening the file.
		print "error opening file";
	}; # if



	# Summary
	#$NB_GENE=(count($variant_list["SYMBOL"])+0)." gene".(((count($variant_list["SYMBOL"])+0)>1)?"s":"");
	#$NB_VARIANT=($variant_nb+0)." variant".((($variant_nb+0)>1)?"s":"");
	$NB_GENE_SHOW=($NB_GENE+0)." gene".((($NB_GENE+0)>1)?"s":"");
	$NB_VARIANT_SHOW=($NB_VARIANT+0)." variant".((($NB_VARIANT+0)>1)?"s":"");
	#$NB_PASS_SHOW=($NB_VARIANT+0)." variant".((($NB_VARIANT+0)>1)?"s":"");
	#$NB_FILTERED_SHOW=($NB_VARIANT+0)." variant".((($NB_VARIANT+0)>1)?"s":"");

	$summary="<!--<div class=' variant_header' style='clear:left;width:100%'>-->
			<table border=0 width='100%' class='  variant_header block' >
			<TR>
			<TD width='20px'></TD>
			<TD>
			<BR>
			".$NB_GENE_SHOW."&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
			".$NB_VARIANT_SHOW."&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;
			".($NB_PASS+0)." PASS &nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp; ".($NB_FILTERED+0)." FILTERED
			<BR><BR></TD></TR>
			</TABLE>

		";
		# $variant_list["FILTER_FLAG"]["FILTERED"]+$variant_list["FILTER_FLAG"]["FILTERED,PASS"]
	# ASC/DESC sorting
	if (0) {
	if ($options["ascdesc"]=="DESC") {
		krsort($Variants_HTMLContent);
	} elseif ($options["ascdesc"]=="ASC") {
		ksort($Variants_HTMLContent);
	};#if
	};#if

	$v_nb=0;
	$return=$summary."";
	$limited=0;
	foreach ($Variants_HTMLContent as $korder=>$Variants_list) {
		foreach ($Variants_list as $Variant_id=>$Variant_profile) {
			$v_nb++;
			if ($v_nb<=$options["limit"] || $options["limit"]=="") {
				$return.=$Variant_profile;
			} elseif (!$limited) {
				$return.="
						<table border=0 width='100%' class='  variant_header block' >
						<TR>
						<TD width='20px'></TD>
						<TD>

						... List of variants limited to ".$options["limit"]."
						</TD></TR>
						</TABLE>


					";
				# (out of ".($variant_nb+0).")
				$limited=1;
			};#if
		};#foreach
	};#foreach

	#$return=$summary.implode("",$Variants_HTMLContent);

	return $return;

	#return $Variants_HTMLContent;

}#function VCFFiletoHTML

function ListCaller($Directory,$ListFiles){
  	$output="";
  	$MyDirectory = opendir($Directory); # or die('Erreur');
	$BAM="";
	$Callers="";
	while($Entry = @readdir($MyDirectory)) {
		$mod1 = preg_grep("/^.*\.bam$/",array($Entry));
		if (count($mod1)>0) {
			$BAM.="<A>$Entry</A>";
		};#if
		if(is_dir($Directory.'/'.$Entry)&& $Entry != '.' && $Entry != '..') {
            		$Callers.= '<ul>'.$Entry;
			if ($ListFiles) {
				$Callers.= '<ul>';
				$Callers.=ListFiles($Directory.'/'.$Entry);
				$Callers.= '</ul>';
			};#if
			$Callers.= '</ul>';
		}
		else {
			#echo '<li>'.$Entry.'</li>';
        	};#if
	}
	closedir($MyDirectory);
	return "<ul>".$BAM."</ul>".$Callers;
}

function ListFiles($Directory,$ListFiles){
  	$output="";
	$MyDirectory = opendir($Directory); # or die('Erreur');
	$VCF="";
	$Files="";
	while($Entry = @readdir($MyDirectory)) {
		$mod1 = preg_grep("/^.*\.vcf$/",array($Entry));
		if (count($mod1)>0) {
			$VCF.="<A>$Entry</A> ";
		}
		#$VCF.=$Directory.''.$Entry;
		if(is_dir($Directory.'/'.$Entry)&& $Entry != '.' && $Entry != '..') {
            		#$output.= '<ul>'.$Entry.'</ul>';
			$MyDirectory2 = opendir($Directory.'/'.$Entry); # or die('Erreur');
			while($Entry2 = @readdir($MyDirectory2)) {
				$mod2 = preg_grep("/^.*\.vcf$|txt$/",array($Entry2));
				if (count($mod2)>0) {
					$VCF.="<A>$Entry/$Entry2</A> ";
				}
			}
		}
		else {
			#$output.= '<li>'.$Entry.'</li>';
        	};#if
	}
	closedir($MyDirectory);
	return "<ul>".$VCF."</ul>";
}

function match($what,$where) {
  $x=file_get_contents($where);
  $a=explode("\n",$x);
  $b=preg_grep("/^{$what}$/",$a);
  if (count($b)!=1) return FALSE;
  $b=implode("",$b);
  $b=explode(":",$b);
  if (count($b)!=2) return FALSE;
  $b=$b{1};
  return $b;
}

function find_thing($what,$where)
 {
 $temp=file_get_contents($where);
 if (!$start=strpos($temp,$what))
  return false;
 else
  $start+=strlen($what)+1; //Add +1 because of the :
 $end=strpos($temp,"\n",$start);
 return substr($temp,$start,$end);
 }

function standard_deviation($aValues, $sample = false)
{
	if ($sample && $aValues[$sample]!="") {
		$fMean = (array_sum($aValues) - $aValues[$sample]) / (count($aValues) - 1);
	} else {
		$fMean = array_sum($aValues) / count($aValues);
	};#if
	$fVariance = 0.0;
	foreach ($aValues as $key=>$i)
	{
		#echo "$key<BR>";
		if ($sample && $sample==$key) {

		} else {
			$fVariance += pow($i - $fMean, 2);
		};#if
	}
	$fVariance /= ( $sample ? count($aValues) - 1 : count($aValues) );
	return (float) sqrt($fVariance);
}

function standard_deviation_old($aValues, $bSample = false)
{
    $fMean = array_sum($aValues) / count($aValues);
    $fVariance = 0.0;
    foreach ($aValues as $i)
    {
        $fVariance += pow($i - $fMean, 2);
    }
    $fVariance /= ( $bSample ? count($aValues) - 1 : count($aValues) );
    return (float) sqrt($fVariance);
}

?>
