<?PHP

require "config.php"; 

if (USERNAME!="USERNAME") {
	$username=USERNAME;
} else {
	$username="-";
};#if

$text.="

<TABLE class=' ' width='100%' border=0>
	<TR width='100%'>
		<TD class='suptitle ' width='50px'>
			$NGS_name
		</TD>
		<TD width='1px'></TD>
		<TD width='*' class=' subtitle'>
			<span class='' title='$NGS_description_title'><B>$NGS_description</B></span> [$NGS_release]
			<span style=''><BR>$NGS_comment<BR></span>
		</TD>
		<TD width='1px'></TD>
		<TD width='80px' class=' subtitle' style='text-align:center;'>
			<span style=''><B>".$username."</B><BR>".date("d-m-Y")."<BR>".date("H:i:s")."</span>
		</TD>
	</TR>
</TABLE>
<BR>
";
#-$NGS_release_date
?>


