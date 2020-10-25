<?php
/*
 *
 * Released under the terms and conditions of the
 * GNU General Public License (http://www.gnu.org/licenses/gpl.txt)
 *
*/

#
# Versioning
##############

$release="0.9.1";
$date="20131025";
$authors="Antony Le Béchec";

#
# HEADER
##########

if(!defined('e107_INIT'))
{
	require_once('../../class2.php');
}
$e107 = e107::getInstance();
$tp = e107::getParser();
$sql = e107::getDb();

#require_once(HEADERF);
require("functions.inc.php");
require "connect.php"; 
require "config.php"; 



$mime=($_REQUEST["mime"]!="")?$_REQUEST["mime"]:"image/png";
$type=($_REQUEST["type"]!="")?$_REQUEST["type"]:"snapshot";
$id=($_REQUEST["id"]!="")?$_REQUEST["id"]:"1";

$query="SELECT snapshot AS blob_image
	FROM $type
	WHERE id=$id
	"; 
#echo "$query";
$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
$summary= "<table class=' '>\n"; 
while ($row = mysql_fetch_assoc($result)) {
	$blob_image=$row['blob_image'];
	$summary.= "<tr>"
		."<td>blob_image</td><td>$blob_image</td>"
	."</tr>"; 
};#while
$summary.= "</table>\n"; 
#echo $summary;

#$filetmp="/media/IRCV2/NGSEnv/www/e107_plugins/NGSV2_plugin/tmp/filetmp_snapshot.png";
#file_put_contents($filetmp,$blob_image);
#echo $filetmp;

header("Content-Type:$mime");
echo $blob_image;
#echo file_get_contents($filetmp); #$blob_image;

?>
